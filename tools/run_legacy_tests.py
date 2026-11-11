#!/usr/bin/env python3
"""Run the DISCUS-suite macro tests unattended and report pass/fail.

Each test listed in tests/manifest.toml runs in its own fresh copy of its
top-level test folder, in a separate discus_suite process. The driver macro
sets LTEST_BATCH = TRUE, so tests/check.mac appends "PASS <name>" or
"FAIL <name>" to a result file instead of stopping. The driver writes "END"
as its last line; a missing END means the macro run aborted (for example on
a command error).

Statuses
  PASS     all expected checks passed and the run reached its end
  FAIL     at least one check wrote FAIL
  ABORTED  the run stopped before its end, or fewer checks than expected
  TIMEOUT  the run was killed after the test's timeout
  RAN      check = "none": the macro finished, but nothing was verified
  SKIPPED  enabled = false in the manifest, or filtered out by --tags

Exit code: 0 if no test is FAIL, ABORTED or TIMEOUT, else 1. 2 = usage error.

Only the Python standard library is used (Python >= 3.11 for tomllib).

Examples
  tools/run_legacy_tests.py --exe ../DiffuseBuild/suite/prog/discus_suite
  tools/run_legacy_tests.py --tags fast -j 4 --junit results.xml
  tools/run_legacy_tests.py --id '04.*' --keep
  tools/run_legacy_tests.py --list
  tools/run_legacy_tests.py --check-test-all
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import fnmatch
import json
import os
import queue
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import tomllib
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TESTS = REPO / "tests"
DEFAULT_MANIFEST = TESTS / "manifest.toml"

PASS, FAIL, ABORTED, TIMEOUT, RAN, SKIPPED = (
    "PASS", "FAIL", "ABORTED", "TIMEOUT", "RAN", "SKIPPED")
BAD = {FAIL, ABORTED, TIMEOUT}

ANSI = re.compile(r"\x1b\[[0-9;]*m")
ERROR_LINE = re.compile(r"\*\*\*[A-Z ]{3,5}\*\*\*")
XVFB_FAILED = "Xvfb failed to start"
XVFB_ATTEMPTS = 3


# --------------------------------------------------------------------- manifest
@dataclass
class Test:
    id: str
    dir: str
    macro: str
    args: str = ""
    check: str = "macro"            # macro | internal | none
    checks: int = 1
    tags: list[str] = field(default_factory=list)
    timeout: float = 600.0
    in_test_all: bool = True
    enabled: bool = True
    skip_reason: str = ""


def load_manifest(path: Path) -> list[Test]:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    tests, seen = [], set()
    for raw in data.get("test", []):
        unknown = set(raw) - set(Test.__dataclass_fields__)
        if unknown:
            raise ValueError(f"{raw.get('id')}: unknown field(s) {sorted(unknown)}")
        t = Test(**raw)
        if t.id in seen:
            raise ValueError(f"duplicate test id {t.id}")
        if t.check not in ("macro", "internal", "none"):
            raise ValueError(f"{t.id}: check must be macro, internal or none")
        if not (TESTS / t.dir / t.macro).is_file():
            raise ValueError(f"{t.id}: {t.dir}/{t.macro} not found")
        seen.add(t.id)
        tests.append(t)
    return tests


# ------------------------------------------------------- test_all.mac cross-check
def parse_test_all(path: Path) -> set[tuple[str, str]]:
    """(dir, macro) pairs that test_all.mac runs and then checks.

    Mirrors the structure of test_all.mac: 'cd <dir>' lines, '@<macro> <args>'
    followed by '@../check.mac', and macros that do their own checks
    (all.mac). Commented lines are ignored.
    """
    found: set[tuple[str, str]] = set()
    cwd = ""
    pending = None
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        m = re.match(r"cd\s+(\S+)", s)
        if m:
            cwd = os.path.normpath(os.path.join(cwd, m.group(1)))
            continue
        m = re.match(r"@(\S+)", s)
        if not m:
            continue
        mac = m.group(1)
        if "check.mac" in mac:
            if pending:
                found.add(pending)
                pending = None
            continue
        if mac in ("datum.mac", "kuplot_standard.mac"):
            continue
        if not mac.endswith(".mac"):
            mac += ".mac"
        pending = (cwd, mac)
        if mac == "all.mac":                      # checks itself
            found.add(pending)
            pending = None
        elif mac == "pdf.mac":                    # runs without a check
            found.add(pending)
            pending = None
    return found


def check_test_all(tests: list[Test]) -> int:
    in_mac = parse_test_all(TESTS / "test_all.mac")
    in_man = {(t.dir, t.macro) for t in tests if t.in_test_all}
    only_mac = sorted(in_mac - in_man)
    only_man = sorted(in_man - in_mac)
    for d, m in only_mac:
        print(f"only in test_all.mac : {d}/{m}")
    for d, m in only_man:
        print(f"only in manifest     : {d}/{m}")
    if not (only_mac or only_man):
        print(f"test_all.mac and manifest agree ({len(in_man)} macros)")
        return 0
    return 1


# ------------------------------------------------------------------------ running
@dataclass
class Result:
    id: str
    status: str
    seconds: float = 0.0
    checks: list[str] = field(default_factory=list)
    message: str = ""
    log: str = ""


def find_exe(given: str | None) -> Path | None:
    for cand in (given, os.environ.get("DISCUS_SUITE"), shutil.which("discus_suite")):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return Path(cand).resolve()
    return None


def macro_string(s: str) -> str:
    """Quote a value for a macro string assignment."""
    if '"' in s:
        raise ValueError(f"path contains a double quote: {s}")
    return s


def write_driver(t: Test, run_dir: Path, result_file: Path, root: Path) -> Path:
    args = f" {t.args}" if t.args else ""
    lines = [
        "# generated by tools/run_legacy_tests.py",
        "variable integer,   LTEST_ALL",
        "variable integer,   LTEST_MAC",
        "variable integer,   LTEST_BATCH",
        "variable character, LTEST_FILE",
        "LTEST_ALL   = TRUE",
        "LTEST_MAC   = FALSE",
        "LTEST_BATCH = TRUE",
        f'LTEST_FILE  = "%c", "{macro_string(str(result_file))}"',
        "set wait, off",
        f"@{t.macro}{args}",
    ]
    if t.check == "macro":
        lines.append(f"@{macro_string(str(root / 'check.mac'))} {t.macro}")
    lines += [
        'fopen 1, "%c", LTEST_FILE, append',
        'fput  1, "END"',
        "fclose 1",
        "exit",
    ]
    path = run_dir / "zz_run_legacy_test.mac"
    path.write_text("\n".join(lines) + "\n")
    return path


def copy_test_tree(t: Test, dest: Path) -> None:
    """Copy tests/<top-level folder> and the files directly in tests/."""
    dest.mkdir(parents=True)
    for item in TESTS.iterdir():
        if item.is_file():
            shutil.copy2(item, dest / item.name)
    top = Path(t.dir).parts[0]
    shutil.copytree(TESTS / top, dest / top, symlinks=True)


def kill_group(proc: subprocess.Popen) -> None:
    """Terminate a whole process group; SIGKILL whatever is left after 5 s."""
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    time.sleep(0.5)                     # let xvfb-run clean up its Xvfb
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()


def first_error(log: str) -> str:
    for line in log.splitlines():
        line = ANSI.sub("", line).strip()
        if ERROR_LINE.search(line) or XVFB_FAILED in line:
            return line
    return ""


def run_one(t: Test, exe: Path, work: Path, displays: queue.Queue, opts) -> Result:
    root = work / t.id
    if root.exists():
        shutil.rmtree(root)
    copy_test_tree(t, root)
    run_dir = root / t.dir
    result_file = root / "zz_results.txt"
    driver = write_driver(t, run_dir, result_file, root)
    log_path = work / f"{t.id}.log"

    cmd = [str(exe), "-macro", driver.name]
    if opts.xvfb:
        cmd = ["xvfb-run", "-n", "DISPLAY", "-s",
               "-screen 0 1280x1024x24"] + cmd
    env = dict(os.environ)
    if "PGPLOT_DIR" not in env and Path("/usr/lib/pgplot5").is_dir():
        env["PGPLOT_DIR"] = "/usr/lib/pgplot5"
    if opts.omp_threads:
        env["OMP_NUM_THREADS"] = str(opts.omp_threads)
    # Temporary files of discus_suite and xvfb-run stay inside the work dir.
    env["TMPDIR"] = str(root / "tmp")
    (root / "tmp").mkdir()

    display = displays.get()            # one X display per running test
    t0 = time.time()
    status = None
    try:
        for attempt in range(XVFB_ATTEMPTS):
            if opts.xvfb:
                # First try the display reserved for this job; if Xvfb cannot
                # start there (display taken by another program or runner),
                # try a random high display number.
                cmd[2] = str(display if attempt == 0 else random.randint(1000, 30000))
            with open(log_path, "w") as log:
                # Own process group, so a timeout also stops discus_suite and
                # Xvfb, not only the xvfb-run wrapper.
                proc = subprocess.Popen(cmd, cwd=run_dir, stdin=subprocess.DEVNULL,
                                        stdout=log, stderr=subprocess.STDOUT, env=env,
                                        start_new_session=True)
                try:
                    proc.wait(timeout=t.timeout * opts.timeout_scale)
                except subprocess.TimeoutExpired:
                    status = TIMEOUT
                    kill_group(proc)
            if not (opts.xvfb and XVFB_FAILED in log_path.read_text(errors="replace")):
                break
    finally:
        displays.put(display)
    seconds = round(time.time() - t0, 1)

    log_text = log_path.read_text(errors="replace")
    raw = result_file.read_text().splitlines() if result_file.exists() else []
    finished = "END" in [r.strip() for r in raw]
    checks = [r.strip() for r in raw if r.strip() and r.strip() != "END"]
    n_pass = sum(c.startswith("PASS") for c in checks)
    n_fail = sum(c.startswith("FAIL") for c in checks)

    msg = ""
    if status == TIMEOUT:
        msg = f"killed after {t.timeout * opts.timeout_scale:.0f} s"
    elif n_fail:
        status = FAIL
        msg = "; ".join(c for c in checks if c.startswith("FAIL"))
    elif not finished:
        status = ABORTED
        msg = first_error(log_text) or "macro run ended before its last line"
    elif t.check == "none":
        status = RAN
        msg = "no check defined; the macro finished"
    elif n_pass < t.checks:
        status = ABORTED
        msg = f"expected {t.checks} check(s), got {n_pass}"
    else:
        status = PASS

    if not opts.keep and status not in BAD:
        shutil.rmtree(root, ignore_errors=True)
    return Result(t.id, status, seconds, checks, msg, str(log_path))


# ------------------------------------------------------------------------ reports
def write_junit(results: list[Result], tests: dict[str, Test], path: Path) -> None:
    suite = ET.Element("testsuite", name="discus-legacy-macros",
                       tests=str(len(results)),
                       failures=str(sum(r.status == FAIL for r in results)),
                       errors=str(sum(r.status in (ABORTED, TIMEOUT) for r in results)),
                       skipped=str(sum(r.status == SKIPPED for r in results)),
                       time=f"{sum(r.seconds for r in results):.1f}")
    for r in results:
        t = tests[r.id]
        case = ET.SubElement(suite, "testcase", name=r.id,
                             classname=t.dir.replace("/", "."), time=f"{r.seconds:.1f}")
        if r.status == FAIL:
            ET.SubElement(case, "failure", message=r.message)
        elif r.status in (ABORTED, TIMEOUT):
            ET.SubElement(case, "error", message=f"{r.status}: {r.message}")
        elif r.status == SKIPPED:
            ET.SubElement(case, "skipped", message=r.message)
        if r.status in BAD and r.log and Path(r.log).exists():
            tail = Path(r.log).read_text(errors="replace").splitlines()[-200:]
            ET.SubElement(case, "system-out").text = ANSI.sub("", "\n".join(tail))
    ET.ElementTree(suite).write(path, encoding="utf-8", xml_declaration=True)


def print_summary(results: list[Result]) -> None:
    width = max(len(r.id) for r in results)
    for r in results:
        extra = f"  {r.message}" if r.message else ""
        print(f"{r.id:<{width}}  {r.status:<8} {r.seconds:7.1f} s{extra}")
    counts = {s: sum(r.status == s for r in results)
              for s in (PASS, FAIL, ABORTED, TIMEOUT, RAN, SKIPPED)}
    print("-" * (width + 20))
    print("  ".join(f"{k}: {v}" for k, v in counts.items() if v)
          + f"   total {sum(r.seconds for r in results):.0f} s")


# --------------------------------------------------------------------------- main
def select(tests: list[Test], opts) -> tuple[list[Test], list[Result]]:
    run, skipped = [], []
    want = set(opts.tags.split(",")) if opts.tags else None
    drop = set(opts.exclude_tags.split(",")) if opts.exclude_tags else set()
    for t in tests:
        if opts.id and not any(fnmatch.fnmatch(t.id, p) for p in opts.id):
            continue
        if not t.enabled:
            skipped.append(Result(t.id, SKIPPED, message=t.skip_reason))
        elif want is not None and not want & set(t.tags):
            continue
        elif drop & set(t.tags):
            skipped.append(Result(t.id, SKIPPED, message="excluded by tag"))
        else:
            run.append(t)
    return run, skipped


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", help="discus_suite binary (default: $DISCUS_SUITE, then PATH)")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--id", action="append",
                   help="run only tests whose id matches this glob (repeatable)")
    p.add_argument("--tags", help="run only tests with at least one of these tags (comma list)")
    p.add_argument("--exclude-tags", default="mpi",
                   help="skip tests with any of these tags (default: mpi)")
    p.add_argument("-j", "--jobs", type=int, default=1, help="tests to run in parallel")
    p.add_argument("--omp-threads", type=int,
                   help="OMP_NUM_THREADS per test (default: CPUs / jobs when jobs > 1)")
    p.add_argument("--workdir", type=Path, help="where test copies and logs go (default: temp)")
    p.add_argument("--keep", action="store_true", help="keep copies of passing tests too")
    p.add_argument("--junit", type=Path, help="write JUnit XML here")
    p.add_argument("--json", type=Path, help="write a JSON summary here")
    p.add_argument("--xvfb", choices=("auto", "on", "off"), default="auto",
                   help="run under xvfb-run (auto: when DISPLAY is unset)")
    p.add_argument("--display-base", type=int,
                   help="first X display number used with xvfb "
                        "(default: 100 + a number derived from the process ID)")
    p.add_argument("--timeout-scale", type=float, default=1.0)
    p.add_argument("--list", action="store_true", help="list the selected tests and exit")
    p.add_argument("--check-test-all", action="store_true",
                   help="compare the manifest with test_all.mac and exit")
    opts = p.parse_args(argv)

    try:
        tests = load_manifest(opts.manifest)
    except (OSError, ValueError, tomllib.TOMLDecodeError) as exc:
        print(f"manifest error: {exc}", file=sys.stderr)
        return 2

    if opts.check_test_all:
        return check_test_all(tests)

    to_run, skipped = select(tests, opts)
    if opts.list:
        for t in to_run:
            print(f"{t.id:<20} {t.dir}/{t.macro} {t.args}  [{', '.join(t.tags)}]")
        for r in skipped:
            print(f"{r.id:<20} SKIPPED  {r.message}")
        return 0

    exe = find_exe(opts.exe)
    if exe is None:
        print("discus_suite not found: use --exe, $DISCUS_SUITE or PATH", file=sys.stderr)
        return 2

    if opts.xvfb == "auto":
        opts.xvfb = not os.environ.get("DISPLAY") and shutil.which("xvfb-run") is not None
    else:
        opts.xvfb = opts.xvfb == "on"
    if opts.display_base is None:       # keep concurrent runners apart
        opts.display_base = 100 + (os.getpid() % 100) * 8
    if opts.jobs > 1 and not opts.omp_threads:
        opts.omp_threads = max(1, (os.cpu_count() or 1) // opts.jobs)

    work = opts.workdir or Path(tempfile.mkdtemp(prefix="discus_tests_"))
    work = work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    print(f"discus_suite : {exe}\nwork dir     : {work}\n"
          f"tests        : {len(to_run)} to run, {len(skipped)} skipped\n", flush=True)

    by_id = {t.id: t for t in tests}
    results: dict[str, Result] = {}
    jobs = max(1, opts.jobs)
    displays: queue.Queue = queue.Queue()
    for n in range(jobs):
        displays.put(opts.display_base + n)
    # Tests tagged "exclusive" set their own OpenMP thread count (for example
    # "set parallel, useomp:use, nthreads:all") and slow everything down badly
    # when they share the CPUs. They run one at a time after the others.
    shared = [t for t in to_run if "exclusive" not in t.tags]
    alone = [t for t in to_run if "exclusive" in t.tags]
    for batch, width in ((shared, jobs), (alone, 1)):
        with cf.ThreadPoolExecutor(max_workers=width) as pool:
            futures = [pool.submit(run_one, t, exe, work, displays, opts) for t in batch]
            for fut in cf.as_completed(futures):
                r = fut.result()
                results[r.id] = r
                print(f"  {r.status:<8} {r.id} ({r.seconds} s)", flush=True)

    ordered = [results[t.id] for t in to_run] + skipped
    print()
    print_summary(ordered)
    if opts.junit:
        write_junit(ordered, by_id, opts.junit)
    if opts.json:
        opts.json.write_text(json.dumps([asdict(r) for r in ordered], indent=2))
    bad = [r for r in ordered if r.status in BAD]
    if bad:
        print(f"\nLogs of failed tests are in {work}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

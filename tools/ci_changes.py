#!/usr/bin/env python3
"""Decide which CI jobs a push needs (used by .github/workflows/tests.yml).

Every changed file falls into one group:

  text       plan and notes: migration/, *.md outside docs/, License,
             release notes, .gitignore             -> no job
  docs       docs/                                 -> docs build
  python     the Python package and its tests      -> Python tests, docs
  packaging  pyproject.toml, environment.yml, CMake files, the f2py
             wrapper, the workflow and this script -> everything
  fortran    everything else (Fortran and C sources, legacy macro tests,
             their runner and reference data)      -> build, all tests

Unknown files count as fortran, so a new kind of file never skips tests.

What is compared:

  push to a claude/** branch    the whole branch, against its merge base
                                with claude-python-app (so later pushes
                                keep running what earlier ones needed)
  push to claude-python-app     the push itself (merges: the merged PR)
  pull_request (opened)         if the branch already contains the latest
                                claude-python-app, nothing runs: the push
                                runs tested exactly this code. Otherwise
                                the PR's changes, filtered as for a push,
                                tested on the merge result (catches clashes
                                with what was merged since the branch began)
  schedule, workflow_dispatch   everything runs

Outputs (GITHUB_OUTPUT): build, legacy, conda, docs ("true"/"false") and
package_os, the JSON list of Ubuntu versions for the python-package job.

    python3 tools/ci_changes.py --files a.f90 docs/x.md    # try it locally
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import subprocess
import sys

BASE_BRANCH = "claude-python-app"
ALL_OS = ["ubuntu-24.04", "ubuntu-22.04"]
WSL2_OS = "ubuntu-22.04"          # the main development environment
DOCS_OS = "ubuntu-24.04"          # builds (and publishes) the docs

# First match wins.
GROUPS = [
    ("packaging", ["pyproject.toml", "environment.yml", "CMakeLists.txt",
                   "*/CMakeLists.txt", "cmake/*", "*.cmake",
                   "python/suite_python.f90", ".github/*",
                   "tools/ci_changes.py", "tests/requirements.txt"]),
    ("docs", ["docs/*"]),
    ("text", ["migration/*", "*.md", "License.txt", "release_notes.txt",
              ".gitignore", ".travis.yml"]),
    ("python", ["python/*", "tests/python/*", "tests/golden/*"]),
]


def group(path: str) -> str:
    for name, patterns in GROUPS:
        if any(fnmatch.fnmatch(path, p) for p in patterns):
            return name
    return "fortran"


def decide(files: list[str] | None) -> dict[str, str]:
    """files=None means: run everything."""
    if files is None:
        groups = {"packaging"}
    else:
        groups = {group(f) for f in files}
    full = "packaging" in groups
    code = full or bool(groups & {"python", "fortran"})
    docs = full or bool(groups & {"python", "docs"})
    if full:
        os_list = list(ALL_OS)
    else:
        os_list = ([DOCS_OS] if docs else []) + ([WSL2_OS] if code else [])
    return {
        "build": str(code).lower(),                     # linux job: pytest
        "legacy": str(full or "fortran" in groups).lower(),
        "conda": str(full).lower(),
        "docs": str(docs).lower(),
        "package_os": json.dumps(os_list),
    }


def git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True,
                          text=True).stdout


def is_ancestor(older: str, newer: str) -> bool:
    return subprocess.run(["git", "merge-base", "--is-ancestor", older, newer],
                          capture_output=True).returncode == 0


def changed_files() -> tuple[list[str] | None, str]:
    """The changed files (None: run everything) and a note on why."""
    event = os.environ.get("GITHUB_EVENT_NAME", "")
    ref = os.environ.get("GITHUB_REF_NAME", "")
    try:
        if event == "pull_request":
            base, head = os.environ.get("PR_BASE_SHA", ""), os.environ.get("PR_HEAD_SHA", "")
            if not base or not head:
                return None, "pull request without base/head commits"
            if is_ancestor(base, head):
                return [], (f"The branch already contains {BASE_BRANCH} ({base[:7]}): its "
                            "push runs tested this exact code, so nothing runs here.")
            files = [f for f in git("diff", "--name-only", f"{base}...{head}").splitlines() if f]
            return files, (f"{BASE_BRANCH} moved on since the branch began: testing the "
                           "merge result, for the jobs the PR's changes need.")
        if event != "push":
            return None, f"{event or 'local'}: everything runs"
        if ref == BASE_BRANCH:
            base = os.environ.get("BEFORE", "")
            if not base or set(base) == {"0"}:
                return None, "new branch: everything runs"
        else:
            base = git("merge-base", f"origin/{BASE_BRANCH}", "HEAD").strip()
        return [f for f in git("diff", "--name-only", base, "HEAD").splitlines() if f], ""
    except subprocess.CalledProcessError as exc:   # unknown base: run everything
        return None, f"git failed ({exc.stderr.strip()}); everything runs"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--files", nargs="*", help="classify these paths instead of using git")
    args = ap.parse_args()
    files, note = (args.files, "") if args.files is not None else changed_files()
    out = decide(files)

    report = [note] if note else []
    if files is None:
        report.append("Running everything.")
    else:
        report.append(f"{len(files)} changed file(s):")
        for f in files[:200]:
            report.append(f"  {group(f):9}  {f}")
    report += [f"{k} = {v}" for k, v in out.items()]
    print("\n".join(report))

    for var, text in (("GITHUB_OUTPUT", "".join(f"{k}={v}\n" for k, v in out.items())),
                      ("GITHUB_STEP_SUMMARY", "```\n" + "\n".join(report) + "\n```\n")):
        if os.environ.get(var):
            with open(os.environ[var], "a") as fh:
                fh.write(text)


if __name__ == "__main__":
    main()

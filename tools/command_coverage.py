#!/usr/bin/env python3
"""Which documented commands do the macro tests actually run?

Reads the documented commands from the help files (level 2 = command,
level 3 = subcommand or option) and the commands echoed in the logs of a
tools/run_legacy_tests.py run, and reports what is and is not exercised.

    python3 tools/run_legacy_tests.py --exe .../discus_suite --workdir /tmp/run
    python3 tools/command_coverage.py --logs /tmp/run --markdown coverage.md

Limits (it is a guide, not an exact measure):
  * Macro echo also prints lines of if-branches that are not taken, so a
    command can count as "reached" without being executed.
  * Commands are matched by unique prefix, as the command language allows
    abbreviations. Inside do-loops in a sub-menu the log shows the outer
    prompt; such lines are attributed to the most recently entered menu.
  * "Reached" says nothing about whether the result is checked.

mixscat (deprecated) is not included.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
HELP = {
    "suite": "suite/prog/appl_suite.hlp",
    "discus": "discus/prog/appl_dis.hlp",
    "kuplot": "kuplot/prog/appl_kup.hlp",
    "diffev": "diffev/prog/appl_dif.hlp",
    "refine": "refine/prog/appl_ref.hlp",
    "experi": "experi/prog/appl_exp.hlp",
    "lib": "lib_f90/lib_f90.hlp",           # command language, all programs
}
# level-2 headings that are topics, not commands
TOPICS = {"summary", "news", "manual", "menu", "unit_cell", "functions", "variables",
          "bugs", "errors", "space_groups", "wavelength", "atom_names", "synopsis",
          "description", "example", "internal", "help", "command_language", "options",
          "expressions", "filenames", "commands", "2dm", "hpackage"}
SUB_TOPICS = {"commands", "summary", "news", "examples", "example", "general", "options",
              "appl", "pgplot", "special"}
ALIAS = {"else": "if", "elseif": "if", "endif": "if", "enddo": "do"}
PROMPT = re.compile(r"^\s*([a-z_]+(?:/[a-z_]+)*)\s?>\s?(.*)$")
ASSIGN = re.compile(r"^[A-Za-z_]\w*(\[[^\]]*\])?\s*=")
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def expand(name: str) -> list[str]:
    name = name.lower()
    return {"achx/y/z": ["achx", "achy", "achz"], "tit1/2": ["tit1", "tit2"]}.get(name, [name])


def read_help() -> dict[str, dict[str, set[str]]]:
    doc: dict[str, dict[str, set[str]]] = {}
    for prog, rel in HELP.items():
        cmds: dict[str, set[str]] = {}
        cur = None
        for line in (REPO / rel).read_text(errors="replace").splitlines():
            m = re.match(r"^([234])\s+(\S+)", line)
            if not m:
                continue
            lvl, name = int(m.group(1)), m.group(2)
            if prog == "lib":
                if lvl == 3 and name.lower() not in TOPICS and name not in "#@=" \
                        and not name[0].isdigit():
                    cmds[name.lower()] = set()
                continue
            if lvl == 2:
                cur = None
                if name.lower() not in TOPICS:
                    for n in expand(name):
                        cmds[n] = set()
                    cur = expand(name)[0]
            elif lvl == 3 and cur and not name[0].isdigit() \
                    and name.lower() not in SUB_TOPICS:
                cmds[cur].add(name.lower())
        doc[prog] = cmds
    return doc


def match(token: str, names) -> str | None:
    """Exact name, or the single name that the token abbreviates."""
    if token in names:
        return token
    cands = [n for n in names if n.startswith(token) and len(token) >= 2]
    return cands[0] if len(cands) == 1 else None


def analyse(doc, logdir: Path):
    used_cmd = defaultdict(set)                       # (prog, cmd) -> tests
    used_sub = defaultdict(lambda: defaultdict(set))  # (prog, cmd) -> sub -> tests
    setparams = defaultdict(lambda: defaultdict(set))
    unknown = defaultdict(set)
    level3 = defaultdict(list)
    for prog, cmds in doc.items():
        for c, subs in cmds.items():
            level3[prog] += [(c, s) for s in subs]

    def lib(tok):
        return match(ALIAS.get(tok, tok), doc["lib"])

    def sub(prog, cmd, tok, test):
        s = match(tok, doc[prog][cmd]) if doc[prog][cmd] else None
        if s:
            used_sub[(prog, cmd)][s].add(test)
        return s

    for log in sorted(logdir.glob("*.log")):
        test, last_menu = log.stem, {}
        for raw in log.read_text(errors="replace").splitlines():
            m = PROMPT.match(ANSI.sub("", raw))
            if not m:
                continue
            prompt, line = m.group(1), m.group(2).strip()
            if not line or line[0] in "#!@" or ASSIGN.match(line):
                continue
            words = [w.lower() for w in re.split(r"[\s,]+", line) if w]
            tok = words[0]
            if not re.match(r"^[a-z][a-z0-9_]*$", tok):
                continue
            parts = prompt.split("/")
            prog = parts[0]
            if prog not in doc:
                continue
            if len(parts) > 1:                                 # inside a menu
                menu = match(parts[1], doc[prog])
                if not menu:
                    unknown[f"{prompt}> {tok}"].add(test)
                    continue
                last_menu[prog] = menu
                used_cmd[(prog, menu)].add(test)
                if not sub(prog, menu, tok, test):
                    if (lh := lib(tok)):
                        used_cmd[("lib", lh)].add(test)
                    else:
                        unknown[f"{prompt}> {tok}"].add(test)
                if tok == "set" and len(words) > 1:
                    setparams[f"{prog}/{menu}"][words[1].split(":")[0]].add(test)
                continue
            hit = match(tok, doc[prog])                         # top level
            if hit:
                used_cmd[(prog, hit)].add(test)
                if len(words) > 1:
                    sub(prog, hit, words[1], test)              # e.g. "load xy"
                continue
            if (lh := lib(tok)):
                used_cmd[("lib", lh)].add(test)
                continue
            lm = last_menu.get(prog)
            if lm and sub(prog, lm, tok, test):                 # do-loop in a menu
                if tok == "set" and len(words) > 1:
                    setparams[f"{prog}/{lm}"][words[1].split(":")[0]].add(test)
                continue
            cands = {(c, s) for c, s in level3[prog]
                     if s == tok or (s.startswith(tok) and len(tok) >= 3)}
            if len({c for c, _ in cands}) == 1:                 # e.g. kuplot "nfra"
                c, s = sorted(cands)[0]
                used_cmd[(prog, c)].add(test)
                used_sub[(prog, c)][s].add(test)
                continue
            unknown[f"{prog}> {tok}"].add(test)

    result = {}
    for prog, cmds in doc.items():
        rows = []
        for c in sorted(cmds):
            subs = sorted(cmds[c])
            rows.append({
                "cmd": c,
                "tests": sorted(used_cmd.get((prog, c), [])),
                "subs": subs,
                "subs_used": [s for s in subs if used_sub[(prog, c)].get(s)],
                "subs_unused": [s for s in subs if not used_sub[(prog, c)].get(s)],
            })
        result[prog] = rows
    return {
        "commands": result,
        "set_params": {k: {p: sorted(v) for p, v in d.items()} for k, d in setparams.items()},
        "unmatched": {k: sorted(v) for k, v in unknown.items()},
    }


def summary_rows(data):
    for prog, rows in data["commands"].items():
        n, c = len(rows), sum(bool(r["tests"]) for r in rows)
        ns = sum(len(r["subs"]) for r in rows)
        cs = sum(len(r["subs_used"]) for r in rows)
        yield prog, c, n, cs, ns


def markdown(data) -> str:
    out = ["| Program | Commands reached | Sub-commands / options reached |",
           "|---|---|---|"]
    for prog, c, n, cs, ns in summary_rows(data):
        sub = f"{cs} / {ns}" if ns else "–"
        out.append(f"| {prog} | {c} / {n} | {sub} |")
    for prog, rows in data["commands"].items():
        out += ["", f"#### {prog}", ""]
        miss = [r["cmd"] for r in rows if not r["tests"]]
        out.append("**Not reached:** " + (", ".join(f"`{m}`" for m in miss) or "none"))
        reached = [r for r in rows if r["tests"]]
        if reached:
            out += ["", "| Command | Tests | Sub-commands reached | Sub-commands not reached |",
                    "|---|---:|---|---|"]
            for r in reached:
                used = f"{len(r['subs_used'])} / {len(r['subs'])}" if r["subs"] else "–"
                unused = ", ".join(r["subs_unused"]) or "–"
                out.append(f"| `{r['cmd']}` | {len(r['tests'])} | {used} | {unused} |")
    if data["set_params"]:
        out += ["", "#### `set` parameters used inside menus", "",
                "| Menu | Parameters (number of tests) |", "|---|---|"]
        for k, v in sorted(data["set_params"].items()):
            out.append(f"| {k} | " + ", ".join(f"{p} ({len(t)})" for p, t in sorted(v.items()))
                       + " |")
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--logs", type=Path, required=True,
                   help="work dir of a run_legacy_tests.py run (contains <test>.log)")
    p.add_argument("--json", type=Path, help="write full results as JSON")
    p.add_argument("--markdown", type=Path, help="write tables as Markdown")
    opts = p.parse_args(argv)
    if not any(opts.logs.glob("*.log")):
        print(f"no *.log files in {opts.logs}", file=sys.stderr)
        return 2
    data = analyse(read_help(), opts.logs)
    for prog, c, n, cs, ns in summary_rows(data):
        print(f"{prog:7s} commands {c:3d}/{n:3d}   sub-commands/options {cs:3d}/{ns:3d}")
    print(f"unmatched command words: {len(data['unmatched'])}")
    if opts.json:
        opts.json.write_text(json.dumps(data, indent=1))
    if opts.markdown:
        opts.markdown.write_text(markdown(data))
    return 0


if __name__ == "__main__":
    sys.exit(main())

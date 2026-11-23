#!/usr/bin/env python3
"""How much of the Fortran code is physics, and how much is the rest?

Sorts every Fortran file of the suite (code lines, without comments and
blank lines) into what the Python migration replaces and what would have to
be replicated if the Fortran were rewritten (migration/WP2_Binding_Decision.md,
option E):

  A  command language, menus, program frame, macro and expression language
  B  plotting (KUPLOT graphics, structure plots)
  C  fitting and optimisation (DIFFEV, REFINE, KUPLOT fit, weights)
  D  file formats and export
  E  general numerics (FFT, splines, matrices, random numbers, helpers)
  F  system (operating system, MPI, sockets, timers)
  G  physics core
  H  data tables (scattering factors, profile functions)

Most files fall into one group by folder or name. DISCUS files that hold
both a menu and the physics behind it are split by subroutine: a routine
counts as menu (A) if its name says so (*menu*, do_set*, *_show, *setup*,
*_kdo, *_loop, *reset*, ...) or if more than 4 % of its lines call the
command parser (str_comp, get_params, get_optional, ber_params, cpara).

This is an estimate, good to perhaps +-10 % for the physics core: some
routines mix parsing and calculation, and error messages and checks inside
physics routines count as physics.

    python3 tools/fortran_inventory.py              # summary
    python3 tools/fortran_inventory.py --files 20   # plus the largest files per group
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

GROUPS = {
    "A": "command language, menus, program frame",
    "B": "plotting",
    "C": "fitting and optimisation",
    "D": "file formats and export",
    "E": "general numerics",
    "F": "system (OS, MPI, sockets, timers)",
    "G": "physics core",
    "H": "data tables",
}

PHYSICS_AREAS = [
    ("Fourier, powder, PDF, Patterson, 3D-PDF",
     r"four|powder|pdf|patters|debye|diffuse_mod|intens|qval|bragg|phases|discamb|exp2pdf|exp3_pdf|"
     r"mrc_|discus_fit|estimate|mixscat|thermal|prep_anis|tensors|quad"),
    ("structure: reading, symmetry, space groups, generators",
     r"structur|spcgr|symm|generate|sym_add|wyckoff|unitcell|celltoindex|pointgrp|metric|trafo|transfrm|"
     r"trans_|inverse|update_cr|perioditize|discus_super|discus_reduce|discus_crystal|crystal_mod|class_|"
     r"atom_|get_iscat|check_|find_top|do_find|private|zone|insert|recipro|refine_mod|gen_add"),
    ("Monte Carlo (MMC, RMC)", r"mmc|rmc|mc_mod|extrmc"),
    ("chemistry, correlations, connectivity, properties",
     r"chem|conn_|prop_|atom_env|mole_env|bv_data"),
    ("molecules, domains, stacking, surfaces, waves, modify",
     r"molecule|mole|place_|domain|stack|surface|deco|demolec|waves|modify|shear|micro|external|guess|"
     r"storage|build_molecule|interpret|read_internal"),
]

CMD = re.compile(r"str_comp\s*\(|get_params\s*\(|get_optional\s*\(|ber_params\s*\(|cpara\s*\(", re.I)
START = re.compile(r"^\s*(?:recursive\s+|pure\s+|elemental\s+|module\s+)*"
                   r"(?:(?:integer|real|logical|character|complex|double\s+precision)[^!]*?\s+)?"
                   r"(subroutine|function)\s+(\w+)", re.I)
END = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
MENU_NAME = re.compile(r"menu|^do_set|_set$|set_sub|_show|^show|setup|_kdo|kdo_|_loop|upd_par|err_appl|"
                       r"blk_appl|branch|_exit|^exit|_reset|^rese|_help|^do_hel|_alloc|allocate|_init$|_default",
                       re.I)


def rx(*patterns: str) -> re.Pattern:
    return re.compile("|".join(patterns))


def code_lines(path: Path) -> list[str]:
    return [line for line in path.read_text(errors="replace").splitlines()
            if line.strip() and not line.strip().startswith("!")]


def split_menu_physics(path: Path) -> tuple[int, int]:
    """(menu lines, physics lines) of a file, routine by routine."""
    lines = code_lines(path)
    routines, current, depth = [], None, 0
    for line in lines:
        m = START.match(line)
        if m and not line.strip().lower().startswith("end"):
            if depth == 0:
                current = [m.group(2), 0, 0]
            depth += 1
        if current is not None:
            current[1] += 1
            if CMD.search(line):
                current[2] += 1
        if END.match(line):
            depth -= 1
            if depth == 0 and current is not None:
                routines.append(current)
                current = None
    outside = len(lines) - sum(r[1] for r in routines)   # module data: goes with the physics
    menu = sum(n for name, n, c in routines if MENU_NAME.search(name) or (n and c / n > 0.04))
    return menu, len(lines) - menu, outside


# --- rules by folder ---------------------------------------------------------------
DISCUS_D = rx(r"^discus_save", r"inter_readstru", r"read_internal", r"discus_export", r"_out_file", r"_s_.*_file",
              r"discus_pdf_file", r"discus_powder_file", r"powder_write", r"^output_", r"^graphic", r"discus_mrc",
              r"vtk_mod", r"discus_xplor", r"discus_nexus", r"hdf_write", r"nipl_header", r"^storage",
              r"trans_to_short", r"discus_exp2pdf_(load|write)", r"discus_exp3_pdf_load", r"powder_out_partial")
DISCUS_B = rx(r"^discus_plot")
DISCUS_A = rx(r"discus_prepare_refine", r"discus_upd_par", r"discus_set_sub", r"discus_setup", r"_blk_appl",
              r"_err_appl", r"discus_kdo", r"discus_loop", r"discus_branch", r"discus_exit", r"discus_reset_all",
              r"_reset\.f90", r"stack_rese", r"discus_show", r"discus_py", r"^discus\.f90", r"allocate_appl",
              r"discus_config_mod", r"discus_init")
DISCUS_E = rx(r"^utilities", r"discus_fft")
FRAME = r"kdo|loop|setup|upd_par|err_appl|blk_appl|branch|exit|reset|show|_py|config|set_sub"
LIB = [  # first match wins
    ("C", rx(r"wichtung", r"lib_weights", r"refine_control", r"refine_data", r"refine_params")),
    ("H", rx(r"element_data", r"gauss_lorentz", r"profile_tof")),
    ("F", rx(r"appl_unix", r"op_", r"operating", r"sys_", r"parallel", r"mpi", r"sockets", r"remote", r"send_get",
             r"times", r"timer", r"lib_f90_profile", r"exit_para", r"dummy_loop")),
    ("D", rx(r"lib_detect", r"hdf", r"nx_", r"unified", r"lib_load", r"lib_conver_shelx", r"mrc_data", r"lib_h5fortran",
             r"lib_forpython", r"data_struc", r"map_1dtofield")),
    ("E", rx(r"fast_fourier", r"lib_trans", r"lib_use_coor", r"lib_conv\.f90", r"fft", r"spline", r"gaussj",
             r"matrix", r"sorting", r"random", r"gamma", r"lib_lanczos", r"lib_math", r"math_sup", r"lib_functions",
             r"sine_table", r"precision_mod", r"arrays_mod", r"trig_degree", r"wink", r"lib_metric")),
    ("A", re.compile(r".")),   # everything else in lib_f90: command language, variables, macros, terminal ...
]


def classify() -> tuple[dict[str, int], dict[str, list[tuple[int, str]]]]:
    totals = {g: 0 for g in GROUPS}
    files = {g: [] for g in GROUPS}

    def add(group: str, n: int, name: str) -> None:
        if n:
            totals[group] += n
            files[group].append((n, name))

    for f in sorted((REPO / "discus/prog").glob("*.f90")):
        n, name = len(code_lines(f)), "discus/" + f.name
        if DISCUS_D.search(f.name):
            add("D", n, name)
        elif DISCUS_B.search(f.name):
            add("B", n, name)
        elif DISCUS_A.search(f.name):
            add("A", n, name)
        elif DISCUS_E.search(f.name):
            add("E", n, name)
        else:
            menu, physics, _ = split_menu_physics(f)
            add("A", menu, name + " (menu part)")
            add("G", physics, name)
    for f in sorted((REPO / "kuplot/prog").glob("*.f90")):
        n, name = len(code_lines(f)), "kuplot/" + f.name
        if re.search(r"^fit|kuplot_wichtung|kuplot_extrema|^math|^calc|^para", f.name):
            add("C", n, name)
        elif re.search(r"load|save|nexus|gsas|_h5|shelx|_rmc", f.name):
            add("D", n, name)
        elif re.search(FRAME + r"|kuplot\.f90|_top|words|global|toglobal|kuplot_mod|kuplot_low", f.name):
            add("A", n, name)
        else:
            add("B", n, name)
    for folder in ("diffev", "refine"):
        for f in sorted((REPO / folder / "prog").glob("*.f90")):
            n, name = len(code_lines(f)), f"{folder}/{f.name}"
            add("F" if "mpi" in f.name else "A" if re.search(FRAME + r"|^diffev\.f90", f.name) else "C", n, name)
    for folder in ("suite", "experi"):
        for f in sorted((REPO / folder / "prog").glob("*.f90")):
            add("C" if "execute_cost" in f.name else "A", len(code_lines(f)), f"{folder}/{f.name}")
    for f in sorted(list((REPO / "mixscat/prog").glob("*.f90")) + list((REPO / "mixscat/prog").glob("*.f"))):
        n, name = len(code_lines(f)), "mixscat/" + f.name
        group = ("A" if re.search(r"kdo|blk_appl|err_appl|exit|set|show|upd_par|config|^mixscat\.f90|lazy", f.name)
                 else "D" if re.search(r"load|save", f.name) else "G")
        add(group, n, name)
    for f in sorted((REPO / "lib_f90").glob("*.f90")):
        if "forpy" in f.name:          # bundled third-party module (Fortran calling Python)
            continue
        group = next(g for g, pattern in LIB if pattern.search(f.name))
        add(group, len(code_lines(f)), "lib_f90/" + f.name)
    return totals, files


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--files", type=int, default=0, metavar="N", help="also list the N largest files per group")
    args = ap.parse_args()
    totals, files = classify()
    grand = sum(totals.values())
    print(f"Fortran code lines (no comments, no blank lines): {grand}\n")
    for g, label in GROUPS.items():
        print(f"  {g}  {totals[g]:7d}  {100 * totals[g] / grand:5.1f} %  {label}")
    print("\nPhysics core (G) by area:")
    areas = {label: 0 for label, _ in PHYSICS_AREAS}
    for n, name in files["G"]:
        for label, pattern in PHYSICS_AREAS:
            if re.search(pattern, name):
                areas[label] += n
                break
        else:
            areas.setdefault("other", 0)
            areas["other"] = areas.get("other", 0) + n
    for label, n in areas.items():
        print(f"     {n:7d}  {label}")
    if args.files:
        for g, label in GROUPS.items():
            print(f"\n{g} {label}:")
            for n, name in sorted(files[g], reverse=True)[:args.files]:
                print(f"     {n:6d}  {name}")


if __name__ == "__main__":
    main()

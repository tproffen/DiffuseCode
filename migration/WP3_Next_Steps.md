# WP3 next steps: asymmetric unit from Python, and the option E pilot

Status (Nov 2026) and two pieces of planned work in detail: building a
crystal in DISCUS from an asymmetric unit and a space group given in Python,
and the option E pilot (a core kernel rebuilt outside Fortran), with its
benchmarks against the Fortran.

## Where WP3 stands

Through the C interface (`python/diffuse_capi.f90`, `python/diffuse_core.h`,
version 4; Python in `python/diffuse/_core.py` and `crystal.py`):

| | out of DISCUS | into DISCUS |
|---|---|---|
| cell, space group, atoms, types, B, occupancy | `Crystal.atoms()`, `atom_types()`, `space_group` | `dc.make_crystal()` |
| anisotropic displacement parameters | `Crystal.adps()` (per atom, symmetry-rotated) | `make_crystal(u_aniso=...)` |
| molecules | `Crystal.molecules()` | `make_crystal(molecules=...)` |
| ASE | `Crystal.to_ase()` (with names and tensors) | `dc.from_ase()` |

Every direction is checked against DISCUS reading the same data from a file:
identical arrays in round trips, and identical powder patterns, Fourier maps
and PDFs. `Crystal.save()` (`.stru`, `.h5`) and `Crystal.export()` (CIF,
SHELX, RMCProfile, VASP, Scatty) are in place, each tested by reading back
what it writes (`tests/python/test_crystal_save_export.py`).

### DISCUS bugs found so far (crystal I/O)

Kept as strict xfail tests (they turn red once fixed), the last one as a
test that documents the current behaviour.

| Bug | Effect | Test |
|---|---|---|
| `read stru` sets the displacement tensors of the first 5 × (atoms per cell) atoms only (`prep_anis`); later atoms keep tensor 1 | a saved supercell with anisotropic ADPs reads back with wrong tensors: a Fourier map changes by up to 3 %; isotropic crystals only report it wrong (`adps()`), their calculations use B | `test_save_anisotropic_supercell` |
| `save` writes one atom type per distinct tensor, `read stru` assigns atoms by name, B and occupancy | symmetry-rotated tensors are lost (236 of 384 atoms, anatase 4×4×2) | `test_save_and_read_keep_rotated_tensors` |
| `export cif` with `spcgr:original` lists the whole unit cell with the space group's operations | readers (DISCUS, gemmi) generate 80 atoms for anatase instead of 12 | `test_export_cif_original_space_group` |
| `export shelx` writes the whole unit cell but keeps LATT/SYMM and occupancy 0.5 | re-imported CuO has 32 atoms instead of 8 | `test_shelx_export_import_round_trip` |
| `.cell` files in the old format read fractions such as 1/3 as 1. | silently wrong positions | `test_cell_legacy_format_misreads_fractions` |

`dc.make_crystal` works around the first one. The first two together mean:
until fixed, keep anisotropic supercells as Python arrays (or as the
`.cell` file they came from), not as saved `.stru` files.

`make_crystal` takes **every atom**: the space group is recorded (Fourier
symmetry, export), not applied. That is right for anything built or changed
in Python, but it leaves out the crystallographer's usual description.

## 1. An asymmetric unit and a space group from Python

### What it is

A function that takes what a CIF or a `.cell` file holds, from Python:

```python
cr = dc.crystal_from_sites(
    lattice=(3.785, 3.785, 9.514, 90, 90, 90),
    space_group="I41/amd", origin_choice=1,          # or setting="abc", "cab" ...
    sites=[  # name, x, y, z, and optionally B or U, occupancy, label
        dc.Site("TI4+", 0, 0, 0,      u=(0.004, 0.004, 0.006, 0, 0, 0)),
        dc.Site("O2-",  0, 0, 0.2081, u=(0.003, 0.008, 0.005, 0, 0, 0)),
    ],
    supercell=(10, 10, 4),
)
```

DISCUS applies the space group (equivalent positions, rotated tensors), takes
out duplicates on special positions and builds the supercell, exactly as
`read cell` does for a file. Names in the method: `crystal_from_sites` (or
`dc.Structure(...).build(supercell)` if a separate asymmetric-unit object
turns out useful, see below).

### What it gives that `make_crystal` does not

* **The crystallographic description itself:** sites with their Wyckoff
  multiplicity, the space group applied by DISCUS's own routines (the ones
  the golden tests check against gemmi), and the symmetry operation each
  atom came from (`cr_iscat(2,i)`), which DISCUS uses for the anisotropic
  ADPs and for site-resolved analyses.
* **Symmetry-correct tensors:** a tensor given once per site is rotated for
  each equivalent atom (with `make_crystal`, Python would have to do it).
* **CIF without a file:** gemmi reads the CIF (any dialect gemmi handles),
  Python passes the asymmetric unit to DISCUS: no temporary file, full
  control over occupancies, labels and settings, and a clear error if a
  setting is not one DISCUS knows.
* **Sites as DISCUS sees them:** the per-site atom types that `.cell` files
  give (see the `types` label of `make_crystal`), and `as_*` (DISCUS's
  asymmetric-unit list), which some menus use.
* **A starting point for disorder:** site occupancies below 1 and mixed
  sites, kept per site, then decorated or relaxed by DISCUS (MMC, RMC).

### How DISCUS does it, and what the C function must repeat

`read cell` (`do_readcell` → `readcell` in `discus/prog/structur.f90`):

1. Header (`stru_readheader`): name, space group (`spcgr_no`: number,
   crystal system, origin choice, setting `cr_set`), lattice, atom types,
   B, occupancy, anisotropic tensors per type.
2. `setup_lattice` and `get_symmetry_matrices` (the operations of the space
   group in the chosen setting).
3. For each site: `firstcell` (into the first unit cell), then `symmetry`
   (`spcgr_apply.f90`): the equivalent positions, each with the index of its
   operation; a molecule in the asymmetric unit is expanded as a whole.
4. `test_identical`: atoms closer than a tolerance (special positions)
   are merged.
5. The supercell: the unit cell copied `cr_icc` times; `cr_ncatoms`,
   periodic boundaries (`chem_period`, `chem_quick`).
6. `prep_anis` (tensors per atom), `guess_atom_all`, `test_mole_gap`,
   stack reset.

A C function `dc_set_cell(lattice, spcgr, origin, setting, ntypes, names,
biso, occ, has_anis, anis, nsites, pos, type, icc)` calls these routines in
memory (no file), the same way `dc_set_crystal` follows `read stru`. The
difference: steps 3–5 are DISCUS's own symmetry code, not ours.

### Points to settle

* **Origin choice and settings** (e.g. `Fd-3m` origin 1 or 2, monoclinic
  `b`/`c` unique, rhombohedral `R` axes): DISCUS's `spcgr_no` decides from
  the symbol, the lattice and an optional origin parameter; the Python side
  must pass them explicitly and fail clearly otherwise.
* **Tolerance for special positions** (`test_identical`) and what to report
  when a site is duplicated (CIFs with rounded coordinates).
* **Partial occupancy:** DISCUS keeps one occupancy per atom type; mixed
  sites (two elements on one site) are two types on the same position.
* **Molecules in the asymmetric unit:** `read cell` expands them; the first
  version takes atoms only, molecules follow.
* **An asymmetric-unit object** (`dc.Structure`: lattice, space group,
  sites, with `.build(supercell)`, `.to_cif()`, `from_cif()` through gemmi):
  useful as the editable description users keep; decide after the first
  version.
* **spglib, the other direction:** finding the space group and asymmetric
  unit of a P1 crystal (`Crystal.find_symmetry()`), e.g. after ASE or a
  relaxation. Separate, small, optional.

### Tests

* The nine golden `.cell` structures: `crystal_from_sites` with the sites of
  the file gives identical arrays (positions, types, tensors, symmetry
  operation per atom) and identical powder, Fourier and PDF results.
* The same sites expanded by gemmi (already used in `test_wp3_structures`):
  the same set of positions for every space group in the test list.
* Anisotropic tensors on special positions: DISCUS's Wyckoff check
  (`test_symm`) accepts allowed tensors and refuses others with a clear
  error.
* CIF → gemmi → `crystal_from_sites` against `read_structure("x.cif")`.

**Effort:** about 2–3 days for the C function, the Python API and the tests;
the asymmetric-unit object and spglib about 1–2 days more.

## 2. The option E pilot, with benchmarks against the Fortran

Option E (`WP2_Binding_Decision.md`) is to move DISCUS's physics out of
Fortran gradually, behind the C interface. Before deciding anything larger,
the pilot measures what it costs and what it gives, on one kernel.

### Questions it answers

1. Can a rebuild give **the same results** as the Fortran (to rounding),
   checked by the golden tests?
2. How **fast** is it, on one core and on all cores, against the Fortran's
   OpenMP code, for small to very large crystals?
3. How much **effort** was it (hours, lines, build complexity), and how
   readable and maintainable is the result?
4. What does a **GPU** offer for this kind of kernel (later step)?

### The kernel: the Debye powder pattern

`powder_debye_hist_cart` (`discus/prog/powder_pdf_hist.f90`), which
`dc.powder(method="debye")` uses:

1. For every pair of atoms (O(N²), OpenMP over the first atom), the
   distance goes into a histogram per pair of atom types (bin width
   `pow_del_hist`).
2. The histograms are transformed: I(Q) is the sum over type pairs of
   f_i(Q) f_j(Q) times the displacement factors times Σ_r n_ij(r) sin(Qr)/(Qr),
   plus the self term. The Fortran takes sin from a lookup table
   (`sinetab`), not from the exact function.

Why this one: compact (a few hundred lines), compute-bound, already
parallel in the Fortran, covered by golden tests (`02.powder`, the
`dc.powder` tests), central for nanoparticles, and a classic GPU win. The
PDF (`dc.pdf`, the same histograms) comes almost for free once the
histogram is fast; it is the optional second kernel.

### What is built

| Version | How | Built / run |
|---|---|---|
| **Fortran** (reference) | the existing code through `dc.powder` | `libdiffuse_core`, OpenMP |
| **C++** | a C function `dc_debye_cxx(...)` with the same arrays in and out, OpenMP, in the same library | CMake, no new dependency |
| **NumPy** | vectorised over blocks of pairs; the readable reference | pure Python |
| **Numba** | the same loop as the Fortran, `@njit(parallel=True)`, `prange` | `pip install numba` |
| *later:* **GPU** | Numba CUDA or CuPy (same histogram), or JAX | needs a GPU machine |

All versions take the **same input**: the crystal's arrays from
`Crystal.atoms()` / `adps()` (now available), the same Q grid and the same
form factors (taken from DISCUS once, so the comparison is about the
kernel, not about tables). Double precision everywhere.

### Correctness first

* The histogram, bin by bin, against the Fortran's: the integer counts
  must be equal (same bins, same distances in double precision).
* I(Q) against the Fortran: equal to the accuracy of the Fortran's sine
  table (measured first, by running the Fortran transform with the exact
  sine once); each rebuild offers both, the table (for the comparison) and
  the exact sine. The rebuilds agree with each other to about 1e-10.
* The golden tests (`02.powder`, `test_diffuse_api`) run on each version.
* Edge cases: one atom type, many types, vacancies (type 0), partial
  occupancy, a single atom, atoms closer than one bin.

### The benchmarks

A script, `tools/benchmark_pilot.py`, kept in the repository:

* **Crystal sizes:** N = 10², 10³, 10⁴, 10⁵, 10⁶ atoms (built with
  `make_crystal`: fcc Ni and anatase supercells, and a random alloy), up to
  what finishes in reasonable time for O(N²) (10⁶ atoms: 5·10¹¹ pairs, run
  only for the fast versions and noted otherwise).
* **Grids:** Q points 1 000 and 10 000; histogram bins 0.01 and 0.001 Å.
* **Threads:** 1, and all cores (`OMP_NUM_THREADS`, `NUMBA_NUM_THREADS`
  set to the same value); a scaling curve 1, 2, 4, … cores on the largest
  size.
* **Repeats:** each case 5 times after one warm-up; report the median and
  the spread (quartiles). Numba's compile time is measured and reported
  separately, not hidden in the first run.
* **What is timed:** the kernel alone (histogram + transform) and the whole
  call (`dc.powder` including DISCUS's set-up), so set-up costs stay
  visible.
* **Memory:** peak resident memory per case, measured in a separate
  process (`resource.getrusage(RUSAGE_CHILDREN).ru_maxrss`).
* **Recorded with every result:** CPU model, cores, memory, OS, compilers
  and their flags (gfortran, g++), Python, NumPy, Numba versions, the
  DiffuseCode commit. Results go to a JSON file; the script also writes a
  Markdown table and plots (time against N on log axes, speed-up against
  the Fortran, scaling with cores).
* **Where:** correctness on every CI run (small sizes only); the
  benchmarks on the owner's machine(s) and, if available, a cluster node;
  later the same script on a GPU machine. CI runners are shared and noisy,
  so CI does not benchmark.

### Fairness rules

* The same algorithm (pair histogram, then transform) in every version; an
  algorithmically better version (e.g. a cut-off with periodic images, or
  an FFT-based method) is reported separately, not mixed into the
  comparison.
* The same compiler optimisation level for Fortran and C++ (`-O3`; with and
  without `-march=native` reported).
* No I/O inside the timing; warm caches; the machine otherwise idle.

### Effort record

For each version: hours spent, lines of code, new build or install
dependencies, and a short note on what was hard. This is what decides
option E as much as speed does.

### Deciding

The pilot ends with a short report (`migration/OptionE_Pilot_Report.md`)
and a recommendation. A rebuild of more kernels is worth starting if a
version

* matches the Fortran to rounding on all tests,
* is no more than about 1.5× slower than the Fortran on one core, and
  scales as well with cores (or is faster), and
* took effort in line with the estimate (about a week per kernel of this
  size), with code the owner finds clearer to maintain.

A GPU result changes the picture for large crystals and is reported on its
own.

### Plan and effort

| Days | Work |
|---|---|
| 1 | inputs (form factors, Q grid, histogram parameters) from DISCUS through the C interface; NumPy reference; correctness harness |
| 2–3 | C++ version behind a C function, OpenMP; tests |
| 4 | Numba version; tests |
| 5–6 | benchmark script, runs on the owner's machine, plots |
| 7 | report and recommendation |
| later | GPU version on a GPU machine; PDF kernel |

About 1.5 weeks in all, in line with the estimate in the binding decision.

## 3. Crystal methods from the Fortran

DISCUS does most things to "the crystal" through commands and menus; in
Python they belong on `Crystal` (or next to it). A scan of the DISCUS
command table (`discus/prog/discus_kdo.f90`), the command summary in
`discus/prog/appl_dis.hlp` and the sub-menus (save, export, import, chem,
plot) gives the list below. P1: most useful next; P3: later or niche.

### Writing the crystal

| DISCUS | Method | What it gives | Priority |
|---|---|---|---|
| `save` (`.stru`, `internal`, NeXus/HDF5) | `crystal.save(path)` | DISCUS's own format: everything (molecules, properties, ADPs, cells) | **done** |
| `export cif / shelx / rmcprofile / vasp / scatty` | `crystal.export(path, format=...)` | exchange formats; RMCProfile and Scatty are hard to get elsewhere | **done** |
| `plot` (xyz, cif, diamond, jmol, xbs, frames, kupl; selection, bonds, polyhedra) | `crystal.export_view(...)` | files for structure viewers; `view()` covers the interactive case | P3 |
| `storage` (internal, in memory) | `crystal.snapshot()` or internal | fast copy/restore without a file | P2 |

`to_ase()` already writes dozens of formats; the DISCUS writers matter for
`.stru`, RMCProfile and Scatty, and for keeping DISCUS's ADPs and molecules.

### Questions about the crystal (methods returning data)

| DISCUS | Method | What it gives | Priority |
|---|---|---|---|
| `asym` | `crystal.asymmetric_unit()` | sites: name, position, B, occupancy | **P1** (with `crystal_from_sites`) |
| `wyckoff` | `crystal.site_symmetry(position)` | symmetry operations at a position, multiplicity | P2 |
| `show` (cell, metric ...) | `metric_tensor`, `reciprocal_lattice` | (`volume`, `lattice_vectors()` exist) | P2 |
| chem `elem` | `crystal.composition()` | counts and fractions per atom type | **P1** |
| `find` | `crystal.find_atoms(center, rmin, rmax, types=)` | atoms around a point | P2 |
| chem `neig`, `env`, `coordination` | `crystal.neighbours(...)`, `crystal.coordination(...)` | neighbour lists and coordination numbers, DISCUS's definitions | **P1** |
| chem `blen`, `bang` | `crystal.bond_lengths(...)`, `crystal.bond_angles(...)` | bond-length and bond-angle distributions | **P1** |
| chem `bval` | `crystal.bond_valence(...)` | bond-valence sums | P2 |
| chem `corr`, `homo`, `disp` | `crystal.correlations(...)`, `crystal.displacement_correlations(...)` | occupancy and displacement correlations: the core of short-range-order analysis | **P1** |
| chem `aver` | `crystal.average_structure()` | the average unit cell of a disordered crystal | P2 |
| `connectivity` | `crystal.connectivity(...)` | DISCUS connectivity definitions (used by MMC and chem) | P2 |
| `property` | `crystal.properties()`, selection by property | the per-atom flags (normal, molecule, surface, inside, outside ...) | P2 |
| `d2r`, `r2d`, `vprod` | `crystal.to_reciprocal(v)`, `to_direct(v)`, `cross(u, v)` | vector maths in the crystal metric (could be pure Python) | P3 |
| `project` | `crystal.project(...)` | projections in real and reciprocal space | P3 |

### Changing the crystal

| DISCUS | Method | What it gives | Priority |
|---|---|---|---|
| `boundary` (hkl planes, sphere, ellipsoid, cylinder) | `crystal.cut(shape=...)` | nanoparticles of any shape from a supercell | **P1** |
| `replace` | `crystal.replace(old, new, fraction=)` | random substitution: alloys, doping | **P1** |
| `remove` | `crystal.remove(selection)` | vacancies; removing atoms or molecules | **P1** |
| `therm` | `crystal.thermal_displacements(seed=)` | displace atoms according to B / U: a thermal snapshot for diffuse scattering | **P1** |
| `transform` | `crystal.transform(matrix, origin)` | a new unit cell (supercells, other settings) | **P1** |
| `symm` | `crystal.apply_symmetry(...)` | rotate / mirror / translate the crystal or parts of it | P2 |
| `insert`, `append`, `copy`, `switch`, `kick` | `crystal.insert(...)` ... | add, copy, swap atoms and molecules | P2 |
| `perioditize` | `crystal.perioditize()` | make the crystal periodic | P2 |
| `spacegroup` | `crystal.set_space_group(...)` | change the recorded space group | P2 |
| `molecularize`, `demolecularize`, `change` | `crystal.group_molecules(...)` ... | build or dissolve molecules; molecule character | P2 |
| `domain` | `crystal.insert_domains(...)` | domains of another structure | P3 |
| `surface` | `crystal.surface(...)` | find and mark surface atoms | P3 |
| `decorate` | `crystal.decorate(...)` | decorate with molecules or ligands | P3 |
| `waves`, `shear` | `crystal.modulate(...)`, `crystal.shear(...)` | modulation waves, shear | P3 |
| `super` | `crystal.superspace(...)` | superspace description | P3 |
| `purge` | – | DISCUS says "DO NOT USE" | skip |

### Better outside `Crystal`

* **Making a crystal** (functions, next to `read_structure`, `make_crystal`,
  `from_ase`): `read` (cell, stru, lcell, free), `import` (cif, shelx,
  rmcprofile, cmaker, lammps, nexus), `free` → `dc.empty_crystal(lattice,
  space_group)`.
* **Calculations** (functions taking a crystal): `fourier`, `powder`, `pdf`
  (exist), `patterson`, `inverse`, `diff-four`, `exp2pdf`; `output` writes
  results, not the crystal. Methods such as `crystal.powder(...)` could be
  added as shortcuts.
* **Simulations** (classes returning a new `Crystal`): `mmc`, `rmc`,
  `stack` (stacking faults), `place_molecule`.
* **Session settings**: `reset`, `seed`, `define`, `echo`.

### Design questions to settle first

1. **New crystal or change in place?** DISCUS changes its one crystal; in
   Python, `cut()`, `replace()` and the like should return a new `Crystal`
   and leave the original as it was. Possible now: after the DISCUS
   command the result comes back through the arrays (`make_crystal` in
   reverse).
2. **Selections:** DISCUS's `select` / `deselect` state per menu becomes
   arguments: atom names, a boolean mask, or properties.
3. **Analysis results** as objects with arrays and `plot()`, like
   `PowderPattern` and `FourierMap`.
4. **Python or DISCUS?** Some are easy in NumPy now (remove or replace by
   mask, composition, vector maths). Wrap DISCUS where it keeps molecules,
   properties and ADPs consistent, or where the physics lives
   (correlations, `therm`, `boundary`, MMC).

**First batch (proposal):** `save()` and `export()` (done, Nov 2026);
`composition()` and `asymmetric_unit()`; `cut()`, `replace()`, `remove()`,
`thermal_displacements()`; then the chem analysis (neighbours, bond
lengths, correlations) as one block.

## Order of the next steps (proposal)

1. `crystal_from_sites` (asymmetric unit + space group), with the gemmi
   CIF path; `asymmetric_unit()` and `composition()` with it.
2. The rest of the first batch of crystal methods (section 3): `cut()`,
   `replace()`, `remove()`, `thermal_displacements()`, then the chem
   analysis.
3. The option E pilot as above (it needs `make_crystal`, `atoms()` and
   `adps()`, which are in place).
4. Then, depending on the pilot, WP4 (single-crystal Fourier as the first
   vertical slice) on the Fortran or on a rebuilt kernel.

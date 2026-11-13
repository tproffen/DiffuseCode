<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DiffuseCode Python Migration</title>
<style>
:root{
  --bg:#fbfaf7; --surface:#ffffff; --text:#1d1f23; --muted:#5d636e; --border:#e2e0da;
  --accent:#2f5d8a; --accent-soft:#e8eff6; --good:#2e7d4f; --bad:#a6432f; --code:#f3f1ec;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#16181c; --surface:#1e2126; --text:#e6e6e3; --muted:#9aa1ab; --border:#33373e;
    --accent:#7fb0de; --accent-soft:#22303e; --good:#6cc191; --bad:#e08a75; --code:#25292f;
  }
}
:root[data-theme="dark"]{
  --bg:#16181c; --surface:#1e2126; --text:#e6e6e3; --muted:#9aa1ab; --border:#33373e;
  --accent:#7fb0de; --accent-soft:#22303e; --good:#6cc191; --bad:#e08a75; --code:#25292f;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);
  font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:920px;margin:0 auto;padding:40px 16px 80px}
h4{font-size:.98rem;margin:18px 0 6px}
.wp{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:14px 18px;margin:14px 0}
.wp h3{margin:0 0 6px}
.wp .meta{color:var(--muted);font-size:.88rem;margin:0 0 8px}
.wp dl{margin:6px 0 0;display:grid;grid-template-columns:max-content 1fr;gap:4px 14px}
.wp dt{font-weight:600;color:var(--muted);font-size:.9rem}
.wp dd{margin:0}
@media (max-width:560px){.wp dl{grid-template-columns:1fr}.wp dt{margin-top:6px}}
.bar{display:inline-block;height:10px;border-radius:3px;background:var(--accent);vertical-align:middle}
.checklist{padding-left:0}
.checklist li{list-style:none}
.checklist li::before{content:"☐ ";color:var(--accent)}
header{border-bottom:1px solid var(--border);padding-bottom:20px;margin-bottom:28px}
h1{font-size:2rem;line-height:1.2;margin:0 0 8px}
.sub{color:var(--muted);margin:0}
h2{font-size:1.35rem;margin:44px 0 12px;padding-top:8px;border-top:1px solid var(--border)}
h3{font-size:1.08rem;margin:24px 0 8px}
p,li{max-width:75ch}
ul{padding-left:1.3em}
li{margin:4px 0}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.88em;
  background:var(--code);padding:1px 5px;border-radius:4px}
pre{background:var(--code);border:1px solid var(--border);border-radius:8px;padding:16px;
  overflow-x:auto;font-size:.85rem;line-height:1.45}
pre code{background:none;padding:0}
.table-wrap{overflow-x:auto;border:1px solid var(--border);border-radius:8px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:.93rem}
th,td{text-align:left;padding:9px 12px;border-bottom:1px solid var(--border);vertical-align:top}
th{background:var(--accent-soft);font-weight:600}
tr:last-child td{border-bottom:none}
td.num{text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}
nav.toc{background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:14px 18px}
nav.toc ol{margin:0;padding-left:1.3em}
nav.toc a{color:var(--accent);text-decoration:none}
nav.toc a:hover{text-decoration:underline}
.options{display:grid;gap:14px;grid-template-columns:repeat(auto-fit,minmax(260px,1fr))}
.card{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:16px 18px}
.card.rec{border:2px solid var(--accent)}
.card h3{margin-top:0}
.tag{display:inline-block;font-size:.75rem;font-weight:600;padding:2px 8px;border-radius:999px;
  background:var(--accent-soft);color:var(--accent);margin-left:6px;vertical-align:middle}
.pro{color:var(--good);font-weight:600}
.con{color:var(--bad);font-weight:600}
.callout{background:var(--accent-soft);border-left:4px solid var(--accent);border-radius:6px;
  padding:14px 18px;margin:18px 0}
footer{margin-top:56px;color:var(--muted);font-size:.85rem;border-top:1px solid var(--border);padding-top:14px}
@media print{body{background:#fff;color:#000} .card,.table-wrap{break-inside:avoid}}
</style>
</head>
<body>
<main>
<header>
  <h1>DiffuseCode: moving to a Python framework</h1>
  <p class="sub">A high-level look at the current Fortran code and options for separating compute from the command language, graphics and I/O, with a web front end as a later goal.</p>
</header>

<nav class="toc">
<ol>
  <li><a href="#today">What the code is today</a></li>
  <li><a href="#target">Target architecture</a></li>
  <li><a href="#options">High-level options</a></li>
  <li><a href="#route">Suggested route</a></li>
  <li><a href="#risks">Risks and early decisions</a></li>
  <li><a href="#bottom">Bottom line</a></li>
  <li><a href="#planb">Detailed plan for option B</a>
    <ol>
      <li><a href="#b-goals">Goals and non-goals</a></li>
      <li><a href="#b-principles">Guiding principles</a></li>
      <li><a href="#b-layout">Package layout</a></li>
      <li><a href="#b-interface">Interface design</a></li>
      <li><a href="#b-api">Target Python API (sketch)</a></li>
      <li><a href="#b-wp">Work packages</a></li>
      <li><a href="#b-timeline">Timeline</a></li>
      <li><a href="#b-testing">Testing strategy</a></li>
      <li><a href="#b-parallel">Parallelism and performance</a></li>
      <li><a href="#b-risks">Risks and mitigations</a></li>
      <li><a href="#b-90">First 90 days</a></li>
      <li><a href="#b-decisions">Decisions needed</a></li>
    </ol>
  </li>
</ol>
</nav>

<h2 id="today">1. What the code is today</h2>
<p><strong>Size and makeup</strong> (Fortran line counts, not counting the bundled <code>finufft</code> and <code>h5fortran</code>)</p>
<div class="table-wrap">
<table>
<thead><tr><th>Component</th><th>Files</th><th>Lines</th><th>What it does</th></tr></thead>
<tbody>
<tr><td><code>discus</code></td><td class="num">~200</td><td class="num"><strong>~156k</strong></td><td>Builds structures, simulates defects (MC/RMC/MMC, domains, stacking faults, molecules, surfaces, waves), and computes Fourier, powder, PDF, 3D-PDF and Patterson</td></tr>
<tr><td><code>lib_f90</code></td><td class="num">~150</td><td class="num">~70k</td><td>Shared base: the command interpreter, variables, macros, loops, expression evaluator, error lists, HDF5/NeXus I/O, FFT, random numbers, MPI/OpenMP glue, sockets, readline, and a Fortran↔Python bridge (<code>forpy</code>)</td></tr>
<tr><td><code>kuplot</code></td><td class="num">66</td><td class="num">~38k</td><td>Plotting (PGPLOT on X11, Windows or PNG) plus data loading, maths and least-squares fitting</td></tr>
<tr><td><code>diffev</code></td><td class="num">39</td><td class="num">~10k</td><td>Differential-evolution refinement, with MPI for parallel cost functions</td></tr>
<tr><td><code>refine</code></td><td class="num">31</td><td class="num">~8k</td><td>Least-squares refinement that drives <code>discus</code> macros</td></tr>
<tr><td><code>suite</code></td><td class="num">16</td><td class="num">~3k</td><td>One umbrella program that joins all of the above and switches between their menus</td></tr>
<tr><td><code>mixscat</code>, <code>experi</code>, <code>nexus</code></td><td class="num">small</td><td class="num">~7k</td><td>Extra tools</td></tr>
</tbody>
</table>
</div>
<p>There is also C/C++: <code>finufft</code>, plus readline, sockets and glob helpers.</p>

<h3>Architecture</h3>
<ul>
  <li><strong>Everything runs through the command language.</strong> Each program and each sub-menu has its own <code>*_loop</code> → <code>get_cmd</code> → <code>*_kdo</code>. In <code>discus</code>, 30 files contain such a menu loop, and the top-level <code>discus_kdo</code> dispatches about 70 commands. In <code>discus/prog</code> alone the code makes about 1,700 <code>str_comp</code> keyword comparisons, about 1,100 <code>get_params</code>/<code>ber_params</code> calls and about 3,400 <code>ier_num</code> error checks.</li>
  <li><strong>The command language is a whole small language.</strong> It has variables (<code>i[]</code>, <code>r[]</code>, <code>res[]</code>), an expression evaluator, <code>do</code>/<code>if</code> blocks, macros with parameters, <code>branch</code> between programs, <code>system</code>, <code>socket</code> and <code>learn</code>.</li>
  <li><strong>All state is global.</strong> Crystal, Fourier, powder, MC and similar data live in roughly 75 <code>*_mod</code> modules (for example <code>crystal_mod</code> holds <code>cr_pos</code>, <code>cr_iscat</code> and so on). Compute routines read their inputs from these modules and write results back into them, not through arguments.</li>
  <li><strong>Results come back as text and through <code>res_para</code>.</strong> There are about 1,250 <code>WRITE(output_io…)</code> statements in <code>discus</code>, and numeric results go into <code>res[]</code> for macros to use.</li>
  <li><strong>The programs talk to each other through macros.</strong> <code>refine</code> and <code>diffev</code> work by running user macros that call <code>discus</code> and <code>kuplot</code>, then reading the results back from variables or files. MPI spreads those macro runs across processes.</li>
  <li><strong>Graphics live in <code>kuplot</code>, which calls PGPLOT directly.</strong> A no-graphics build already exists (<code>*_nographic.f90</code> and the <code>DIFFUSE_GRAPHICS=OFF</code> option).</li>
  <li><strong>Some useful pieces already exist:</strong>
    <ul>
      <li><code>python/suite_python.f90</code> (built with f2py) offers <code>initialize_suite</code>, <code>execute_macro</code>, <code>get_data</code> and <code>set_value</code>. Python can call it, but only by passing command strings.</li>
      <li><code>discus_py.f90</code>, <code>kuplot_py.f90</code> and <code>diffev_py.f90</code> expose <code>command()</code>, <code>get</code> and <code>send</code>.</li>
      <li><code>forpy</code> lets the Fortran side call Python.</li>
      <li>There is a CUDA kernel for the structure factor, OpenMP in the Fourier, PDF and MMC code, and FINUFFT for non-uniform FFTs.</li>
    </ul>
  </li>
  <li><strong>There are almost no tests.</strong> The repo has only <code>sysmac/systest.mac</code> files, so there is no regression suite.</li>
</ul>

<div class="callout">
<strong>What this means for a migration.</strong> The compute kernels themselves are mostly decoupled. For example, <code>four_strucf</code> only uses <code>discus_config_mod</code>, <code>diffuse_mod</code> and <code>precision_mod</code>. The coupling sits in the menu layer around them: it parses commands, sets global module variables, calls the kernel, then prints the results or stores them in <code>res[]</code>. So the usual entry point is not "the algorithm" but "the global state that menu commands have built up". That global state is the main obstacle, more than the Fortran language itself.
</div>

<h2 id="target">2. Target architecture (shared by all options)</h2>
<pre><code> Web UI / Jupyter / scripts / CLI
            │
   Python API layer   (Crystal, Structure, Fourier, Powder, PDF, MC, Fit objects)
            │                     │                    │
   I/O (CIF, HDF5/NeXus,     Plotting (matplotlib,   Refinement (scipy.optimize,
   ASE/pymatgen, xarray)     plotly, a web viewer)   lmfit, DE, dask/MPI)
            │
   Compute core: no global I/O and no parsing; arrays in, arrays out
   (the existing Fortran kernels behind a thin interface, or Python versions)</code></pre>
<p>The rules this architecture follows:</p>
<ul>
  <li>The core never prints, never parses text and never draws.</li>
  <li>The core reports errors through return codes or exceptions, not <code>ier_num</code> and a global error state.</li>
  <li>Plotting, file formats and refinement become ordinary Python libraries.</li>
  <li>The command language becomes optional: a compatibility layer at most, not the main interface.</li>
</ul>

<h2 id="options">3. High-level options</h2>
<div class="options">
  <div class="card">
    <h3>A. Wrap the current interpreter <span class="tag">weeks</span></h3>
    <p>Extend the existing <code>suite_python</code> and f2py bridge. Python sends command strings, and data comes back as NumPy arrays through <code>get_data</code>, <code>get_r</code> and <code>set_value</code>.</p>
    <p><span class="pro">Good:</span> almost no Fortran changes, works now, and is useful for Jupyter.</p>
    <p><span class="con">Bad:</span> the command language remains the real API, global state stays, and you can't run two structures at once. Errors and output are still text. A web UI on top of this would be fragile.</p>
    <p><strong>Best use:</strong> a short-term stopgap, and a way to build a regression-test corpus.</p>
  </div>
  <div class="card rec">
    <h3>B. Fortran compute behind a clean Python API <span class="tag">recommended</span></h3>
    <ol>
      <li>Find the computational entry points: structure factor and Fourier, powder, PDF, 3D-PDF, Patterson, MC/RMC/MMC moves and energies, symmetry and space-group expansion, stacking faults, domains, and so on.</li>
      <li>For each one, write a <code>bind(C)</code> wrapper that takes explicit arguments (positions, types, cell, q-grid, form factors) and returns arrays. At first, the wrapper can copy its arguments into the existing global modules and call the current routine. Later, refactor each routine so it takes arguments directly.</li>
      <li>Build it as a Python extension, either with f2py or with <code>bind(C)</code> + nanobind/pybind11/ctypes, packaged by scikit-build-core or meson-python. The second route handles derived types and allocatable arrays more cleanly and gives proper wheels.</li>
      <li>Write a Pythonic object model on top: <code>Crystal</code>, <code>Fourier(…)</code>, <code>Powder(…)</code> and similar, using NumPy and xarray.</li>
      <li>Replace <code>kuplot</code> plotting with matplotlib/plotly, <code>refine</code> with scipy/lmfit, and <code>diffev</code> with scipy's DE or a small Python DE that uses dask/mpi4py/Ray for parallel runs. Move file I/O to h5py, gemmi/ASE/pymatgen (CIF) and nexusformat.</li>
      <li>Optional: a thin "macro compatibility" translator for existing user scripts, or keep the old <code>discus_suite</code> binary as a legacy build alongside.</li>
    </ol>
    <p><span class="pro">Good:</span> keeps decades of validated physics and the existing speed (OpenMP, CUDA, FINUFFT). Each step can be tested against the old binary. It also leads naturally to a web or service backend.</p>
    <p><span class="con">Bad:</span> removing the global state needs care. Several modules, especially <code>mmc</code> (~9k lines), <code>structur</code>, <code>chem</code> and <code>place_molecule</code>, are tightly tied to their menus. Building and packaging Fortran plus HDF5 and FFTW for wheels takes real work, though cibuildwheel can handle it.</p>
  </div>
  <div class="card">
    <h3>C. Full rewrite in Python <span class="tag">years</span></h3>
    <p>Reimplement the physics in pure Python, with Numba, JAX or CuPy for the hot loops: structure factors, PDF histograms, MC.</p>
    <p><span class="pro">Good:</span> one language, simplest packaging, the easiest code for others to contribute to, and GPU support and automatic differentiation (JAX) for refinement.</p>
    <p><span class="con">Bad:</span> about 160k lines of domain logic to re-derive and re-validate, with no test suite as a safety net. It would take years of effort and carries a high risk of subtle physics regressions.</p>
    <p><strong>Best use:</strong> only for selected parts, for example new features, or kernels that are short and easy to check, such as the Debye sum, PDF histogramming and the structure factor on a grid.</p>
  </div>
  <div class="card">
    <h3>D. Hybrid that grows into a rewrite <span class="tag">long term</span></h3>
    <p>Start with B. Port each module to Python when it is cheaper to rewrite than to untangle: I/O, plotting, refinement, set-up and geometry tools, symmetry (using spglib/gemmi). Keep Fortran, C++ or CUDA for the heavy kernels (Fourier, MC energy loops, powder/PDF) for as long as they are faster.</p>
    <p>This is option B followed by C, done as a long-term plan.</p>
  </div>
</div>

<h2 id="route">4. Suggested route (B, then D)</h2>
<div class="table-wrap">
<table>
<thead><tr><th>Phase</th><th>Work</th><th>Outcome</th></tr></thead>
<tbody>
<tr><td><strong>0. Safety net</strong></td><td>Use the existing <code>suite_python</code> (option A) to run a set of representative macros (Manual examples, user workflows) and save the outputs as golden HDF5/NumPy files. Add CI with pytest.</td><td>Regression tests: the precondition for everything else</td></tr>
<tr><td><strong>1. Packaging</strong></td><td>Add a <code>pyproject.toml</code> with scikit-build-core around the existing CMake. Build one importable <code>diffuse._core</code> with <code>DIFFUSE_GRAPHICS=OFF</code>, no X11 and no readline.</td><td><code>pip install</code> works</td></tr>
<tr><td><strong>2. Data model</strong></td><td>Python <code>Crystal</code>/<code>Structure</code> classes holding NumPy arrays, plus Fortran getters and setters for <code>crystal_mod</code> arrays (positions, types, cell, properties, molecules).</td><td>Structures move between Python and Fortran without files</td></tr>
<tr><td><strong>3. Kernel APIs</strong></td><td>Fourier/structure factor, then powder (Debye and FFT), then PDF, then 3D-PDF/Patterson. Each gets a pure-argument entry point and is checked against the golden data.</td><td>The main calculations can be called from Python</td></tr>
<tr><td><strong>4. Remove pieces from Fortran</strong></td><td>Replace plotting with matplotlib/plotly, refinement with lmfit/scipy plus Python DE, and file I/O with h5py/gemmi/ASE. Stop maintaining <code>kuplot</code>, <code>refine</code> and <code>diffev</code> as Fortran menus.</td><td>Smaller Fortran code, Python-native workflows</td></tr>
<tr><td><strong>5. Simulation engines</strong></td><td>MC/MMC/RMC, domains, stacking, surfaces, molecules. These are the hardest because they hold the most state. Consider a callback or "move/energy plug-in" interface, or port parts to Numba.</td><td>Defect modelling from Python</td></tr>
<tr><td><strong>6. Front ends</strong></td><td>Jupyter widgets first. Then a service layer (FastAPI with a job queue such as Celery, RQ or Dask, results in HDF5/Zarr) and a web UI (React with plotly or three.js for 3D reciprocal space).</td><td>Web front end</td></tr>
<tr><td><strong>7. Legacy</strong></td><td>Freeze <code>discus_suite</code> as the classic interface, or offer a macro→Python translator for common commands.</td><td>Existing users keep working</td></tr>
</tbody>
</table>
</div>

<h2 id="risks">5. Main risks and decisions to make early</h2>
<ul>
  <li><strong>Global state vs. multiple instances.</strong> A web service needs several independent sessions. With global Fortran modules, each session needs its own process, which is acceptable at first and naturally handled by workers. True re-entrancy requires refactoring into derived types.</li>
  <li><strong>Error handling.</strong> <code>ier_num</code>/<code>ier_typ</code> and some <code>STOP</code> paths must become return codes turned into Python exceptions. A Fortran <code>STOP</code> would kill the Python process.</li>
  <li><strong>Output.</strong> The many <code>WRITE(output_io…)</code> calls need a quiet or logging mode.</li>
  <li><strong>Precision and compiler matrix.</strong> You currently support gfortran, ifort, PGI, OpenMP and optional CUDA. Wheels mean picking a supported subset.</li>
  <li><strong>User base.</strong> Macros are how existing users work, so decide early whether you need backward compatibility (a translator, or the old binary frozen) or will make a clean break.</li>
  <li><strong>Licensing and dependencies.</strong> PGPLOT disappears (a good thing). FFTW, HDF5 and FINUFFT stay as build dependencies.</li>
</ul>

<h2 id="bottom">6. Bottom line</h2>
<div class="callout">
<p style="margin:0">The physics code is valuable and fairly modular, but it is wrapped in a command interpreter and global-state layer that makes up a large share of the code. The recommendation is <strong>option B</strong>, a Python API over the existing Fortran kernels, built in phases. Start by building the regression-test corpus and Python packaging, replace plotting, I/O and refinement with standard Python libraries, and port kernels to Python only where it pays off. That gives a scriptable Python package in months rather than years, and a natural backend for a web front end later.</p>
</div>

<h2 id="planb">7. Detailed plan for option B</h2>
<p>Option B keeps the Fortran physics and replaces everything around it, one piece at a time. Each step must leave a working <code>discus_suite</code> and a passing test suite. The plan is written for a small team (about 1–2 developers plus a domain expert to review results). All durations are rough estimates.</p>

<h3 id="b-goals">7.1 Goals and non-goals</h3>
<div class="table-wrap">
<table>
<thead><tr><th>Goals</th><th>Non-goals (for option B)</th></tr></thead>
<tbody>
<tr><td>
<ul>
  <li><code>pip install diffuse</code> gives a Python package with no X11, PGPLOT or readline.</li>
  <li>Structures, Fourier, powder, PDF and defect simulations can be run from Python with NumPy arrays in and out. No command strings and no temporary files.</li>
  <li>Plotting, file I/O and refinement are handled by standard Python libraries.</li>
  <li>Errors become Python exceptions, and output goes to Python <code>logging</code>.</li>
  <li>The package works as the compute backend for a later web service.</li>
</ul>
</td><td>
<ul>
  <li>Rewriting the numerical kernels in Python (that is option C/D, case by case).</li>
  <li>Making the Fortran core fully re-entrant in the first pass. Isolation comes from processes at first.</li>
  <li>Reproducing every command-language feature (loops, <code>branch</code>, <code>socket</code>, <code>learn</code>). Python replaces them.</li>
  <li>Dropping <code>discus_suite</code>. It stays as a legacy build until users have moved.</li>
</ul>
</td></tr>
</tbody>
</table>
</div>

<h3 id="b-principles">7.2 Guiding principles</h3>
<ul>
  <li><strong>Strangler pattern.</strong> New Python entry points wrap existing routines. The command menus keep calling the same routines, so old and new paths share one implementation and both stay correct.</li>
  <li><strong>Golden data first.</strong> No refactor lands without a regression test that compares against output from the current release.</li>
  <li><strong>Thin, flat C boundary.</strong> The only thing that crosses from Fortran to Python is plain C types and contiguous arrays (<code>iso_c_binding</code>). No Fortran derived types, strings or allocatables cross it.</li>
  <li><strong>Python owns the data model.</strong> The truth lives in Python objects. Fortran modules are a cache that gets filled before a computation runs. Later refactors can remove that cache one module at a time.</li>
  <li><strong>Keep the menu layer, stop growing it.</strong> New features are added as Python API first, and as a menu command only if really needed.</li>
</ul>

<h3 id="b-layout">7.3 Package layout</h3>
<pre><code>DiffuseCode/
├── (existing Fortran trees: discus/ lib_f90/ kuplot/ diffev/ refine/ suite/ …)
├── capi/                        # NEW: Fortran bind(C) shims + C header
│   ├── dc_capi_core.f90         #   init, reset, errors, log redirection
│   ├── dc_capi_crystal.f90      #   get/set cell, atoms, scattering types, molecules
│   ├── dc_capi_fourier.f90      #   single-crystal Fourier / structure factor
│   ├── dc_capi_powder.f90
│   ├── dc_capi_pdf.f90
│   └── diffuse_capi.h
├── pyproject.toml               # NEW: scikit-build-core, drives the existing CMake
└── python/diffuse/              # NEW: the Python package
    ├── _core.cpp                #   nanobind module over diffuse_capi.h
    ├── __init__.py
    ├── errors.py                #   DiffuseError hierarchy mapped from ier_num/ier_typ
    ├── crystal.py               #   Crystal / Structure / Molecule classes
    ├── io/                      #   CIF (gemmi), DISCUS .stru, HDF5/NeXus (h5py), xyz, …
    ├── scattering/              #   fourier.py, powder.py, pdf.py, patterson.py
    ├── defects/                 #   replace, domains, stacking, surface, waves, …
    ├── mc/                      #   MMC / RMC drivers
    ├── refine/                  #   lmfit + differential-evolution wrappers
    ├── plot/                    #   matplotlib / plotly helpers (replaces kuplot)
    └── legacy.py                #   run_macro() via suite_python, for migration</code></pre>
<p>CMake already builds <code>lib_f90</code>, <code>discus_all</code> and the others as libraries, so the new extension can link them directly. One build step needs untangling: <code>discus</code> currently includes and links <code>kuplot</code>, so the dependency on <code>kuplot</code>'s plotting code has to be cut or replaced with stubs for the Python build. The existing <code>*_nographic.f90</code> files give a head start.</p>

<h3 id="b-interface">7.4 Interface design</h3>
<h4>The C-ABI shim</h4>
<p>Each shim function takes plain arguments, fills the existing global modules, calls the existing routine and copies results back out. Example for the structure:</p>
<pre><code>module dc_capi_crystal
  use iso_c_binding
  use crystal_mod
  use discus_allocate_appl_mod
  implicit none
contains
  integer(c_int) function dc_set_atoms(n, pos, iscat, prop, surf) &amp;
       bind(C, name="dc_set_atoms")
    integer(c_int), value         :: n
    real(c_double), intent(in)    :: pos(3, n)     ! fractional coordinates
    integer(c_int), intent(in)    :: iscat(n), prop(n), surf(n)
    call dc_clear_error()
    if (n &gt; NMAX) call alloc_crystal_nmax(n)       ! grow existing arrays
    cr_natoms        = n
    cr_pos(:, 1:n)   = pos
    cr_iscat(1:n,1)  = iscat
    cr_prop(1:n)     = prop
    ! … surface, molecule info, then update derived quantities
    dc_set_atoms = dc_status()                     ! 0 = ok, else ier_num
  end function
end module</code></pre>
<p>The same pattern applies to <code>dc_set_cell</code>, <code>dc_set_scat_types</code>, <code>dc_get_atoms</code> and others. The existing <code>discus_allocate_appl_mod</code> already contains the resizing routines (<code>alloc_crystal_nmax</code>, <code>alloc_crystal_scat</code>, <code>alloc_diffuse_four</code> and so on).</p>

<h4>Calls into compute routines</h4>
<ul>
  <li>Each <code>dc_*_run</code> function sets the module variables that the menu would normally set (for example the corners, increments and radiation in <code>diffuse_mod</code>), then calls <code>four_run</code>, <code>powder_run</code>, <code>pdf_run</code> and so on.</li>
  <li>Results are copied into caller-supplied NumPy buffers (<code>csf</code>, <code>dsi</code>, powder and PDF arrays), or exposed through a getter that reports their shape.</li>
  <li>The menu <code>set</code> commands contain validation logic. Where it matters, move that logic into small shared routines so both the menu and the shim call them, instead of copying it.</li>
</ul>

<h4>Errors</h4>
<ul>
  <li>Every shim function starts by clearing <code>ier_num</code>/<code>ier_typ</code> and returns them at the end. The C++ layer turns non-zero values into a Python exception and takes the text from the existing error lists (<code>errlist</code>).</li>
  <li>About 30 <code>STOP</code> statements exist, for example in <code>four_strucf_mod.f90</code> and the <code>assert*</code> helpers in <code>utilities.f90</code>. Each one must become an error return. Otherwise a single bad input ends the Python process.</li>
  <li>Any path that would drop into an interactive prompt on error (as the current <code>discus_py</code> <code>command()</code> does) must be disabled when called from the API.</li>
</ul>

<h4>Output and logging</h4>
<ul>
  <li>Add a global "API mode" flag in <code>prompt_mod</code>. When it is set, <code>output_io</code> goes to a buffer or <code>/dev/null</code>, and selected messages are passed to a C callback that forwards them to Python's <code>logging</code>.</li>
  <li>This needs no changes to the ~1,250 <code>WRITE</code> statements; only the unit they write to changes.</li>
</ul>

<h4>State and sessions</h4>
<ul>
  <li>One Python process holds one Fortran state. <code>dc_reset()</code> calls the existing reset routines (<code>discus_reset_all</code>, <code>fourier_reset</code>, <code>powder_reset</code>, <code>pdf_reset</code>) so that each calculation starts clean.</li>
  <li>A Python <code>Session</code> object can run the core in a worker process (multiprocessing or loky). That gives isolation for notebooks running several models and for the web service, without making the Fortran re-entrant.</li>
  <li>Longer term, the existing <code>cl_cryst</code> class (<code>class_crystal.f90</code>) is a starting point for passing a crystal as a derived type instead of through <code>crystal_mod</code>.</li>
</ul>

<h3 id="b-api">7.5 Target Python API (sketch)</h3>
<pre><code>import diffuse as dc

cr = dc.io.read_cif("PbTiO3.cif")                 # or dc.io.read_stru("old.stru")
sc = cr.expand(supercell=(20, 20, 20))            # space-group expansion + supercell

rng = dc.Random(seed=42)
sc = dc.defects.replace(sc, "Ti", "Zr", fraction=0.1, rng=rng)
sc = dc.mc.mmc(sc, targets=[dc.mc.Correlation("Ti", "Zr", neighbour=1, value=-0.3)],
               cycles=50, rng=rng)

f = dc.scattering.single_crystal(sc, radiation="xray", wavelength=0.71,
                                 plane=dc.Plane(corner=(-5,-5,0), right=(5,-5,0),
                                                top=(-5,5,0), npoints=(401,401)))
f.intensity.plot()                                # xarray.DataArray with hkl coordinates

p = dc.scattering.powder(sc, qmin=0.5, qmax=25, dq=0.01, method="debye")
g = dc.scattering.pdf(sc, rmax=30, dr=0.01)

fit = dc.refine.differential_evolution(model=my_model, data=obs,
                                       params={"frac": (0.0, 0.3)}, workers=8)</code></pre>
<p>The results are <code>xarray</code> objects with coordinates and metadata, and they can be saved straight to HDF5/NeXus.</p>

<h3 id="b-wp">7.6 Work packages</h3>

<div class="wp">
  <h3>WP0 · Regression corpus</h3>
  <p class="meta">Months 0–3 · prerequisite for everything else</p>
  <dl>
    <dt>Tasks</dt><dd>Collect 30–60 representative macros (Manual examples, tutorials, real user workflows). Run them with the current release through <code>suite_python</code> and save structures, Fourier maps, powder patterns and PDFs as reference HDF5/NumPy files. Set fixed random seeds. Set up pytest and CI (GitHub Actions on Linux and macOS).</dd>
    <dt>Deliverable</dt><dd><code>tests/golden/</code> plus a CI job that runs it on every push.</dd>
    <dt>Done when</dt><dd>All main command areas (read, symmetry, replace, domain, stack, fourier, powder, pdf, mmc, rmc) are covered by at least one test.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP1 · Build and packaging</h3>
  <p class="meta">Months 1–4</p>
  <dl>
    <dt>Tasks</dt><dd>Add <code>pyproject.toml</code> with scikit-build-core on top of the existing CMake. Add a <code>DIFFUSE_PYTHON_CORE</code> option that builds without X11, PGPLOT, readline or the <code>kuplot</code> plotting code. Build wheels with cibuildwheel for Linux (manylinux) and macOS, with FFTW and HDF5 bundled. Look at Windows support later.</dd>
    <dt>Deliverable</dt><dd><code>pip install .</code> and <code>import diffuse</code> both work, and CI builds wheels.</dd>
    <dt>Done when</dt><dd>A clean virtual environment can install the wheel and run one smoke test.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP2 · Core infrastructure</h3>
  <p class="meta">Months 2–6</p>
  <dl>
    <dt>Tasks</dt><dd>Add the API-mode flag, output redirection and logging callback. Map errors to exceptions. Remove or convert every <code>STOP</code>. Add <code>dc_init</code> and <code>dc_reset</code>. Expose the random seed and state (<code>random_mod</code>) so runs are reproducible. Add a <code>legacy.run_macro()</code> helper for running old macros during the migration.</dd>
    <dt>Deliverable</dt><dd><code>diffuse._core</code> with init, reset, errors, logging and RNG control.</dd>
    <dt>Done when</dt><dd>A deliberately bad input raises a Python exception. The process survives, prints nothing to the terminal, and the next call works.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP3 · Crystal data model and I/O</h3>
  <p class="meta">Months 3–8</p>
  <dl>
    <dt>Tasks</dt><dd>Python <code>Crystal</code> class (cell, space group, atoms, ADPs, occupancies, properties, molecules). Two-way transfer to and from <code>crystal_mod</code>. Readers and writers for CIF (gemmi), DISCUS <code>.stru</code>, xyz and HDF5. Space-group expansion through the existing <code>spcgr_*</code> routines, cross-checked with spglib.</dd>
    <dt>Deliverable</dt><dd><code>dc.Crystal</code>, <code>dc.io</code>.</dd>
    <dt>Done when</dt><dd>Round-tripping all golden structures (Python → Fortran → Python) gives identical data, and the files written match the legacy <code>save</code> output.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP4 · Single-crystal Fourier (first vertical slice)</h3>
  <p class="meta">Months 5–10</p>
  <dl>
    <dt>Tasks</dt><dd>Shim over <code>four_run</code> and <code>four_strucf</code>: radiation, wavelength, corners and increments, form factors, ADPs, lots, Friedel symmetry and output arrays. Results returned as <code>xarray</code> with hkl coordinates. Optional builds with OpenMP, CUDA and FINUFFT.</dd>
    <dt>Deliverable</dt><dd><code>dc.scattering.single_crystal()</code></dd>
    <dt>Done when</dt><dd>Results match the golden data to a set tolerance (for example a relative difference of 1e-10 in double precision), and speed is within 5% of the binary.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP5 · Powder</h3>
  <p class="meta">Months 8–13</p>
  <dl>
    <dt>Tasks</dt><dd>Shims over <code>powder_run</code>, <code>powder_complete</code> and <code>powder_nufft</code> (Debye, complete, FFT). Profile functions, preferred orientation, instrumental settings. Output as Q, 2θ and d. Comparison with experimental data in Python.</dd>
    <dt>Deliverable</dt><dd><code>dc.scattering.powder()</code></dd>
    <dt>Done when</dt><dd>All powder golden tests pass.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP6 · PDF, 3D-PDF and Patterson</h3>
  <p class="meta">Months 10–15</p>
  <dl>
    <dt>Tasks</dt><dd>Shims over <code>pdf_setup</code>, <code>pdf_run</code> and <code>pdf_determine</code>, plus the 3D-PDF (<code>discus_3dpdf</code>), <code>exp2pdf</code> and Patterson routines. Keep the PDF-fitting code paths (<code>pdf_rmc_*</code>) for WP8.</dd>
    <dt>Deliverable</dt><dd><code>dc.scattering.pdf()</code>, <code>pdf3d()</code>, <code>patterson()</code>, <code>exp2pdf()</code></dd>
    <dt>Done when</dt><dd>All PDF golden tests pass.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP7 · Structure building and defects</h3>
  <p class="meta">Months 9–18 · broad, so done in several releases</p>
  <dl>
    <dt>Tasks</dt><dd>Python functions over the existing routines for symmetry operations, replace and remove, insert, domains (<code>domain</code>), stacking faults (<code>stack</code>), surfaces, waves, shear, molecules (molecularize, place molecules, decorate), and connectivity and neighbour analysis (<code>chem</code>, <code>conn</code>). Priority goes to the commands that users rely on most.</dd>
    <dt>Deliverable</dt><dd><code>dc.defects.*</code>, <code>dc.analysis.*</code></dd>
    <dt>Done when</dt><dd>Each wrapped command has a golden test. Stacking faults and domains are written up as notebook tutorials.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP8 · Monte Carlo (MMC/RMC)</h3>
  <p class="meta">Months 14–22 · the most stateful part</p>
  <dl>
    <dt>Tasks</dt><dd>Wrap <code>mmc</code> (~9k lines) and <code>rmc</code>. Express targets (correlations, distances, energies) as Python objects and send them to the existing set-up code. Report progress through a callback so the caller can follow and stop long runs. Later option: Python-defined energy terms through callbacks.</dd>
    <dt>Deliverable</dt><dd><code>dc.mc.mmc()</code>, <code>dc.mc.rmc()</code></dd>
    <dt>Done when</dt><dd>With fixed seeds, results are statistically the same as the golden runs, and runs can be stopped and resumed.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP9 · Refinement and plotting replaced in Python</h3>
  <p class="meta">Months 6–14 · can run alongside WP4–WP7</p>
  <dl>
    <dt>Tasks</dt><dd>Least-squares refinement with lmfit (replaces <code>refine</code>). Differential evolution (replaces <code>diffev</code>), either scipy's or a ported DE that keeps the current DIFFEV options, with parallel workers through concurrent.futures, mpi4py or Dask. Plot helpers in matplotlib and plotly for 1D, 2D maps and 3D slices (replaces <code>kuplot</code> for typical uses). Readers for the formats <code>kuplot</code> loads.</dd>
    <dt>Deliverable</dt><dd><code>dc.refine</code>, <code>dc.plot</code></dd>
    <dt>Done when</dt><dd>Two or three published refinement examples are reproduced in notebooks.</dd>
  </dl>
</div>

<div class="wp">
  <h3>WP10 · Service layer and web front end (first version)</h3>
  <p class="meta">Months 16–24</p>
  <dl>
    <dt>Tasks</dt><dd>A FastAPI service with a job queue (Celery, RQ or Dask) and one worker process per job. Pydantic schemas for model and calculation requests. Results stored as HDF5/Zarr. A first web UI: load a structure, set up a single-crystal or powder calculation, view the results interactively (plotly; three.js for 3D).</dd>
    <dt>Deliverable</dt><dd>A deployable container image and a demo UI.</dd>
    <dt>Done when</dt><dd>A user can go from a CIF file to a diffuse-scattering map in a browser.</dd>
  </dl>
</div>

<h3 id="b-timeline">7.7 Timeline</h3>
<p>Months from the start, assuming about 1.5 full-time developers. Some packages overlap because they can be worked on in parallel.</p>
<div class="table-wrap">
<table>
<thead><tr><th style="width:38%">Work package</th><th>Months 0 → 24</th></tr></thead>
<tbody>
<tr><td>WP0 Regression corpus</td><td><span class="bar" style="margin-left:0%;width:12.5%"></span> 0–3</td></tr>
<tr><td>WP1 Build and packaging</td><td><span class="bar" style="margin-left:4.2%;width:12.5%"></span> 1–4</td></tr>
<tr><td>WP2 Core infrastructure</td><td><span class="bar" style="margin-left:8.3%;width:16.7%"></span> 2–6</td></tr>
<tr><td>WP3 Crystal model and I/O</td><td><span class="bar" style="margin-left:12.5%;width:20.8%"></span> 3–8</td></tr>
<tr><td>WP4 Single-crystal Fourier</td><td><span class="bar" style="margin-left:20.8%;width:20.8%"></span> 5–10</td></tr>
<tr><td>WP9 Refinement and plotting</td><td><span class="bar" style="margin-left:25%;width:33.3%"></span> 6–14</td></tr>
<tr><td>WP5 Powder</td><td><span class="bar" style="margin-left:33.3%;width:20.8%"></span> 8–13</td></tr>
<tr><td>WP7 Structure building and defects</td><td><span class="bar" style="margin-left:37.5%;width:37.5%"></span> 9–18</td></tr>
<tr><td>WP6 PDF, 3D-PDF, Patterson</td><td><span class="bar" style="margin-left:41.7%;width:20.8%"></span> 10–15</td></tr>
<tr><td>WP8 Monte Carlo</td><td><span class="bar" style="margin-left:58.3%;width:33.3%"></span> 14–22</td></tr>
<tr><td>WP10 Service and web UI</td><td><span class="bar" style="margin-left:66.7%;width:33.3%"></span> 16–24</td></tr>
</tbody>
</table>
</div>
<p><strong>Milestones:</strong></p>
<ul>
  <li><strong>M1 (month 6):</strong> installable package, error handling and logging, crystal round trip. Alpha release for internal users.</li>
  <li><strong>M2 (month 12):</strong> single-crystal Fourier, powder, Python refinement and plotting. First public beta, usable for most diffuse-scattering workflows from notebooks.</li>
  <li><strong>M3 (month 18):</strong> PDF and the main defect tools. Most of what users do no longer needs macros.</li>
  <li><strong>M4 (month 24):</strong> Monte Carlo and the first web front end. Python becomes the primary interface; <code>discus_suite</code> moves to maintenance only.</li>
</ul>

<h3 id="b-testing">7.8 Testing strategy</h3>
<ul>
  <li><strong>Three levels.</strong> (1) Golden tests against the legacy binary. (2) Unit tests for the shims (shapes, errors, reset). (3) Physics checks that don't depend on legacy output: Bragg positions, Friedel symmetry, the Debye result vs. the Fourier result for the same structure, PDF normalisation, sum rules.</li>
  <li><strong>Tolerances.</strong> Use tight tolerances for deterministic calculations. Use statistical comparisons (means and spreads over seeds) for MC/RMC.</li>
  <li><strong>CI matrix.</strong> gfortran on Linux and macOS for every push; ifort/ifx and the CUDA build nightly or on demand.</li>
  <li><strong>Performance tracking.</strong> Benchmark key kernels (asv or pytest-benchmark) so the wrapper layer never adds noticeable overhead.</li>
</ul>

<h3 id="b-parallel">7.9 Parallelism and performance</h3>
<ul>
  <li><strong>OpenMP</strong> continues to run inside the kernels. The nanobind layer releases the Python GIL during long calls, so notebooks stay responsive and threads can run other work.</li>
  <li><strong>MPI</strong>, used today by <code>diffev</code> and <code>refine</code> to run macros in parallel, is replaced by Python process pools (concurrent.futures, Dask, Ray or mpi4py on clusters). Each worker holds its own Fortran state, which fits the global-module design.</li>
  <li><strong>CUDA and FINUFFT</strong> become optional build variants, for example <code>diffuse[cuda]</code>. Python picks the backend at run time and falls back to CPU.</li>
  <li><strong>Data copies.</strong> Arrays cross the boundary as contiguous NumPy arrays (Fortran order where needed). Large outputs such as 3D Fourier volumes are written straight into buffers that Python allocates, so there is no extra copy.</li>
</ul>

<h3 id="b-risks">7.10 Risks and mitigations</h3>
<div class="table-wrap">
<table>
<thead><tr><th>Risk</th><th>Impact</th><th>Mitigation</th></tr></thead>
<tbody>
<tr><td>Hidden coupling: kernels depend on module state set by unrelated commands</td><td>Wrong results that look plausible</td><td>Golden tests per command. Start each shim call from <code>dc_reset()</code> and set every relevant variable explicitly. Write down which module variables each kernel reads.</td></tr>
<tr><td>Validation logic lives only in menu <code>set</code> handlers</td><td>The API accepts invalid input</td><td>Move validation into shared routines that both the menu and the shim call. Add pydantic validation in Python.</td></tr>
<tr><td>Packaging Fortran with HDF5 and FFTW as wheels</td><td>Installation problems, especially on macOS and Windows</td><td>cibuildwheel with pinned dependencies. Offer a conda-forge package as well, which is often easier for scientific users.</td></tr>
<tr><td>Development continues in the Fortran code during the migration</td><td>Merge conflicts, a moving target</td><td>Keep shims in a separate <code>capi/</code> tree. Golden tests catch behaviour changes. Agree that new features are added through the API first.</td></tr>
<tr><td>Users depend on macros</td><td>Slow adoption</td><td>Keep <code>discus_suite</code>. Provide <code>legacy.run_macro()</code>. Write side-by-side "macro → Python" guides, and possibly a translator for common commands.</td></tr>
<tr><td>Size of MMC/RMC</td><td>WP8 overruns</td><td>Schedule it last, wrap rather than refactor, and limit the first release to the most-used energy terms.</td></tr>
<tr><td>Limited developer time and few people with Fortran experience</td><td>The plan slips</td><td>Release in milestones so value arrives early. Make the Python layer welcoming to new contributors. The Fortran work is concentrated in the shim and infrastructure packages.</td></tr>
</tbody>
</table>
</div>

<h3 id="b-90">7.11 First 90 days</h3>
<ul class="checklist">
  <li>Agree on package name, licence, supported platforms and minimum Python version.</li>
  <li>Choose and document 30+ golden macros. Run them with the current release and store the outputs.</li>
  <li>Set up pytest and GitHub Actions, run the golden tests through <code>suite_python</code>.</li>
  <li>Add the <code>DIFFUSE_PYTHON_CORE</code> CMake option and cut the <code>discus</code> → <code>kuplot</code> plotting dependency for that build.</li>
  <li>Add <code>pyproject.toml</code> (scikit-build-core) and a minimal nanobind module that only calls <code>dc_init</code> and returns the version.</li>
  <li>Prototype the API-mode output redirect and error-to-exception mapping.</li>
  <li>Make an inventory of all <code>STOP</code> statements and interactive fallback paths, and a plan to remove each one.</li>
  <li>Prototype <code>dc_set_cell</code>, <code>dc_set_atoms</code> and <code>dc_get_atoms</code>, plus a Python <code>Crystal</code> that round-trips a golden structure.</li>
  <li>Write down which module variables <code>four_run</code> reads, as preparation for WP4.</li>
</ul>

<h3 id="b-decisions">7.12 Decisions needed</h3>
<ul>
  <li><strong>Binding technology:</strong> <code>bind(C)</code> + nanobind (recommended: cleaner arrays, exceptions and GIL control) or f2py (already used by <code>suite_python</code>, less code at first, weaker for complex interfaces).</li>
  <li><strong>Distribution:</strong> PyPI wheels, conda-forge, or both.</li>
  <li><strong>Supported platforms and compilers</strong> for official builds. Is Windows native, or WSL only?</li>
  <li><strong>Backward-compatibility promise</strong> for macros and file formats, and how long <code>discus_suite</code> is maintained.</li>
  <li><strong>Data formats:</strong> NeXus/HDF5 as the standard result format, and which version of the DISCUS structure format to keep.</li>
  <li><strong>Web goals:</strong> local single-user app, shared group server or public service. This decides how much work goes into authentication, quotas and job management in WP10.</li>
</ul>

<footer>Based on a read-only review of the DiffuseCode repository (version 6.23.01). Line counts are approximate.</footer>
</main>
</body>
</html>

"""Demo: read a CIF, calculate powder patterns with DISCUS, plot with matplotlib.

    DIFFUSE_PYTHON_PATH=<build>/python PYTHONPATH=python \\
        python3 python/examples/powder_demo.py [--show] [--out powder_demo.png]

The structure is anatase (TiO2), COD entry 1010942 (anatase.cif here).
"""

import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt

import diffuse as dc

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cif", default=HERE / "anatase.cif", type=Path)
    parser.add_argument("--out", default="powder_demo.png", type=Path)
    parser.add_argument("--show", action="store_true", help="open a plot window")
    args = parser.parse_args()

    # 1. Read the structure: DISCUS applies the space group and builds a supercell
    crystal = dc.read_structure(args.cif, supercell=(10, 10, 4))
    print(crystal)

    # 2. Powder patterns: X-rays and neutrons, Debye equation
    t0 = time.perf_counter()
    xray = dc.powder(crystal, qmin=0.5, qmax=10, dq=0.005, radiation="xray")
    neutron = dc.powder(crystal, qmin=0.5, qmax=10, dq=0.005, radiation="neutron")
    print(f"two patterns of {len(xray.q)} points in {time.perf_counter() - t0:.2f} s")

    # 3. Finite-size broadening: the same structure as a much smaller crystal
    small = dc.read_structure(args.cif, supercell=(3, 3, 1))
    xray_small = dc.powder(small, qmin=0.5, qmax=10, dq=0.005, radiation="xray")

    # 4. Plot: these are plain NumPy arrays, so any plotting library works
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                      constrained_layout=True)
    top.plot(xray.q, xray.intensity / xray.intensity.max(), lw=1, label="X-ray")
    top.plot(neutron.q, neutron.intensity / neutron.intensity.max(), lw=1, label="neutron")
    top.set_title(f"Anatase, {crystal.supercell[0]}×{crystal.supercell[1]}×"
                  f"{crystal.supercell[2]} unit cells ({crystal.natoms} atoms)")
    top.set_ylabel("Intensity (normalised)")
    top.legend()

    for p, label in ((xray, f"{crystal.natoms} atoms"), (xray_small, f"{small.natoms} atoms")):
        bottom.plot(p.q, p.intensity / p.intensity.max(), lw=1, label=label)
    bottom.set_title("X-ray: crystal size broadens the peaks")
    bottom.set_xlabel(r"$Q$ (Å$^{-1}$)")
    bottom.set_ylabel("Intensity (normalised)")
    bottom.set_xlim(xray.q[0], xray.q[-1])
    bottom.legend()

    fig.savefig(args.out, dpi=150)
    print(f"plot written to {args.out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

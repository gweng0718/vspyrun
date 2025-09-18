# vspyrun

Small, composable helpers for building and managing **VASP** inputs (INCAR, KPOINTS, POSCAR, POTCAR) and simple file utilities. This repo is organized with the modern **src/** layout and exposes a flat public API for convenience.

## Features

- **INCAR builders**: predefined blocks for SCF, bands, relax, SOC, vdW, geometry & parallel settings.
- **KPOINTS helpers**: uniform mesh for SCF, line-mode splitter (`K1`, `K2`, …) for band paths.
- **POSCAR utilities**: mirror/inversion builders, vacuum & centering, element filtering, selective dynamics, quick writers.
- **POTCAR**: assemble a combined POTCAR from a directory of element subfolders.
- **File helpers**: tiny utilities for listing/handling files in your run tree.

## Install (editable)

```bash
# from the repository root
python -m venv .venv
source .venv/bin/activate  # or: conda activate <env>
pip install -e .
```

> Requires Python 3.9+ and the packages: `numpy`, `ase`.

## Project layout

```
VSPYRUN/
├─ pyproject.toml
├─ README.md
└─ src/
   └─ vspyrun/
      ├─ __init__.py
      ├─ file.py
      ├─ incar.py
      ├─ kpoints.py
      ├─ poscar.py
      └─ potcat.py
```

## Quick start

```python
from vspyrun import (
    # INCAR
    incar_sys_scf, incar_sys_bs, incar_sys_rel,
    incar_elec, incar_soc, incar_ecorr, incar_geo, incar_para, incar_solv,
    # KPOINTS
    kpts_scf, kpts_bs_split,
    # POSCAR
    build_mir_sym, add_vacuum, center_atoms, write_pos_cart, write_pos_direct,
    # POTCAR
    generate_potcar,
    # File helpers
    get_file_name,
)

# 1) INCAR: compose a simple SCF input
incar_text = "\n".join([
    incar_sys_scf(),
    incar_elec(algo="Normal", amin="0.1", encut="520", prec="Accurate", sigma="0.1", ispin="1", xcf="PE"),
    incar_geo(ibrion="2", isif="2", ediffg="-0.02"),
])
open("INCAR", "w").write(incar_text)

# 2) KPOINTS: uniform mesh for SCF
open("KPOINTS", "w").write(kpts_scf(kcen="Gamma", nx="8", ny="8", nz="8"))

# 3) KPOINTS: split a line-mode file into segments K1, K2, ...
kpts_scf_text = open("KPOINTS", "r").read()  # your line-mode KPOINTS
kpts_bs_split(inkpts="KPOINTS", outdir="segments", nkpts=101)

# 4) POTCAR: assemble from element folders (e.g., POT_GGA_PAW_PBE/Ag/POTCAR, Sn/POTCAR, ...)
generate_potcar(inpos="POSCAR", potcar_dir="/path/to/paw/PBE", outpot="POTCAR")

# 5) POSCAR helpers: mirror, add vacuum, center, then write
# (Requires ASE Atoms, not shown here.)
# from ase.io import read
# atoms = read("POSCAR")
# atoms2 = build_mir_sym(atoms, z_mir=0.0, tol=0.1)
# atoms2 = add_vacuum(atoms2, vacuum_space=15.0)
# atoms2 = center_atoms(atoms2)
# write_pos_cart(atoms2, "POSCAR.new")
```

## API overview (selected)

- `incar.py`: `incar_sys_scf`, `incar_sys_bs`, `incar_sys_rel`, `incar_elec`, `incar_soc`, `incar_ecorr`, `incar_geo`, `incar_para`, `incar_solv`, plus small helpers for composing blocks.
- `kpoints.py`: `kpts_scf`, `kpts_bs_split`.
- `poscar.py`: `build_mir_sym`, `add_vacuum`, `center_atoms`, `write_pos_cart`, `write_pos_direct`, and other geometry utilities.
- `potcat.py`: `generate_potcar`.
- `file.py`: `get_file_name` (and any other simple file helpers included).

## Tips

- Keep your **scripts thin**: put reusable logic into functions in these modules and import them where needed.
- Prefer absolute imports in your own code (`from vspyrun import kpts_scf`) to avoid path confusion on HPC.
- Use version control for any changes to templates you commonly use (INCAR/KPOINTS).

## License

This project is intended for internal research workflows. Add a license of your choice (e.g., MIT) if you plan to publish or share.

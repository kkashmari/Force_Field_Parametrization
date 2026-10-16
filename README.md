# Force-field parameterization tool

This package combines the corrected `.frc` builder with automatic native
GAFF2 and UFF pipelines for MOL2 and POSCAR inputs.

Author: Khatereh Kashmari  
Date: 2026-08-10  
Version: 8.0.0

### Required packages

- GUI: Tkinter (included in the `ffbuild` Miniforge environment).
- GAFF2: Open Babel, AmberTools, ParmEd, NumPy, and AMBER2LAMMPS.
- UFF: Open Babel and `lammps-interface`.
- Existing `.frc`: Python standard library; Open Babel/AmberTools/RDKit are
  optional charge engines.

## Interactive use

Run without arguments and answer the questions (you can drag files from
Finder into Terminal when a path is requested):

```bash
python3 parameterize.py
```

The command-line arguments described below remain available for batch jobs.

## Desktop GUI

On macOS, double-click `launch_gui.command`, or run:

```bash
/usr/local/Caskroom/miniforge/base/envs/ffbuild/bin/python parameterize_gui.py
```

Choose the input, force field, and output folder, then click **Run
parameterization**. GAFF2-specific fields and existing `.frc` fields appear
only when those workflows are selected. For `.frc`, atom typing and `.nta`
generation run first automatically; that new `.nta` is then used to build the
LAMMPS data file. To keep manually checked atom types, select **Use an existing
.nta file (manual atom types)** and choose a LAMMPS `.data` input plus its
matching `.nta`. Output is shown live in the window.

## Existing `.frc` force fields

```bash
python parameterize.py molecule.data --ff /path/to/pcff.frc \
  --nta molecule.nta -o molecule_pcff
```

Without `--nta`, atom types are assigned automatically. With `--nta`, the
selected types are preserved and the completed data file is written inside the
requested output directory.

The builder reports exact, equivalence, wildcard/similarity fallback, and
missing counts. Missing terms remain visible in the report and must not be
called successful parameterization.

## GAFF2

Dependencies: Open Babel, AmberTools (`antechamber`, `parmchk2`, `tleap`),
ParmEd/NumPy, and AMBER2LAMMPS.

```bash
python parameterize.py molecule.mol2 --ff gaff2 -o molecule_gaff2 --charge 0 \
  --amber-converter /path/to/AMBER2LAMMPS/amber_to_lammps.py
```

The output directory contains `.data`, `.params`, a runnable `in.*` file,
AMBER intermediates, `pipeline.log`, and `audit.json`. The pipeline stops if
`parmchk2` marks a term `ATTN/needs revision` or `tleap` reports missing terms.
The script automatically searches the common Miniforge `ffbuild` environment
locations, so activation is optional once the dependencies are installed.

## UFF

```bash
python3 -m pip install lammps-interface
python parameterize.py structure.mol2 --ff uff -o structure_uff
```

UFF conversion uses `lammps-interface`, which is intended mainly for P1
periodic crystals and MOFs. MOL2/POSCAR bonds are inferred from geometry and
must be inspected. The pipeline stops if it detects untyped/missing terms.

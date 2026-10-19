# Force-field parameterization tool

This package combines a native LigParGen OPLS-AA importer, the corrected
`.frc` builder, and automatic GAFF2 and UFF pipelines.

Author: Khatereh Kashmari  
Date: 2026-08-10  
Version: 8.0.0

### Required packages

- GUI: Tkinter (included in the `ffbuild` Miniforge environment).
- LigParGen import and web generation: Python standard library only; internet
  access is required when generating directly from SMILES.
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

## LigParGen OPLS-AA

The program can submit a SMILES string to LigParGen, download its native
results ZIP, select the LAMMPS file, validate it, and write the final OPLS-AA
data file:

```bash
python parameterize.py --smiles 'CCO' --ff ligpargen --name ethanol \
  --charge 0 --ligpargen-charge-model CM1A-LBCC --ligpargen-opt 0 \
  -o ethanol_opls
```

`CM1A-LBCC` is restricted to neutral molecules. For a charged molecule use,
for example, `--charge -1 --ligpargen-charge-model CM1A`. Optimization can be
`0`, `1`, `2`, or `3`. The default runner is LigParGen's public web service.
The submitted SMILES is sent to that service and its job/rate limits apply.

To use a local LigParGen Docker image built with your licensed BOSS install,
add `--ligpargen-runner docker`. The default image is `ligpargen:dev`; override
it with `--ligpargen-docker-image IMAGE`.

You can still generate the molecule manually on LigParGen and import its
downloaded LAMMPS data file or complete ZIP bundle:

```bash
python parameterize.py molecule.lmp --ff ligpargen -o molecule_opls
```

This workflow does not use a Tinker `.frc`, retype atoms, recalculate charges,
or alter coefficients. It validates the native LigParGen/BOSS LAMMPS topology
and coefficient tables, preserves the data, and creates a `run 0` input
with the required OPLS styles. Outputs are the copied `.data`,
`in.*_ligpargen_opls`, and `audit.json`. ZIP input is supported when it
contains exactly one recognizable LigParGen LAMMPS data file.

Some legacy LigParGen 2.1 files contain an incorrect header count for improper
interactions. When detected, the importer corrects only that count to match the
existing `Impropers` rows and records the repair in `audit.json`; it never
creates or changes a force-field term.

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

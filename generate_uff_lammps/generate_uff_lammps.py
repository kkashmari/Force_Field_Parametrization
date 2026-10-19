#!/usr/bin/env python3
# Author: Khatereh Kashmari
"""Convert a POSCAR structure to UFF LAMMPS data and an equilibration input."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


# Defaults for running without command-line arguments.

INPUT_DIRECTORY = Path(
    "/Users/khatereh/Desktop/Research/Polymer_CoScientist/generate_uff_lammps"
)
INPUT_POSCAR_NAME = "amor_model.vasp"
OUTPUT_DIRECTORY = INPUT_DIRECTORY / "amor_model_UFF_generated"
SYSTEM_NAME = "amor_model"


class GenerationError(RuntimeError):
    pass


def clean_path(value: str) -> str:
    return value.strip().strip("'\"").replace("\\ ", " ")


def find_command(name: str) -> str:
    found = shutil.which(name)
    if found:
        return found
    candidates = [
        Path("/usr/local/Caskroom/miniforge/base/envs/ffbuild/bin") / name,
        Path("/opt/homebrew/Caskroom/miniforge/base/envs/ffbuild/bin") / name,
        Path.home() / "miniforge3/envs/ffbuild/bin" / name,
        Path.home() / "mambaforge/envs/ffbuild/bin" / name,
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise GenerationError(
        f"Required command '{name}' was not found. Install Open Babel and "
        "lammps-interface or activate the ffbuild environment."
    )


def run(command: list[str], cwd: Path, log: Path) -> str:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    record = f"$ {' '.join(command)}\n{result.stdout}{result.stderr}"
    with log.open("a", encoding="utf-8") as handle:
        handle.write(record)
        if not record.endswith("\n"):
            handle.write("\n")
    if result.returncode:
        raise GenerationError(
            f"Command failed ({result.returncode}): {' '.join(command)}\nSee {log}"
        )
    return result.stdout + result.stderr


def validate_poscar(path: Path) -> None:
    if not path.is_file():
        raise GenerationError(f"POSCAR file not found: {path}")
    if path.suffix.lower() not in {
        ".vasp",
        ".poscar",
        "",
    } and path.name.upper() not in {"POSCAR", "CONTCAR"}:
        raise GenerationError("Input must be POSCAR, CONTCAR, .vasp, or .poscar")


def make_p1_cif(source: Path, destination: Path, log: Path) -> None:
    run(
        [find_command("obabel"), "-iposcar", str(source), "-O", str(destination), "-h"],
        destination.parent,
        log,
    )
    if not destination.is_file() or destination.stat().st_size == 0:
        raise GenerationError(f"Open Babel did not create {destination}")
    text = destination.read_text(encoding="utf-8", errors="replace")
    if not all(
        key in text for key in ("_cell_length_a", "_cell_length_b", "_cell_length_c")
    ):
        raise GenerationError("Open Babel did not preserve the POSCAR lattice")


def lammps_input(system_id: str, data_name: str) -> str:
    stages = [
        ("nvt", 1000, None, 50000),
        ("nvt", 300, None, 50000),
        ("npt", 300, 986.92326, 50000),
        ("nvt", 1000, None, 50000),
        ("nvt", 300, None, 100000),
        ("npt", 300, 29607.6978, 50000),
        ("nvt", 1000, None, 50000),
        ("nvt", 300, None, 100000),
        ("npt", 300, 49346.163, 50000),
        ("nvt", 1000, None, 50000),
        ("nvt", 300, None, 100000),
        ("npt", 300, 24673.0815, 5000),
        ("nvt", 1000, None, 5000),
        ("nvt", 300, None, 10000),
        ("npt", 300, 4934.6163, 5000),
        ("nvt", 1000, None, 5000),
        ("nvt", 300, None, 10000),
        ("npt", 300, 493.46163, 5000),
        ("nvt", 1000, None, 5000),
        ("nvt", 300, None, 10000),
        ("npt", 300, 1.0, 800000),
    ]
    protocol: list[str] = []
    for number, (ensemble, temperature, pressure, steps) in enumerate(stages, 1):
        if ensemble == "nvt":
            command = (
                f"fix step{number} all nvt temp {temperature} {temperature} $(100.0*dt)"
            )
        else:
            command = (
                f"fix step{number} all npt temp {temperature} {temperature} "
                f"$(100.0*dt) iso {pressure:g} {pressure:g} $(100.0*dt)"
            )
        protocol.extend((command, f"run {steps}", f"unfix step{number}", ""))

    header = f"""# UFF relaxation, 21-stage equilibration, and density profiles
# Keep this input beside {data_name}; run: lmp -in in.{system_id}


#------------Initialization------------

units           real
dimension       3
boundary        p p p
newton          on
variable myid string {system_id}

#------------Force Field------------------------------------

atom_style full
pair_style lj/cut 12.5
bond_style harmonic
angle_style fourier
dihedral_style harmonic
special_bonds lj/coul 0.0 0.0 1.0
pair_modify tail yes mix arithmetic

neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes one 5000 page 100000

#------------Read DataFile---------------

read_data ${{input_dir}}/{data_name}

#-------------------------------------Settings----------------------------------------------

variable ptemp equal 300
timestep 1.0
velocity all create ${{ptemp}} 7667654 dist gaussian
thermo 1000
thermo_style custom step temp press density ke pe etotal evdwl epair ebond eangle edihed eimp
log ${{output_dir}}/${{myid}}.log.lammps


min_style cg
minimize 1.0e-4 1.0e-6 1000 100000
fix overlap_relax all nve/limit 0.01
run 10000
unfix overlap_relax

#------------------Write DataFile------------------

write_data ${{output_dir}}/${{myid}}_minimized.data

#--------------------------------------------Equilibration------------------------------------------------

#21 stages total 1,560,000 steps (1.56 ns)

reset_timestep 0
log ${{output_dir}}/${{myid}}_Densified.log.lammps
dump equil_dump all custom 10000 ${{output_dir}}/${{myid}}_equil.lammpstrj id mol type q xs ys zs ix iy iz
restart 100000 ${{output_dir}}/${{myid}}.restart.1 ${{output_dir}}/${{myid}}.restart.2

"""
    footer = """undump equil_dump

#---------------------------------------Density-profile sampling------------------------------------------

compute chunk_x all chunk/atom bin/1d x lower 0.02 units reduced
compute chunk_y all chunk/atom bin/1d y lower 0.02 units reduced
compute chunk_z all chunk/atom bin/1d z lower 0.02 units reduced
fix profile_x all ave/chunk 100 100 10000 chunk_x density/mass file ${output_dir}/${myid}_densityprofile_x.txt
fix profile_y all ave/chunk 100 100 10000 chunk_y density/mass file ${output_dir}/${myid}_densityprofile_y.txt
fix profile_z all ave/chunk 100 100 10000 chunk_z density/mass file ${output_dir}/${myid}_densityprofile_z.txt
fix profile_nvt all nvt temp 300 300 $(100.0*dt)

run 100000

unfix profile_nvt
unfix profile_x
unfix profile_y
unfix profile_z

#------------Write DataFile----------------------

write_data ${output_dir}/${myid}_Densified.data

"""
    return header + "\n".join(protocol) + footer


def generate(source: Path, output: Path, name: str) -> dict:
    validate_poscar(source)
    if output.exists() and any(output.iterdir()):
        raise GenerationError(f"Output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    log = output / "pipeline.log"
    cif = output / f"{name}_P1.cif"
    make_p1_cif(source, cif, log)

    before = {item.resolve() for item in output.iterdir()}
    report = run([find_command("lammps-interface"), cif.name], output, log)
    created = [
        item
        for item in output.iterdir()
        if item.is_file() and item.resolve() not in before
    ]
    data_files = [
        item
        for item in created
        if item.suffix.lower() in {".data", ".lmp"} or "data" in item.name.lower()
    ]
    if len(data_files) != 1:
        detected = ", ".join(item.name for item in data_files) or "none"
        raise GenerationError(f"Expected one LAMMPS data file; detected: {detected}")
    problems = [
        line.strip()
        for line in report.splitlines()
        if re.search(r"missing|untyped|could not|error", line, re.I)
    ]
    if problems:
        raise GenerationError(
            "UFF typing/parameter problem:\n  " + "\n  ".join(problems[:20])
        )

    for item in created:
        if item.name.startswith("in.") or item.suffix.lower() == ".in":
            item.unlink()
    input_file = output / f"in.{name}_P1"
    input_file.write_text(
        lammps_input(f"{name}_P1", data_files[0].name), encoding="utf-8"
    )
    audit = {
        "forcefield": "UFF",
        "status": "complete",
        "source": str(source),
        "lammps_data": str(data_files[0]),
        "lammps_input": str(input_file),
        "pipeline_log": str(log),
        "warning": "Bonds and UFF types are inferred from geometry and must be inspected.",
    }
    (output / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate UFF LAMMPS data and a ready-to-run input from POSCAR."
    )
    parser.add_argument("input", nargs="?", help="POSCAR/CONTCAR/.vasp file")
    parser.add_argument("-o", "--output", help="new or empty output directory")
    parser.add_argument("--name", help="output basename; default is input stem")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.input:
        source = Path(clean_path(args.input)).expanduser().resolve()
    else:
        source = (INPUT_DIRECTORY / INPUT_POSCAR_NAME).expanduser().resolve()
    name = args.name or (
        SYSTEM_NAME
        if not args.input
        else (source.stem if source.suffix else source.name.lower())
    )
    if args.output:
        output = Path(clean_path(args.output)).expanduser().resolve()
    else:
        output = OUTPUT_DIRECTORY.expanduser().resolve()
    try:
        audit = generate(source, output, name)
    except (OSError, GenerationError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(audit, indent=2))
    print("\nRun directly:")
    print(f"  cd {output}")
    print(f"  lmp -in {Path(audit['lammps_input']).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

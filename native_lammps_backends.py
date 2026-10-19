#!/usr/bin/env python3
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0
# Required packages: Open Babel; AmberTools; ParmEd; NumPy; AMBER2LAMMPS; lammps-interface
"""Native GAFF2 and UFF pipelines that produce LAMMPS input artifacts."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


class PipelineError(RuntimeError):
    pass


def find_command(name: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    # Desktop apps and newly opened shells often do not inherit `conda
    # activate`. Search the known environment explicitly so the wizard works
    # without requiring activation in the parent process.
    candidates: list[Path] = []
    for variable in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX"):
        prefix = os.environ.get(variable)
        if prefix:
            candidates.extend([
                Path(prefix) / "bin" / name,
                Path(prefix) / "envs" / "ffbuild" / "bin" / name,
            ])
    candidates.extend([
        Path("/usr/local/Caskroom/miniforge/base/envs/ffbuild/bin") / name,
        Path("/opt/homebrew/Caskroom/miniforge/base/envs/ffbuild/bin") / name,
        Path.home() / "miniforge3" / "envs" / "ffbuild" / "bin" / name,
        Path.home() / "mambaforge" / "envs" / "ffbuild" / "bin" / name,
    ])
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise PipelineError(
        f"Required command was not found: {name}. Activate ffbuild or install it there."
    )


def run(command: list[str], cwd: Path, log: Path) -> str:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    combined = f"$ {' '.join(command)}\n{result.stdout}{result.stderr}"
    with log.open("a", encoding="utf-8") as handle:
        handle.write(combined)
        if not combined.endswith("\n"):
            handle.write("\n")
    if result.returncode:
        raise PipelineError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"See {log}"
        )
    return result.stdout + result.stderr


def input_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".mol2":
        return "mol2"
    if path.name.upper() in {"POSCAR", "CONTCAR"} or suffix in {".vasp", ".poscar"}:
        return "poscar"
    raise PipelineError("Input must be MOL2, POSCAR, CONTCAR, .vasp, or .poscar")


def make_dir(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise PipelineError(f"Output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def convert_with_obabel(source: Path, destination: Path, fmt: str, log: Path) -> None:
    obabel = find_command("obabel")
    in_flag = "-imol2" if fmt == "mol2" else "-iposcar"
    run([obabel, in_flag, str(source), "-O", str(destination), "-h"], destination.parent, log)
    if not destination.is_file() or destination.stat().st_size == 0:
        raise PipelineError(f"Open Babel did not create {destination}")


def make_uff_cif(source: Path, destination: Path, fmt: str, log: Path,
                 padding: float = 14.0) -> None:
    """Create a P1 CIF; cell-less MOL2 files receive an orthogonal box."""
    if fmt == "poscar":
        convert_with_obabel(source, destination, fmt, log)
        text = destination.read_text(encoding="utf-8", errors="replace")
        required = ("_cell_length_a", "_cell_length_b", "_cell_length_c")
        if not all(key in text for key in required):
            raise PipelineError("Open Babel did not preserve the POSCAR lattice in CIF output")
        return

    xyz = destination.with_suffix(".xyz")
    obabel = find_command("obabel")
    run([obabel, "-imol2", str(source), "-oxyz", "-O", str(xyz), "-h"],
        destination.parent, log)
    lines = xyz.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        atom_count = int(lines[0].strip())
        records = []
        for line in lines[2:2 + atom_count]:
            fields = line.split()
            records.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    except (IndexError, ValueError) as exc:
        raise PipelineError(f"Could not parse Open Babel XYZ output: {xyz}") from exc
    if len(records) != atom_count or not records:
        raise PipelineError(f"Expected {atom_count} XYZ atoms but parsed {len(records)}")

    xs = [r[1] for r in records]
    ys = [r[2] for r in records]
    zs = [r[3] for r in records]
    minima = (min(xs) - padding, min(ys) - padding, min(zs) - padding)
    lengths = (
        max(xs) - min(xs) + 2.0 * padding,
        max(ys) - min(ys) + 2.0 * padding,
        max(zs) - min(zs) + 2.0 * padding,
    )
    cif_lines = [
        "data_generated_P1",
        "_symmetry_space_group_name_H-M 'P 1'",
        "_symmetry_Int_Tables_number 1",
        f"_cell_length_a {lengths[0]:.10f}",
        f"_cell_length_b {lengths[1]:.10f}",
        f"_cell_length_c {lengths[2]:.10f}",
        "_cell_angle_alpha 90.0", "_cell_angle_beta 90.0", "_cell_angle_gamma 90.0",
        "loop_", "_atom_site_label", "_atom_site_type_symbol",
        "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    counts: dict[str, int] = {}
    for symbol, x, y, z in records:
        counts[symbol] = counts.get(symbol, 0) + 1
        label = f"{symbol}{counts[symbol]}"
        fx = (x - minima[0]) / lengths[0]
        fy = (y - minima[1]) / lengths[1]
        fz = (z - minima[2]) / lengths[2]
        cif_lines.append(f"{label:<8} {symbol:<3} {fx:.10f} {fy:.10f} {fz:.10f} 1.0")
    destination.write_text("\n".join(cif_lines) + "\n", encoding="utf-8")


def locate_amber_converter(explicit: str | None) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    command = shutil.which("amber_to_lammps.py")
    if command:
        candidates.append(Path(command))
    candidates.extend([
        Path.cwd() / "amber_to_lammps.py",
        Path(__file__).resolve().parent / "AMBER2LAMMPS" / "amber_to_lammps.py",
        Path.home() / "AMBER2LAMMPS" / "amber_to_lammps.py",
        Path.home() / "amber2lammps" / "amber_to_lammps.py",
    ])
    for item in candidates:
        if item.is_file():
            return item.resolve()
    raise PipelineError(
        "amber_to_lammps.py was not found. Pass --amber-converter PATH or clone "
        "https://github.com/askforarun/AMBER2LAMMPS.git into this package."
    )


def gaff2(source: Path, out_dir: Path, name: str, net_charge: int,
          converter: str | None, buffer: float) -> dict:
    fmt = input_format(source)
    make_dir(out_dir)
    log = out_dir / "pipeline.log"
    antechamber = find_command("antechamber")
    parmchk2 = find_command("parmchk2")
    tleap = find_command("tleap")
    converter_path = locate_amber_converter(converter)

    if fmt == "mol2":
        # Going MOL2 -> PDB -> antechamber destroys explicit bond orders.  In
        # particular, PDB CONECT can represent a double bond by repeating it,
        # which antechamber may then preserve as duplicate single bonds.
        amber_input = source
        amber_input_format = "mol2"
    else:
        prepared_pdb = out_dir / f"{name}_input.pdb"
        convert_with_obabel(source, prepared_pdb, fmt, log)
        amber_input = prepared_pdb
        amber_input_format = "pdb"
    mol2 = out_dir / f"{name}_gaff2.mol2"
    coordinate_pdb = out_dir / f"{name}_gaff2.pdb"
    frcmod = out_dir / f"{name}_gaff2.frcmod"
    prmtop = out_dir / f"{name}_gaff2.prmtop"
    inpcrd = out_dir / f"{name}_gaff2.inpcrd"
    data = out_dir / f"{name}_gaff2.data"
    params = out_dir / f"{name}_gaff2.params"

    run([
        antechamber, "-i", str(amber_input), "-fi", amber_input_format,
        "-o", str(mol2),
        "-fo", "mol2", "-at", "gaff2", "-c", "bcc", "-nc", str(net_charge),
        "-s", "2", "-j", "4", "-dr", "yes",
    ], out_dir, log)
    # The current AMBER2LAMMPS interface consumes PDB coordinates separately
    # from the PRMTOP. Converting the newly typed MOL2 preserves its atom order.
    convert_with_obabel(mol2, coordinate_pdb, "mol2", log)
    run([parmchk2, "-i", str(mol2), "-f", "mol2", "-o", str(frcmod),
         "-s", "gaff2"], out_dir, log)

    frc_text = frcmod.read_text(encoding="utf-8", errors="replace")
    suspect = [line.strip() for line in frc_text.splitlines()
               if re.search(r"ATTN|needs revision", line, re.I)]
    if suspect:
        raise PipelineError(
            "parmchk2 produced parameters requiring manual review:\n  "
            + "\n  ".join(suspect[:20])
        )

    leap = out_dir / "tleap.in"
    leap.write_text(
        "source leaprc.gaff2\n"
        f"MOL = loadmol2 {mol2.name}\n"
        f"loadamberparams {frcmod.name}\n"
        "check MOL\n"
        f"saveamberparm MOL {prmtop.name} {inpcrd.name}\n"
        "quit\n", encoding="utf-8")
    leap_output = run([tleap, "-f", leap.name], out_dir, log)
    fatal = [line.strip() for line in leap_output.splitlines()
             if re.search(r"FATAL|Errors? = [1-9]|Could not find", line, re.I)]
    if fatal:
        raise PipelineError("tleap reported incomplete parameters:\n  " + "\n  ".join(fatal[:20]))

    amber_python = str(Path(antechamber).parent / "python")
    if not Path(amber_python).is_file():
        amber_python = sys.executable
    help_result = subprocess.run(
        [amber_python, str(converter_path), "--help"], text=True,
        capture_output=True)
    converter_help = help_result.stdout + help_result.stderr
    if "--topologies" in converter_help and "--charges" in converter_help:
        converter_command = [
            amber_python, str(converter_path), str(data), str(params),
            str(coordinate_pdb), "-t", str(prmtop), "-c", "1", "--charges",
            str(net_charge), "--verbose", "-b", str(buffer),
        ]
    else:
        # Compatibility with the original single-topology CLI documented by
        # LAMMPS releases before the converter added mixture support.
        converter_command = [
            amber_python, str(converter_path), str(data), str(params),
            str(prmtop), str(inpcrd), "--charge", str(net_charge),
            "--verbose", "-b", str(buffer),
        ]
    run(converter_command, out_dir, log)
    improper_helper = Path(__file__).resolve().parent / "add_amber_impropers.py"
    run([amber_python, str(improper_helper), str(data), str(params), str(prmtop)],
        out_dir, log)
    for required in (data, params):
        if not required.is_file() or required.stat().st_size == 0:
            raise PipelineError(f"Expected output was not created: {required}")

    lmp_input = out_dir / f"in.{name}_gaff2"
    lmp_input.write_text(
        "units real\natom_style full\nboundary p p p\n"
        f"read_data {data.name}\n"
        "pair_style lj/cut/coul/long 9.0 9.0\n"
        "bond_style harmonic\nangle_style harmonic\n"
        "dihedral_style fourier\n"
        "improper_style cvff\n"
        "special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333\n"
        f"include {params.name}\n"
        "kspace_style pppm 1.0e-5\nthermo 1\nrun 0\n", encoding="utf-8")
    audit = {
        "forcefield": "GAFF2", "status": "complete", "input": str(source),
        "net_charge": net_charge, "charge_method": "AM1-BCC",
        "data_file": str(data), "parameter_file": str(params),
        "lammps_input": str(lmp_input), "pipeline_log": str(log),
        "warnings": (["POSCAR bonds were inferred by Open Babel from geometry; inspect them."]
                     if fmt == "poscar" else []),
    }
    write_json(out_dir / "audit.json", audit)
    return audit


def uff_lammps_input(name: str, data_name: str) -> str:
    """Return a directly runnable UFF optimization/equilibration input file."""
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
    stage_lines = []
    for index, (ensemble, temperature, pressure, steps) in enumerate(stages, 1):
        if ensemble == "nvt":
            fix = f"fix step{index} all nvt temp {temperature} {temperature} $(100.0*dt)"
        else:
            fix = (f"fix step{index} all npt temp {temperature} {temperature} "
                   f"$(100.0*dt) iso {pressure:g} {pressure:g} $(100.0*dt)")
        stage_lines.extend([fix, f"run {steps}", f"unfix step{index}", ""])

    return (
        "# UFF structural relaxation, 21-stage equilibration, and density profiles\n"
        "# Run from this directory with: lmp -in " + f"in.{name}_P1\n\n"
        "units real\n"
        "dimension 3\n"
        "boundary p p p\n"
        "newton on\n"
        f"variable myid string {name}_P1\n"
        "variable input_dir string .\n"
        "variable output_dir string .\n\n"
        "atom_style full\n"
        "pair_style lj/cut 12.5\n"
        "bond_style harmonic\n"
        "angle_style fourier\n"
        "dihedral_style harmonic\n"
        "special_bonds lj/coul 0.0 0.0 1.0\n"
        "pair_modify tail yes mix arithmetic\n"
        "neighbor 2.0 bin\n"
        "neigh_modify delay 0 every 1 check yes one 5000 page 100000\n\n"
        f"read_data ${{input_dir}}/{data_name}\n\n"
        "variable ptemp equal 300\n"
        "timestep 1.0\n"
        "velocity all create ${ptemp} 7667654 dist gaussian\n"
        "thermo 1000\n"
        "thermo_style custom step temp press density ke pe etotal evdwl epair ebond eangle edihed eimp\n"
        "log ${output_dir}/${myid}.log.lammps\n\n"
        "# Remove severe local strain before finite-temperature equilibration.\n"
        "min_style cg\n"
        "minimize 1.0e-4 1.0e-6 1000 100000\n"
        "fix overlap_relax all nve/limit 0.01\n"
        "run 10000\n"
        "unfix overlap_relax\n"
        "write_data ${output_dir}/${myid}_minimized.data\n\n"
        "# The following 21 stages total 1,560,000 steps (1.56 ns).\n"
        "reset_timestep 0\n"
        "log ${output_dir}/${myid}_Densified.log.lammps\n"
        "dump equil_dump all custom 10000 ${output_dir}/${myid}_equil.lammpstrj id mol type q xs ys zs ix iy iz\n"
        "restart 100000 ${output_dir}/${myid}.restart.1 ${output_dir}/${myid}.restart.2\n\n"
        + "\n".join(stage_lines)
        + "undump equil_dump\n\n"
        "# Sample density profiles for 100 ps at fixed equilibrated volume.\n"
        "compute chunk_x all chunk/atom bin/1d x lower 1.0 units box\n"
        "compute chunk_y all chunk/atom bin/1d y lower 1.0 units box\n"
        "compute chunk_z all chunk/atom bin/1d z lower 1.0 units box\n"
        "fix profile_x all ave/chunk 100 100 10000 chunk_x density/mass file ${output_dir}/${myid}_densityprofile_x.txt\n"
        "fix profile_y all ave/chunk 100 100 10000 chunk_y density/mass file ${output_dir}/${myid}_densityprofile_y.txt\n"
        "fix profile_z all ave/chunk 100 100 10000 chunk_z density/mass file ${output_dir}/${myid}_densityprofile_z.txt\n"
        "fix profile_nvt all nvt temp 300 300 $(100.0*dt)\n"
        "run 100000\n"
        "unfix profile_nvt\n"
        "unfix profile_x\n"
        "unfix profile_y\n"
        "unfix profile_z\n\n"
        "write_data ${output_dir}/${myid}_Densified.data\n"
    )


def uff(source: Path, out_dir: Path, name: str) -> dict:
    fmt = input_format(source)
    make_dir(out_dir)
    log = out_dir / "pipeline.log"
    interface = find_command("lammps-interface")
    cif = out_dir / f"{name}_P1.cif"
    make_uff_cif(source, cif, fmt, log)

    before = {p.resolve() for p in out_dir.iterdir()}
    output = run([interface, cif.name], out_dir, log)
    created = [p for p in out_dir.iterdir() if p.resolve() not in before and p.is_file()]
    data_files = [p for p in created if p.suffix.lower() in {".data", ".lmp"}
                  or "data" in p.name.lower()]
    input_files = [p for p in created if p.name.startswith("in.") or p.suffix == ".in"]
    if not data_files:
        raise PipelineError(
            "lammps-interface finished but no LAMMPS data file was detected. "
            f"See {log}. Created: {', '.join(p.name for p in created) or 'none'}"
        )
    bad = [line.strip() for line in output.splitlines()
           if re.search(r"missing|untyped|could not|error", line, re.I)]
    if bad:
        raise PipelineError("UFF generator reported a typing/parameter problem:\n  "
                            + "\n  ".join(bad[:20]))
    if len(data_files) != 1:
        raise PipelineError(
            "Expected exactly one UFF LAMMPS data file, but detected: "
            + ", ".join(p.name for p in data_files)
        )
    # lammps-interface may emit a demonstration input with version-dependent
    # names and incomplete run settings. Replace it with one input whose
    # read_data target is the exact generated data basename.
    for old_input in input_files:
        old_input.unlink()
    lmp_input = out_dir / f"in.{name}_P1"
    lmp_input.write_text(
        uff_lammps_input(name, data_files[0].name), encoding="utf-8")
    input_files = [lmp_input]
    created = [p for p in created if p.exists()]
    if lmp_input not in created:
        created.append(lmp_input)
    audit = {
        "forcefield": "UFF", "status": "complete", "input": str(source),
        "data_files": [str(p) for p in data_files],
        "lammps_inputs": [str(p) for p in input_files],
        "all_created_files": [str(p) for p in created], "pipeline_log": str(log),
        "warnings": [
            "UFF topology is inferred from geometry; inspect bonds and periodic connectivity.",
            "lammps-interface is primarily designed for P1 periodic crystal/MOF structures.",
        ],
    }
    write_json(out_dir / "audit.json", audit)
    return audit

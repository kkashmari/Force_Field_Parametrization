#!/usr/bin/env python3
# Author: Khatereh Kashmari
"""Prepare LAMMPS inputs with LigParGen, a .frc file, GAFF2, or UFF."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from native_lammps_backends import PipelineError, gaff2, uff
from frc_pipeline import run_frc_pipeline
from ligpargen_backend import LigParGenError, import_ligpargen, parameterize_smiles


def clean_path(value: str) -> str:
    """Accept ordinary paths and paths dragged into macOS Terminal."""
    value = value.strip().strip("'\"")
    return value.replace("\\ ", " ")


def ask_existing_file(label: str) -> str:
    while True:
        value = clean_path(input(f"{label}: "))
        if Path(value).expanduser().is_file():
            return value
        print(f"File not found: {value}")


def ask_forcefield() -> str:
    print("\nForce field")
    print("  1. LigParGen OPLS-AA output (native BOSS parameters)")
    print("  2. GAFF2  (recommended for organic molecules/polymers)")
    print("  3. UFF    (especially crystals/MOFs/unusual elements)")
    print("  4. Existing .frc file")
    while True:
        choice = input("Select [1-4]: ").strip().lower()
        if choice in {"1", "ligpargen", "opls"}:
            return "ligpargen"
        if choice in {"2", "gaff2"}:
            return "gaff2"
        if choice in {"3", "uff"}:
            return "uff"
        if choice in {"4", "frc"}:
            return ask_existing_file("Path to .frc file")
        print("Please enter 1, 2, 3, or 4.")


def ask_int(label: str, default: int) -> int:
    while True:
        value = input(f"{label} [{default}]: ").strip()
        if not value:
            return default
        try:
            return int(value)
        except ValueError:
            print("Enter an integer, for example 0, 1, or -1.")


def ask_float(label: str, default: float) -> float:
    while True:
        value = input(f"{label} [{default}]: ").strip()
        if not value:
            return default
        try:
            return float(value)
        except ValueError:
            print("Enter a number.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="?",
        help=("LigParGen .data/.lmp/.zip, MOL2, LAMMPS data, or POSCAR structure"),
    )
    parser.add_argument(
        "--smiles",
        help=("SMILES input for LigParGen; use instead of the positional input file"),
    )
    parser.add_argument("--ff", help="Path to .frc, 'ligpargen', 'gaff2', or 'uff'")
    parser.add_argument("-o", "--output", help="New output directory")
    parser.add_argument(
        "--charge", type=int, help="Integer net charge (GAFF2; default 0)"
    )
    parser.add_argument("--name", help="Output basename; defaults to input basename")
    parser.add_argument(
        "--buffer", type=float, help="GAFF2 box padding in Angstrom (default 5.0)"
    )
    parser.add_argument("--amber-converter", help="Path to amber_to_lammps.py")
    parser.add_argument(
        "--nta",
        help=(
            "Existing .nta file for manual atom types; input must be a LAMMPS .data file"
        ),
    )
    parser.add_argument(
        "--frc-charge",
        default="auto",
        choices=("auto", "gasteiger", "bond_increments", "fixed", "am1bcc", "input"),
        help="Charge method for automatic .frc typing (default: auto)",
    )
    parser.add_argument(
        "--ligpargen-runner",
        choices=("web", "docker"),
        default="web",
        help="Generate from SMILES using the public web service or local Docker",
    )
    parser.add_argument(
        "--ligpargen-charge-model",
        choices=("CM1A-LBCC", "CM1A"),
        default="CM1A-LBCC",
        help="LigParGen partial-charge model",
    )
    parser.add_argument(
        "--ligpargen-opt",
        type=int,
        choices=range(4),
        default=0,
        metavar="{0,1,2,3}",
        help="BOSS optimization iterations",
    )
    parser.add_argument(
        "--ligpargen-timeout",
        type=float,
        default=180.0,
        help="Generation timeout in seconds (default 180)",
    )
    parser.add_argument(
        "--ligpargen-space-url",
        default="https://lsdodda-ligpargen.hf.space",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--ligpargen-docker-image",
        default="ligpargen:dev",
        help="Local LigParGen Docker image (default ligpargen:dev)",
    )
    return parser


def resolve_inputs(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> tuple[Path | None, str]:
    """Prompt for missing options and check that they describe a valid workflow."""
    if args.input and args.smiles:
        parser.error("Use either an input file or --smiles, not both.")
    interactive = (
        (args.input is None and args.smiles is None)
        or args.ff is None
        or args.output is None
    )
    if interactive:
        print("\nLAMMPS FORCE-FIELD PARAMETERIZATION WIZARD")
        print("------------------------------------------")
    if args.ff is None:
        args.ff = ask_forcefield()
    if args.input is None and args.smiles is None:
        if args.ff.lower() == "ligpargen":
            entered = input(
                "SMILES [press Enter to import a downloaded LigParGen file]: "
            ).strip()
            if entered:
                args.smiles = entered
            else:
                args.input = ask_existing_file("LigParGen .data, .lmp, or .zip file")
        else:
            args.input = ask_existing_file("Input MOL2, LAMMPS data, or POSCAR file")

    source: Path | None = None
    if args.input is not None:
        source = Path(args.input).expanduser().resolve()
        if not source.is_file():
            parser.error(f"Input file not found: {source}")
    if args.smiles and args.ff.lower() != "ligpargen":
        parser.error("--smiles is currently supported only with --ff ligpargen.")
    name = args.name or (source.stem if source is not None else "molecule")
    if args.output is None:
        parent = source.parent if source is not None else Path.cwd()
        default_output = parent / f"{name}_{Path(args.ff).stem.upper()}"
        entered = clean_path(input(f"Output directory [{default_output}]: "))
        args.output = entered or str(default_output)
    if args.ff.lower() == "gaff2":
        if args.charge is None:
            args.charge = ask_int("Net molecular charge", 0) if interactive else 0
        if args.buffer is None:
            args.buffer = (
                ask_float("Vacuum/box padding (Angstrom)", 5.0) if interactive else 5.0
            )
        if args.amber_converter is None and interactive:
            entered = clean_path(
                input("Path to amber_to_lammps.py [press Enter for automatic search]: ")
            )
            args.amber_converter = entered or None
    else:
        args.charge = 0 if args.charge is None else args.charge
        args.buffer = 5.0 if args.buffer is None else args.buffer
    return source, name


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    source, name = resolve_inputs(args, parser)
    root = Path(__file__).resolve().parent
    out_dir = Path(args.output).expanduser().resolve()
    try:
        if args.ff.lower() == "ligpargen":
            if args.smiles:
                audit = parameterize_smiles(
                    args.smiles,
                    out_dir,
                    name,
                    args.charge,
                    args.ligpargen_opt,
                    args.ligpargen_charge_model,
                    args.ligpargen_runner,
                    args.ligpargen_timeout,
                    args.ligpargen_space_url,
                    args.ligpargen_docker_image,
                )
            else:
                if source is None:
                    parser.error("LigParGen requires an input file or --smiles.")
                audit = import_ligpargen(source, out_dir, name)
        elif args.ff.lower() == "gaff2":
            if source is None:
                parser.error("GAFF2 requires an input structure file.")
            audit = gaff2(
                source, out_dir, name, args.charge, args.amber_converter, args.buffer
            )
        elif args.ff.lower() == "uff":
            if source is None:
                parser.error("UFF requires an input structure file.")
            audit = uff(source, out_dir, name)
        else:
            audit = None
    except (OSError, PipelineError, LigParGenError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    if audit is not None:
        print(json.dumps(audit, indent=2))
        print(f"\nComplete {audit['forcefield']} artifacts: {out_dir}")
        return 0

    frc = Path(args.ff).expanduser().resolve()
    if not frc.is_file() or frc.suffix.lower() != ".frc":
        parser.error(
            "--ff must be an existing .frc file, 'ligpargen', 'uff', or 'gaff2'"
        )
    if args.nta:
        if source is None:
            parser.error("--nta requires an input LAMMPS data file.")
        nta = Path(args.nta).expanduser().resolve()
        if not nta.is_file() or nta.suffix.lower() != ".nta":
            print(f"[ERROR] Existing .nta file not found: {nta}", file=sys.stderr)
            return 2
        if source.suffix.lower() != ".data":
            print("[ERROR] --nta requires a LAMMPS .data input file.", file=sys.stderr)
            return 2
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f"{name}_{frc.stem}.data"
        code = subprocess.call(
            [
                sys.executable,
                str(root / "build_input_lmp.py"),
                args.input,
                str(nta),
                str(frc),
                str(out),
            ]
        )
        if code == 0:
            print(f"\nUsed existing NTA: {nta}")
            print(f"Generated LAMMPS data: {out}")
        return code
    try:
        if source is None:
            parser.error(
                "An input structure file is required for .frc parameterization."
            )
        audit = run_frc_pipeline(source, frc, out_dir, name, args.frc_charge)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(audit, indent=2))
    print(f"\nGenerated NTA: {audit['nta_file']}")
    print(f"Generated LAMMPS data: {audit['lammps_data']}")
    if audit["status"] != "complete":
        print(
            f"WARNING: coefficient coverage is {audit['total_coverage_percent']:.2f}%"
        )
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

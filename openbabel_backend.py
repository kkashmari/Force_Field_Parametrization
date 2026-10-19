#!/usr/bin/env python3
# Author: Khatereh Kashmari
"""Prepare structures and audit UFF/GAFF setup with Open Babel.

Open Babel evaluates the force field but does not export its bonded parameter
tables in LAMMPS form. This module writes the prepared structure and audit.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


SUPPORTED = {"uff": "UFF", "gaff": "GAFF"}


def executable(name: str) -> str | None:
    return shutil.which(name)


def run_checked(command: list[str], *, stdout_path: Path | None = None) -> str:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n{detail}"
        )
    if stdout_path is not None:
        stdout_path.write_text(result.stdout, encoding="utf-8")
    return result.stdout


def available_forcefields(obabel: str) -> set[str]:
    output = run_checked([obabel, "-L", "forcefields"])
    return {line.split()[0].lower() for line in output.splitlines() if line.strip()}


def charge_policy(ff: str, requested: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    method = requested.lower()
    if ff == "gaff":
        if method == "auto":
            method = "am1-bcc" if executable("antechamber") else "gasteiger"
        if method == "am1-bcc" and not executable("antechamber"):
            raise RuntimeError(
                "GAFF AM1-BCC was requested, but 'antechamber' is not installed."
            )
        if method != "am1-bcc":
            warnings.append(
                "GAFF is normally paired with RESP or AM1-BCC charges; the selected "
                f"'{method}' charges reduce fidelity."
            )
    elif method == "auto":
        method = "gasteiger"
    return method, warnings


def prepare(input_path: str, ff: str, output_prefix: str, charge: str = "auto") -> dict:
    ff_key = ff.lower()
    if ff_key not in SUPPORTED:
        raise ValueError(f"Unsupported Open Babel force field: {ff}")
    obabel = executable("obabel")
    obenergy = executable("obenergy")
    if not obabel or not obenergy:
        raise RuntimeError(
            "Open Babel 3 is required (commands 'obabel' and 'obenergy'). "
            "Install it, for example with: conda install -c conda-forge openbabel"
        )

    installed = available_forcefields(obabel)
    if ff_key not in installed:
        raise RuntimeError(
            f"This Open Babel installation does not advertise {SUPPORTED[ff_key]}. "
            f"Available force fields: {', '.join(sorted(installed)) or 'none'}"
        )

    source = Path(input_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    prefix = Path(output_prefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    typed_mol2 = prefix.with_name(prefix.name + f"_{ff_key}.mol2")
    energy_log = prefix.with_name(prefix.name + f"_{ff_key}_energy.txt")
    audit_path = prefix.with_name(prefix.name + f"_{ff_key}_audit.json")

    charge_method, warnings = charge_policy(ff_key, charge)
    # Add hydrogens before typing; assign Gasteiger charges when requested.
    convert = [obabel, str(source), "-O", str(typed_mol2), "-h"]
    if charge_method == "gasteiger":
        convert.extend(["--partialcharge", "gasteiger"])
    run_checked(convert)

    energy_text = run_checked([obenergy, "-ff", SUPPORTED[ff_key], str(typed_mol2)])
    energy_log.write_text(energy_text, encoding="utf-8")
    setup_ok = (
        "could not setup" not in energy_text.lower()
        and "error" not in energy_text.lower()
    )

    audit = {
        "backend": "openbabel",
        "forcefield": SUPPORTED[ff_key],
        "input": str(source),
        "prepared_structure": str(typed_mol2),
        "energy_log": str(energy_log),
        "charge_method": charge_method,
        "setup_ok": setup_ok,
        "warnings": warnings,
        "lammps_export_complete": False,
        "lammps_export_reason": (
            "Open Babel's public CLI does not expose every assigned bonded/nonbonded "
            "coefficient needed to build an auditable LAMMPS data file. No zero or "
            "guessed coefficients were written."
        ),
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare and audit UFF/GAFF with Open Babel"
    )
    parser.add_argument("input")
    parser.add_argument("--ff", required=True, choices=sorted(SUPPORTED))
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument(
        "--charge", default="auto", choices=("auto", "gasteiger", "am1-bcc")
    )
    args = parser.parse_args()
    try:
        audit = prepare(args.input, args.ff, args.output_prefix, args.charge)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(audit, indent=2))
    return 0 if audit["setup_ok"] else 3


if __name__ == "__main__":
    raise SystemExit(main())

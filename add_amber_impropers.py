#!/usr/bin/env python3
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0
# Required packages: ParmEd (Python); AmberTools topology files
"""Restore AMBER/GAFF improper torsions omitted by AMBER2LAMMPS."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import parmed as pmd


def phase_sign(phase: float) -> int:
    normalized = phase % 360.0
    if math.isclose(normalized, 0.0, abs_tol=1.0e-3):
        return 1
    if math.isclose(normalized, 180.0, abs_tol=1.0e-3):
        return -1
    raise ValueError(
        f"AMBER improper phase {phase} cannot be represented by LAMMPS improper_style cvff"
    )


def add_impropers(data_path: Path, params_path: Path, prmtop_path: Path) -> int:
    topology = pmd.load_file(str(prmtop_path))
    impropers = [item for item in topology.dihedrals if item.improper]
    if not impropers:
        return 0

    signatures: dict[tuple[float, int, int], int] = {}
    rows = []
    for index, item in enumerate(impropers, 1):
        dtype = item.type
        if hasattr(dtype, "__iter__") and not hasattr(dtype, "phi_k"):
            terms = list(dtype)
            if len(terms) != 1:
                raise ValueError("Multi-term AMBER impropers are not supported by cvff export")
            dtype = terms[0]
        signature = (round(float(dtype.phi_k), 10), phase_sign(float(dtype.phase)),
                     int(dtype.per))
        type_id = signatures.setdefault(signature, len(signatures) + 1)
        rows.append((index, type_id, item.atom1.idx + 1, item.atom2.idx + 1,
                     item.atom3.idx + 1, item.atom4.idx + 1))

    data_text = data_path.read_text(encoding="utf-8")
    data_text, n_count = re.subn(
        r"(?m)^(\d+) dihedrals\s*$",
        lambda match: f"{match.group(0)}\n{len(rows)} impropers",
        data_text, count=1)
    data_text, n_types = re.subn(
        r"(?m)^(\d+) dihedral types\s*$",
        lambda match: f"{match.group(0)}\n{len(signatures)} improper types",
        data_text, count=1)
    if n_count != 1 or n_types != 1:
        raise ValueError("Could not locate LAMMPS dihedral header lines for improper insertion")
    data_text = data_text.rstrip() + "\n\nImpropers\n\n"
    data_text += "\n".join("%d %d %d %d %d %d" % row for row in rows) + "\n"
    data_path.write_text(data_text, encoding="utf-8")

    params_text = params_path.read_text(encoding="utf-8").rstrip()
    params_text += "\n\n# Improper coefficients restored from AMBER topology\n"
    for signature, type_id in sorted(signatures.items(), key=lambda item: item[1]):
        force, sign, periodicity = signature
        params_text += f"improper_coeff {type_id} {force:.10g} {sign:d} {periodicity:d}\n"
    params_path.write_text(params_text, encoding="utf-8")
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("params", type=Path)
    parser.add_argument("prmtop", type=Path)
    args = parser.parse_args()
    count = add_impropers(args.data, args.params, args.prmtop)
    print(f"Restored {count} AMBER improper torsion(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

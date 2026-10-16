#!/usr/bin/env python3
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0
# Required packages: Python standard library; local atom_typing.py and build_input_lmp.py
"""Automatic atom typing -> NTA -> LAMMPS data pipeline for .frc files."""

from __future__ import annotations

import json
from pathlib import Path

import atom_typing as typing
import build_input_lmp as builder


CHARGE_LABELS = {
    "auto": "Automatic recommended method", "gasteiger": "Gasteiger",
    "bond_increments": "Bond increments", "fixed": "Fixed charges",
    "am1bcc": "AM1-BCC", "input": "Input charges",
}


def choose_charge_method(frc_path: Path, structure_format: str, requested: str) -> str:
    keys = [item[1] for item in typing.build_charge_menu(str(frc_path), structure_format)]
    if requested == "auto":
        return keys[0]
    if requested not in keys:
        raise ValueError(f"Charge method '{requested}' is unavailable. Available: {', '.join(keys)}")
    return requested


def run_frc_pipeline(source: Path, frc_path: Path, out_dir: Path,
                     name: str, charge_method: str = "auto") -> dict:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ValueError(f"Output directory is not empty: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ff_name = frc_path.stem
    nta_path = out_dir / f"{name}.nta"
    intermediate = out_dir / f"{name}_typed.data"
    final_data = out_dir / f"{name}_{ff_name}.data"
    score_path = out_dir / f"{name}_typing_scores.txt"

    structure_format = typing.detect_input_format(str(source))
    ff_types = typing.parse_frc(str(frc_path))
    typing.precompute_ff_features(ff_types)
    atoms, bonds, metadata, structure_format = typing.parse_structure(str(source))
    neighbours, bond_type_map = typing.build_graph(bonds)
    rings = typing.atom_ring_sizes(neighbours)
    atoms_by_id = {atom["id"]: atom for atom in atoms}
    assignments = typing.assign_types(
        atoms, atoms_by_id, neighbours, bond_type_map, rings, ff_types)
    typing.apply_pcff_iff_hpan_override(
        str(frc_path), atoms, atoms_by_id, neighbours, assignments)
    changed = typing.optimize_types_for_parameter_coverage(
        atoms, bonds, atoms_by_id, neighbours, bond_type_map, rings, ff_types,
        assignments, str(frc_path), report_path=str(score_path))

    selected_charge = choose_charge_method(frc_path, structure_format, charge_method)
    charges = typing.compute_charges(
        selected_charge, atoms, bonds, assignments, ff_types, str(frc_path), str(source))
    charge_label = CHARGE_LABELS.get(selected_charge, selected_charge)
    typing.write_nta(atoms, assignments, ff_name, charge_label, str(nta_path))
    typing.write_data(atoms, bonds, assignments, charges, ff_name, charge_label,
                      atoms_by_id, neighbours, str(intermediate), metadata)

    nta = builder.parse_nta(str(nta_path))
    _header, box, lmp_atoms, lmp_bonds = builder.parse_lammps_data(str(intermediate))
    builder.apply_pcff_iff_hpan_override(str(frc_path), lmp_bonds, nta)
    masses, coeffs, equiv, auto_equiv, ff_class = builder.parse_frc(str(frc_path))
    missing, fallback, stats = builder.write_data(
        str(final_data), str(intermediate), str(frc_path), lmp_atoms, lmp_bonds,
        nta, masses, coeffs, equiv, auto_equiv, box, ff_class)

    totals = {"total": 0, "matched": 0, "missing": 0, "exact": 0, "fallback": 0}
    for item in stats.values():
        for key in totals:
            totals[key] += item.get(key, 0)
    coverage = 100.0 * totals["matched"] / totals["total"] if totals["total"] else 100.0
    audit = {
        "forcefield": ff_name,
        "status": "complete" if totals["missing"] == 0 else "incomplete",
        "input": str(source), "charge_method": selected_charge,
        "typing_assignments_changed_for_coverage": changed,
        "nta_file": str(nta_path), "typing_score_report": str(score_path),
        "typed_intermediate_data": str(intermediate), "lammps_data": str(final_data),
        "coefficient_statistics": stats, "total_coverage_percent": coverage,
        "missing_coefficients": {
            section: ["-".join(key) for key in keys] for section, keys in missing.items()},
        "fallback_summary": {
            method: {section: len(entries) for section, entries in sections.items()}
            for method, sections in fallback.items()},
    }
    (out_dir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit

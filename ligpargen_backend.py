#!/usr/bin/env python3
"""Import and validate native LigParGen LAMMPS OPLS-AA output.

This backend deliberately performs no atom typing or parameter substitution.
LigParGen/BOSS remains the authority for atom types, CM1A(-LBCC) charges, and
all force-field coefficients. The input data file is copied byte-for-byte.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


class LigParGenError(RuntimeError):
    """Raised when a LigParGen artifact is incomplete or inconsistent."""


SECTION_NAMES = (
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Bonds", "Angles",
    "Dihedrals", "Impropers",
)
COUNT_KEYS = ("atoms", "bonds", "angles", "dihedrals", "impropers")
TYPE_KEYS = ("atom", "bond", "angle", "dihedral", "improper")


def _is_ligpargen_data(text: str) -> bool:
    return ("ligpargen" in text[:1000].lower()
            and all(re.search(rf"(?m)^\s*{re.escape(name)}\s*(?:#.*)?$", text)
                    for name in ("Masses", "Pair Coeffs", "Atoms", "Bonds")))


def read_artifact(source: Path) -> tuple[bytes, str | None]:
    """Return a LigParGen data-file payload and optional ZIP member name."""
    if source.suffix.lower() != ".zip":
        payload = source.read_bytes()
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise LigParGenError(f"LigParGen data is not UTF-8 text: {source}") from exc
        if not _is_ligpargen_data(text):
            raise LigParGenError(
                "Input is not recognizable LigParGen LAMMPS data. Download the "
                "LAMMPS data file (or the complete ZIP) from LigParGen."
            )
        return payload, None

    try:
        with zipfile.ZipFile(source) as archive:
            candidates: list[tuple[str, bytes]] = []
            for info in archive.infolist():
                if info.is_dir() or Path(info.filename).suffix.lower() not in {".data", ".lmp"}:
                    continue
                payload = archive.read(info)
                try:
                    text = payload.decode("utf-8")
                except UnicodeDecodeError:
                    continue
                if _is_ligpargen_data(text):
                    candidates.append((info.filename, payload))
    except zipfile.BadZipFile as exc:
        raise LigParGenError(f"Invalid ZIP archive: {source}") from exc
    if not candidates:
        raise LigParGenError("ZIP does not contain recognizable LigParGen LAMMPS data.")
    if len(candidates) > 1:
        names = ", ".join(name for name, _payload in candidates)
        raise LigParGenError(f"ZIP contains multiple LigParGen LAMMPS data files: {names}")
    member, payload = candidates[0]
    return payload, member


def _sections(text: str) -> dict[str, list[list[str]]]:
    found: dict[str, list[list[str]]] = {name: [] for name in SECTION_NAMES}
    current: str | None = None
    headings = {name.lower(): name for name in SECTION_NAMES}
    for raw in text.splitlines():
        body = raw.split("#", 1)[0].strip()
        heading = headings.get(body.lower())
        if heading is not None:
            current = heading
            continue
        if not body:
            continue
        if re.match(r"^[A-Za-z]", body):
            current = None
            continue
        if current is not None:
            found[current].append(body.split())
    return found


def _header_counts(text: str) -> tuple[dict[str, int], dict[str, int]]:
    counts: dict[str, int] = {}
    types: dict[str, int] = {}
    for key in COUNT_KEYS:
        match = re.search(rf"(?m)^\s*(\d+)\s+{key}\s*$", text)
        if match:
            counts[key] = int(match.group(1))
    for key in TYPE_KEYS:
        match = re.search(rf"(?m)^\s*(\d+)\s+{key}\s+types\s*$", text)
        if match:
            types[key] = int(match.group(1))
    missing = [key for key in COUNT_KEYS if key not in counts]
    missing += [f"{key} types" for key in TYPE_KEYS if key not in types]
    if missing:
        raise LigParGenError("Missing LAMMPS header count(s): " + ", ".join(missing))
    return counts, types


def _integer(token: str, context: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise LigParGenError(f"Expected integer in {context}, found {token!r}") from exc


def _validate_ids(rows: list[list[str]], expected: int, context: str) -> None:
    ids = [_integer(row[0], context) for row in rows if row]
    if ids != list(range(1, expected + 1)):
        raise LigParGenError(
            f"{context} IDs are not consecutive 1..{expected}; found {len(ids)} rows"
        )


def validate_ligpargen_text(text: str) -> dict:
    if not _is_ligpargen_data(text):
        raise LigParGenError("The input is not recognizable LigParGen LAMMPS data.")
    counts, types = _header_counts(text)
    sections = _sections(text)
    coefficient_specs = {
        "Masses": (types["atom"], 2),
        "Pair Coeffs": (types["atom"], 3),
        "Bond Coeffs": (types["bond"], 3),
        "Angle Coeffs": (types["angle"], 3),
        "Dihedral Coeffs": (types["dihedral"], 5),
        "Improper Coeffs": (types["improper"], 4),
    }
    for section_name, (expected, columns) in coefficient_specs.items():
        rows = sections[section_name]
        if len(rows) != expected:
            raise LigParGenError(
                f"{section_name}: expected {expected} rows, found {len(rows)}")
        _validate_ids(rows, expected, section_name)
        for row in rows:
            if len(row) < columns:
                raise LigParGenError(f"{section_name} row {row[0]} has too few columns")
            try:
                values = [float(value) for value in row[1:columns]]
            except ValueError as exc:
                raise LigParGenError(
                    f"Non-numeric coefficient in {section_name} row {row[0]}") from exc
            if not all(math.isfinite(value) for value in values):
                raise LigParGenError(f"Non-finite coefficient in {section_name} row {row[0]}")

    topology_specs = {
        "Atoms": (counts["atoms"], 7, "atom", 0),
        "Bonds": (counts["bonds"], 4, "bond", 2),
        "Angles": (counts["angles"], 5, "angle", 3),
        "Dihedrals": (counts["dihedrals"], 6, "dihedral", 4),
        "Impropers": (counts["impropers"], 6, "improper", 4),
    }
    atom_ids = set(range(1, counts["atoms"] + 1))
    charges: list[float] = []
    for section_name, (expected, columns, type_key, atom_columns) in topology_specs.items():
        rows = sections[section_name]
        if len(rows) != expected:
            raise LigParGenError(
                f"{section_name}: expected {expected} rows, found {len(rows)}")
        _validate_ids(rows, expected, section_name)
        for row in rows:
            if len(row) < columns:
                raise LigParGenError(f"{section_name} row {row[0]} has too few columns")
            type_column = 2 if section_name == "Atoms" else 1
            type_id = _integer(row[type_column], f"{section_name} type")
            if not 1 <= type_id <= types[type_key]:
                raise LigParGenError(
                    f"{section_name} row {row[0]} uses type {type_id}, "
                    f"outside 1..{types[type_key]}"
                )
            if section_name == "Atoms":
                try:
                    values = [float(value) for value in row[3:7]]
                except ValueError as exc:
                    raise LigParGenError(
                        f"Invalid charge/coordinate in Atoms row {row[0]}") from exc
                if not all(math.isfinite(value) for value in values):
                    raise LigParGenError(f"Non-finite charge/coordinate in Atoms row {row[0]}")
                charges.append(values[0])
            else:
                referenced = [_integer(value, f"{section_name} atom ID")
                              for value in row[2:2 + atom_columns]]
                invalid = [value for value in referenced if value not in atom_ids]
                if invalid:
                    raise LigParGenError(
                        f"{section_name} row {row[0]} references nonexistent atom(s): {invalid}"
                    )

    net_charge = sum(charges)
    nearest_integer = round(net_charge)
    warnings = []
    if abs(net_charge - nearest_integer) > 1.0e-4:
        warnings.append(
            f"Net charge {net_charge:+.8f} e is not within 1e-4 of an integer; "
            "the LigParGen charges were preserved without modification."
        )
    return {
        "counts": counts,
        "type_counts": types,
        "net_charge": net_charge,
        "nearest_integer_charge": nearest_integer,
        "warnings": warnings,
    }


def repair_known_ligpargen_header_issues(text: str) -> tuple[str, list[str]]:
    """Repair metadata-only defects seen in legacy LigParGen LAMMPS output.

    LigParGen 2.1 can report the number of improper *types* or candidate
    centers as the interaction count while writing a smaller, valid Impropers
    section.  Only that header count is corrected; topology and parameters are
    never created, removed, or altered.
    """
    counts, _types = _header_counts(text)
    actual_impropers = len(_sections(text)["Impropers"])
    if counts["impropers"] == actual_impropers:
        return text, []
    pattern = re.compile(r"(?m)^(\s*)\d+(\s+impropers\s*)$")
    repaired, changes = pattern.subn(
        lambda match: f"{match.group(1)}{actual_impropers}{match.group(2)}",
        text, count=1)
    if changes != 1:
        raise LigParGenError("Could not repair inconsistent LigParGen improper count header")
    note = (f"Corrected LigParGen improper interaction header from "
            f"{counts['impropers']} to {actual_impropers}; coefficient and topology "
            "rows were not changed.")
    return repaired, [note]


def lammps_input(data_name: str) -> str:
    return f"""# Native LigParGen OPLS-AA validation input
# The data file is preserved exactly as generated by LigParGen/BOSS.

units           real
dimension       3
boundary        p p p
atom_style      full

pair_style      lj/cut/coul/long 12.0
pair_modify     mix geometric
bond_style      harmonic
angle_style     harmonic
dihedral_style  opls
improper_style  cvff
special_bonds   lj/coul 0.0 0.0 0.5
kspace_style    pppm 1.0e-5

read_data       {data_name}

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
thermo          1
thermo_style    custom step atoms temp press pe ebond eangle edihed eimp evdwl ecoul
run             0
"""


def import_ligpargen(source: Path, out_dir: Path, name: str) -> dict:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise LigParGenError(f"Output directory is not empty: {out_dir}")
    source_payload, archive_member = read_artifact(source)
    repaired_text, repairs = repair_known_ligpargen_header_issues(
        source_payload.decode("utf-8"))
    payload = repaired_text.encode("utf-8")
    validation = validate_ligpargen_text(repaired_text)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_path = out_dir / f"{name}_ligpargen_opls.data"
    data_path.write_bytes(payload)
    input_path = out_dir / f"in.{name}_ligpargen_opls"
    input_path.write_text(lammps_input(data_path.name), encoding="utf-8")
    audit = {
        "forcefield": "OPLS-AA",
        "parameter_source": "LigParGen/BOSS native LAMMPS output",
        "status": "complete",
        "input": str(source),
        "archive_member": archive_member,
        "source_sha256": hashlib.sha256(source_payload).hexdigest(),
        "output_sha256": hashlib.sha256(payload).hexdigest(),
        "data_file": str(data_path),
        "lammps_input": str(input_path),
        "data_preserved_verbatim": payload == source_payload,
        "metadata_repairs": repairs,
        "tinker_parameters_used": False,
        "atom_typing_performed_locally": False,
        "charges_modified": False,
        "required_lammps_styles": {
            "pair": "lj/cut/coul/long",
            "mixing": "geometric",
            "bond": "harmonic",
            "angle": "harmonic",
            "dihedral": "opls",
            "improper": "cvff",
            "special_bonds": "lj/coul 0.0 0.0 0.5",
        },
        **validation,
    }
    (out_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit

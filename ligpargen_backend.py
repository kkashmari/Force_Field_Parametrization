#!/usr/bin/env python3
"""Validate and import LigParGen LAMMPS output without changing parameters.

Atom types, charges, and coefficient rows come from LigParGen/BOSS. A known
legacy improper-count error is corrected in the header and recorded in the
audit. SMILES jobs can use the web service or a local Docker image.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any


class LigParGenError(RuntimeError):
    """Raised when a LigParGen artifact is incomplete or inconsistent."""


SECTION_NAMES = (
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
)
COUNT_KEYS = ("atoms", "bonds", "angles", "dihedrals", "impropers")
TYPE_KEYS = ("atom", "bond", "angle", "dihedral", "improper")
DEFAULT_SPACE_URL = "https://lsdodda-ligpargen.hf.space"
CHARGE_MODELS = {
    "CM1A-LBCC": "1.14*CM1A-LBCC (neutral molecules)",
    "CM1A": "1.14*CM1A (neutral or charged)",
}


def _is_ligpargen_data(text: str) -> bool:
    header = text[:1000].lower()
    # Older files name LigParGen; newer files name the writer, Leela S. Dodda.
    # This identifies the format; validation below checks the actual contents.
    signature = "ligpargen" in header or "leela s. dodda" in header
    return signature and all(
        re.search(rf"(?m)^\s*{re.escape(name)}\s*(?:#.*)?$", text)
        for name in ("Masses", "Pair Coeffs", "Atoms", "Bonds")
    )


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
                if info.is_dir() or Path(info.filename).suffix.lower() not in {
                    ".data",
                    ".lmp",
                }:
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
        raise LigParGenError(
            f"ZIP contains multiple LigParGen LAMMPS data files: {names}"
        )
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
                f"{section_name}: expected {expected} rows, found {len(rows)}"
            )
        _validate_ids(rows, expected, section_name)
        for row in rows:
            if len(row) < columns:
                raise LigParGenError(f"{section_name} row {row[0]} has too few columns")
            try:
                values = [float(value) for value in row[1:columns]]
            except ValueError as exc:
                raise LigParGenError(
                    f"Non-numeric coefficient in {section_name} row {row[0]}"
                ) from exc
            if not all(math.isfinite(value) for value in values):
                raise LigParGenError(
                    f"Non-finite coefficient in {section_name} row {row[0]}"
                )

    topology_specs = {
        "Atoms": (counts["atoms"], 7, "atom", 0),
        "Bonds": (counts["bonds"], 4, "bond", 2),
        "Angles": (counts["angles"], 5, "angle", 3),
        "Dihedrals": (counts["dihedrals"], 6, "dihedral", 4),
        "Impropers": (counts["impropers"], 6, "improper", 4),
    }
    atom_ids = set(range(1, counts["atoms"] + 1))
    charges: list[float] = []
    for section_name, (
        expected,
        columns,
        type_key,
        atom_columns,
    ) in topology_specs.items():
        rows = sections[section_name]
        if len(rows) != expected:
            raise LigParGenError(
                f"{section_name}: expected {expected} rows, found {len(rows)}"
            )
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
                        f"Invalid charge/coordinate in Atoms row {row[0]}"
                    ) from exc
                if not all(math.isfinite(value) for value in values):
                    raise LigParGenError(
                        f"Non-finite charge/coordinate in Atoms row {row[0]}"
                    )
                charges.append(values[0])
            else:
                referenced = [
                    _integer(value, f"{section_name} atom ID")
                    for value in row[2 : 2 + atom_columns]
                ]
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
        text,
        count=1,
    )
    if changes != 1:
        raise LigParGenError(
            "Could not repair inconsistent LigParGen improper count header"
        )
    note = (
        f"Corrected LigParGen improper interaction header from "
        f"{counts['impropers']} to {actual_impropers}; coefficient and topology "
        "rows were not changed."
    )
    return repaired, [note]


def lammps_input(data_name: str) -> str:
    return f"""# Native LigParGen OPLS-AA validation input
# Force-field parameters and topology come directly from LigParGen/BOSS.

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


def import_ligpargen(
    source: Path, out_dir: Path, name: str, audit_context: dict | None = None
) -> dict:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise LigParGenError(f"Output directory is not empty: {out_dir}")
    source_payload, archive_member = read_artifact(source)
    repaired_text, repairs = repair_known_ligpargen_header_issues(
        source_payload.decode("utf-8")
    )
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
    if audit_context:
        audit.update(audit_context)
    (out_dir / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    return audit


def _json_request(url: str, payload: dict, timeout: float) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "ff-parameterizer/1",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, json.JSONDecodeError) as exc:
        raise LigParGenError(
            f"Could not submit the SMILES to LigParGen: {exc}"
        ) from exc


def _read_gradio_result(url: str, timeout: float) -> list[Any]:
    deadline = time.monotonic() + timeout
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            event = ""
            for raw in response:
                if time.monotonic() > deadline:
                    raise LigParGenError(
                        "LigParGen timed out while generating parameters."
                    )
                line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
                if line.startswith("event:"):
                    event = line.split(":", 1)[1].strip()
                elif line.startswith("data:"):
                    value = line.split(":", 1)[1].strip()
                    if event == "error":
                        raise LigParGenError(
                            f"LigParGen rejected the molecule: {value}"
                        )
                    if event == "complete":
                        try:
                            result = json.loads(value)
                        except json.JSONDecodeError as exc:
                            raise LigParGenError(
                                "LigParGen returned an unreadable result."
                            ) from exc
                        if not isinstance(result, list):
                            raise LigParGenError(
                                "LigParGen returned an unexpected result."
                            )
                        return result
    except urllib.error.URLError as exc:
        raise LigParGenError(f"Lost connection to LigParGen: {exc}") from exc
    raise LigParGenError("LigParGen ended without returning an output ZIP.")


def _find_zip_url(value: Any) -> str | None:
    if isinstance(value, dict):
        url = value.get("url")
        if isinstance(url, str) and ".zip" in url.lower():
            return url
        for child in value.values():
            found = _find_zip_url(child)
            if found:
                return found
    elif isinstance(value, list):
        for child in value:
            found = _find_zip_url(child)
            if found:
                return found
    return None


def _download(url: str, destination: Path, timeout: float) -> None:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https":
        raise LigParGenError(
            "LigParGen returned a non-HTTPS download URL; refusing it."
        )
    try:
        request = urllib.request.Request(
            url, headers={"User-Agent": "ff-parameterizer/1"}
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            destination.write_bytes(response.read())
    except urllib.error.URLError as exc:
        raise LigParGenError(f"Could not download the LigParGen result: {exc}") from exc


def generate_with_ligpargen_space(
    smiles: str,
    destination: Path,
    charge: int = 0,
    opt: int = 0,
    charge_model: str = "CM1A-LBCC",
    timeout: float = 180.0,
    space_url: str = DEFAULT_SPACE_URL,
) -> dict:
    """Submit SMILES to LigParGen's public Gradio API and download its ZIP."""
    smiles = smiles.strip()
    if not smiles:
        raise LigParGenError("SMILES cannot be empty.")
    if charge not in {-2, -1, 0, 1, 2}:
        raise LigParGenError("LigParGen charge must be -2, -1, 0, 1, or 2.")
    if opt not in {0, 1, 2, 3}:
        raise LigParGenError("LigParGen optimization level must be 0, 1, 2, or 3.")
    if charge_model not in CHARGE_MODELS:
        raise LigParGenError("Charge model must be CM1A or CM1A-LBCC.")
    if charge_model == "CM1A-LBCC" and charge != 0:
        raise LigParGenError(
            "CM1A-LBCC is only available for neutral molecules (charge 0)."
        )

    base = space_url.rstrip("/")
    submit_url = f"{base}/gradio_api/call/run_ligpargen"
    submitted = _json_request(
        submit_url,
        {"data": [smiles, None, str(opt), CHARGE_MODELS[charge_model], str(charge)]},
        timeout,
    )
    event_id = submitted.get("event_id")
    if not isinstance(event_id, str) or not event_id:
        raise LigParGenError("LigParGen did not return a job identifier.")
    result = _read_gradio_result(
        f"{submit_url}/{urllib.parse.quote(event_id)}", timeout
    )
    zip_url = _find_zip_url(result)
    if not zip_url:
        raise LigParGenError(
            "LigParGen completed but did not return a downloadable ZIP."
        )
    _download(urllib.parse.urljoin(base + "/", zip_url), destination, timeout)
    if not zipfile.is_zipfile(destination):
        raise LigParGenError("The LigParGen download is not a valid ZIP archive.")
    return {"runner": "web", "service": base, "event_id": event_id}


def generate_with_ligpargen_docker(
    smiles: str,
    destination: Path,
    charge: int = 0,
    opt: int = 0,
    charge_model: str = "CM1A-LBCC",
    timeout: float = 180.0,
    image: str = "ligpargen:dev",
) -> dict:
    """Run the current LigParGen CLI image, which must include the user's BOSS."""
    if charge_model == "CM1A-LBCC" and charge != 0:
        raise LigParGenError(
            "CM1A-LBCC is only available for neutral molecules (charge 0)."
        )
    with tempfile.TemporaryDirectory(prefix="ligpargen_docker_") as tmp:
        work = Path(tmp)
        command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{work}:/work",
            "-w",
            "/work",
            image,
            "-s",
            smiles,
            "-r",
            "LPG",
            "-c",
            str(charge),
            "-o",
            str(opt),
        ]
        if charge_model == "CM1A-LBCC":
            command.append("-l")
        try:
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
        except FileNotFoundError as exc:
            raise LigParGenError("Docker is not installed or is not on PATH.") from exc
        except subprocess.TimeoutExpired as exc:
            raise LigParGenError("The local LigParGen Docker job timed out.") from exc
        if completed.returncode != 0:
            tail = completed.stdout[-3000:].strip()
            raise LigParGenError(
                f"LigParGen Docker failed (exit {completed.returncode}):\n{tail}"
            )
        candidates = [path for path in work.rglob("*.zip") if zipfile.is_zipfile(path)]
        if len(candidates) != 1:
            raise LigParGenError(
                f"Expected one LigParGen ZIP from Docker, found {len(candidates)}."
            )
        destination.write_bytes(candidates[0].read_bytes())
        return {
            "runner": "docker",
            "docker_image": image,
            "log": completed.stdout[-10000:],
        }


def parameterize_smiles(
    smiles: str,
    out_dir: Path,
    name: str,
    charge: int = 0,
    opt: int = 0,
    charge_model: str = "CM1A-LBCC",
    runner: str = "web",
    timeout: float = 180.0,
    space_url: str = DEFAULT_SPACE_URL,
    docker_image: str = "ligpargen:dev",
) -> dict:
    """Generate native LigParGen output from SMILES, then validate its LAMMPS data."""
    if out_dir.exists() and any(out_dir.iterdir()):
        raise LigParGenError(f"Output directory is not empty: {out_dir}")
    with tempfile.TemporaryDirectory(prefix="ligpargen_smiles_") as tmp:
        archive = Path(tmp) / "ligpargen_result.zip"
        if runner == "web":
            generation = generate_with_ligpargen_space(
                smiles, archive, charge, opt, charge_model, timeout, space_url
            )
        elif runner == "docker":
            generation = generate_with_ligpargen_docker(
                smiles, archive, charge, opt, charge_model, timeout, docker_image
            )
        else:
            raise LigParGenError("LigParGen runner must be 'web' or 'docker'.")
        context = {
            "input": "SMILES",
            "input_smiles": smiles,
            "ligpargen_generation": {
                **generation,
                "formal_charge": charge,
                "optimization_iterations": opt,
                "charge_model": charge_model,
            },
        }
        return import_ligpargen(archive, out_dir, name, context)

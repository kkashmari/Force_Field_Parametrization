#!/usr/bin/env python3
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0
# Required packages: Python standard library
"""
Build a LAMMPS data file from: Required input files
  - a LAMMPS data file containing Atoms and Bonds
  - an .nta atom-type assignment file
  - a .frc force-field file
"""

import itertools
import os
import re
import sys
from collections import defaultdict


CLASS2_COEFF_SPECS = [
    ('Pair Coeffs', 'nonbond(9-6)', 1, 'lj/class2/coul/long'),
    ('Bond Coeffs', ('quartic_bond', 'quadratic_bond'), 2, 'class2'),
    ('Angle Coeffs', ('quartic_angle', 'quadratic_angle'), 3, 'class2'),
    ('Dihedral Coeffs', ('torsion_3', 'torsion_1'), 4, 'class2'),
    ('Improper Coeffs', ('wilson_out_of_plane', 'out_of_plane'), 4, 'class2'),
    ('BondBond Coeffs', 'bond-bond', 3, 'class2'),
    ('BondAngle Coeffs', 'bond-angle', 3, 'class2'),
    ('AngleAngleTorsion Coeffs', 'angle-angle-torsion_1', 4, 'class2'),
    ('EndBondTorsion Coeffs', 'end_bond-torsion_3', 4, 'class2'),
    ('MiddleBondTorsion Coeffs', 'middle_bond-torsion_3', 4, 'class2'),
    ('BondBond13 Coeffs', 'bond-bond_1_3', 4, 'class2'),
    ('AngleTorsion Coeffs', 'angle-torsion_3', 4, 'class2'),
    ('AngleAngle Coeffs', 'angle-angle', 4, 'class2'),
]

CLASS1_COEFF_SPECS = [
    ('Pair Coeffs', ('nonbond(12-6)', 'nonbond(9-6)'), 1, 'lj/cut/coul/long'),
    ('Bond Coeffs', ('quadratic_bond',), 2, 'harmonic'),
    ('Angle Coeffs', ('quadratic_angle', 'quartic_angle'), 3, 'harmonic'),
    ('Dihedral Coeffs', ('torsion_1',), 4, 'harmonic'),
    ('Improper Coeffs', ('out_of_plane', 'wilson_out_of_plane'), 4, 'harmonic'),
]

OPLS_COEFF_SPECS = [
    ('Pair Coeffs', ('nonbond(12-6)',), 1, 'lj/cut/coul/long'),
    ('Bond Coeffs', ('quadratic_bond',), 2, 'harmonic'),
    ('Angle Coeffs', ('quadratic_angle',), 3, 'harmonic'),
    ('Dihedral Coeffs', ('torsion_1',), 4, 'opls'),
    ('Improper Coeffs', ('out_of_plane',), 4, 'harmonic'),
]

ALL_COEFF_SPECS = CLASS2_COEFF_SPECS + [
    (title, frc_section, ntypes, style)
    for title, frc_sections, ntypes, style in CLASS1_COEFF_SPECS
    for frc_section in frc_sections
]

SECTION_BY_FRC = {}
for _out, _frc_sections, _ntypes, _style in ALL_COEFF_SPECS:
    if isinstance(_frc_sections, str):
        _frc_sections = (_frc_sections,)
    for _frc_section in _frc_sections:
        SECTION_BY_FRC[_frc_section] = (_out, _ntypes, _style)
COEFF_VALUE_COUNTS = {
    'Pair Coeffs': 2,
    'Bond Coeffs': 4,
    'Angle Coeffs': 4,
    'Dihedral Coeffs': 6,
    'Improper Coeffs': 2,
    'BondBond Coeffs': 3,
    'BondAngle Coeffs': 4,
    'AngleAngleTorsion Coeffs': 3,
    'EndBondTorsion Coeffs': 8,
    'MiddleBondTorsion Coeffs': 4,
    'BondBond13 Coeffs': 3,
    'AngleTorsion Coeffs': 8,
    'AngleAngle Coeffs': 6,
}
CLASS1_COEFF_VALUE_COUNTS = {
    'Pair Coeffs': 2,
    'Bond Coeffs': 2,
    'Angle Coeffs': 2,
    'Dihedral Coeffs': 3,
    'Improper Coeffs': 2,
}
OPLS_COEFF_VALUE_COUNTS = {
    'Pair Coeffs': 2, 'Bond Coeffs': 2, 'Angle Coeffs': 2,
    'Dihedral Coeffs': 4, 'Improper Coeffs': 2,
}
CLASS2_MARKER_SECTIONS = {
    'bond-bond',
    'bond-angle',
    'angle-angle-torsion',
    'end-bond-torsion',
    'middle-bond-torsion',
    'bond-bond_1_3',
    'angle-torsion',
    'angle-angle',
}
HPAN_HC_FF_FILE = 'pcff_iff_v1_5_CNT_poly_solv_hpan_hc_mod.frc'
GENERAL_BOND_EQUIV = {
    ('C', 'H'): ('c_', 'h_'),
    ('CC', 'CT'): ('c', 'c'),
    ('CC', 'CC'): ('c', 'c'),
    ('CC', 'HC'): ('c', 'h'),
    ('C', 'C'): ('c_', 'c_'),
    ('C', 'O'): ('c', 'o'),
    ('CC', 'OH'): ('c', 'o'),
    ('O', 'H'): ('o_', 'h_'),
    ('C', 'N'): ('c_', 'n_'),
    ('C', 'S'): ('c_', 's_'),
}
NUMERIC_RE = re.compile(r'^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?$')


def prompt_path(label):
    while True:
        raw = input(f'  Path to {label}: ').strip()
        if raw and os.path.isfile(raw):
            return raw
        print(f'  File not found: {raw}')


def prompt_output():
    while True:
        raw = input('  Output .data file name: ').strip()
        if raw:
            return raw if raw.endswith('.data') else raw + '.data'


def clean_num(token):
    return float(token.replace('D', 'E').replace('d', 'e'))


def is_number(token):
    return bool(NUMERIC_RE.match(token))


def split_comment(line):
    if '!' in line:
        line = line.split('!', 1)[0]
    return line.strip()


def section_name(line):
    s = line.strip()
    if s.startswith('#'):
        return s[1:].split()[0].lower()
    return None


def parse_nta(path):
    atom_types = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.split('#', 1)[0].strip()
            if not s or s.lower().startswith(('new ', 'style ')):
                continue
            parts = s.split()
            if len(parts) >= 2 and parts[0].isdigit():
                atom_types[int(parts[0])] = parts[1]
    if not atom_types:
        raise ValueError(f'No atom assignments were found in {path}')
    return atom_types


def parse_lammps_data(path):
    header = []
    atoms = {}
    bonds = []
    box = []
    current = None

    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            stripped = line.strip()
            head = stripped.split('#', 1)[0].strip()

            if not current:
                if re.match(r'^[+-]?\d', head) and any(k in head for k in ('xlo', 'ylo', 'zlo')):
                    box.append(line)
                elif stripped:
                    header.append(line)

            if stripped.startswith('Atoms'):
                current = 'Atoms'
                continue
            if stripped.startswith('Bonds'):
                current = 'Bonds'
                continue
            if stripped.startswith(('Angles', 'Dihedrals', 'Impropers')):
                current = None
                continue
            if stripped.endswith('Coeffs') or ' Coeffs' in stripped:
                current = None
                continue
            if not head:
                continue

            parts = head.split()
            if current == 'Atoms' and len(parts) >= 7 and parts[0].isdigit():
                aid = int(parts[0])
                atoms[aid] = {
                    'id': aid,
                    'mol': int(parts[1]),
                    'old_type': int(parts[2]),
                    'charge': clean_num(parts[3]),
                    'x': clean_num(parts[4]),
                    'y': clean_num(parts[5]),
                    'z': clean_num(parts[6]),
                    'tail': parts[7:],
                }
            elif current == 'Bonds' and len(parts) >= 4 and parts[0].isdigit():
                bonds.append({
                    'id': int(parts[0]),
                    'old_type': int(parts[1]),
                    'a1': int(parts[2]),
                    'a2': int(parts[3]),
                })

    if not atoms or not bonds:
        raise ValueError('Input data file must contain Atoms and Bonds sections')
    return header, box, atoms, bonds


def parse_external_equivalence(path):
    """Read standard/automatic mappings from an equivalence reference file."""
    standard, automatic = {}, {}
    current = None
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            s = raw.strip()
            lower = s.lower()
            if lower.startswith('# standard equivalence table'):
                current = 'standard'; continue
            if lower.startswith('# automatic equivalence table'):
                current = 'automatic'; continue
            if s.startswith('#'):
                if lower.startswith(('# atom metadata', '# canonicalization',
                                     '# general_assumed', '# parameter sections')):
                    current = None
                continue
            if not s or s.startswith('!'):
                continue
            parts = s.split()
            if current == 'standard' and len(parts) >= 6:
                standard.setdefault(parts[0], parts[1:])
            elif current == 'automatic' and len(parts) >= 10:
                automatic.setdefault(parts[0], parts[1:10])
    return standard, automatic


def find_external_equivalence(frc_path):
    stem = os.path.splitext(os.path.basename(frc_path))[0]
    name = f'equivalence_{stem}.txt'
    frc_dir = os.path.dirname(os.path.realpath(frc_path))
    script_dir = os.path.dirname(os.path.realpath(__file__))
    candidates = (
        os.path.join(frc_dir, name),
        os.path.join(frc_dir, 'equivalence_refs', name),
        os.path.join(script_dir, 'frc_files', 'equivalence_refs', name),
    )
    return next((p for p in candidates if os.path.isfile(p)), None)


def parse_frc(path):
    masses = {}
    coeffs = defaultdict(dict)
    equiv = {}
    auto_equiv = {}
    current = None
    seen_sections = set()

    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            s = split_comment(raw)
            if not s:
                continue
            sec = section_name(s)
            if sec is not None:
                current = sec
                seen_sections.add(sec)
                continue
            if s.startswith('>') or current is None:
                continue

            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                clean_num(parts[0])
                int(float(parts[1]))
            except ValueError:
                continue

            body = parts[2:]
            if current == 'atom_types' and len(body) >= 3:
                try:
                    masses[body[0]] = clean_num(body[1])
                except ValueError:
                    pass
                continue

            if current in ('equivalence', 'auto_equivalence') and len(body) >= 2:
                table = auto_equiv if current == 'auto_equivalence' else equiv
                table[body[0]] = body[1:]
                continue

            if current not in SECTION_BY_FRC:
                continue

            _, ntypes, _ = SECTION_BY_FRC[current]
            if len(body) <= ntypes:
                continue
            key = tuple(body[:ntypes])
            vals = []
            for tok in body[ntypes:]:
                if is_number(tok):
                    vals.append(clean_num(tok))
            if vals:
                coeffs[current].setdefault(key, vals)

    # Infer the functional form from the parameter sections, not the filename.
    # CVFF/OPLS files may contain optional cross-term sections but remain class I;
    # quartic bonds/angles, torsion_3, and 9-6 nonbond terms identify class II.
    class2_primary = bool(
        {'quartic_bond', 'quartic_angle', 'torsion_3', 'nonbond(9-6)'}
        & seen_sections)
    opls_primary = (not class2_primary and
                    any(len(vals) >= 4 for vals in coeffs.get('torsion_1', {}).values()))
    ff_class = 'class2' if class2_primary else ('opls' if opls_primary else 'class1')

    external_path = find_external_equivalence(path)
    if external_path:
        external_equiv, external_auto = parse_external_equivalence(external_path)
        for atom_type, row in external_equiv.items():
            equiv.setdefault(atom_type, row)
        for atom_type, row in external_auto.items():
            auto_equiv.setdefault(atom_type, row)
    return masses, coeffs, equiv, auto_equiv, ff_class


def build_graph(bonds):
    neighbours = defaultdict(list)
    seen = set()
    for b in bonds:
        a, c = b['a1'], b['a2']
        if a == c:
            continue
        key = tuple(sorted((a, c)))
        if key in seen:
            continue
        seen.add(key)
        neighbours[a].append(c)
        neighbours[c].append(a)
    for nbs in neighbours.values():
        nbs.sort()
    return neighbours


def apply_pcff_iff_hpan_override(frc_path, bonds, nta):
    if os.path.basename(frc_path) != HPAN_HC_FF_FILE:
        return 0
    neighbours = build_graph(bonds)
    changed = 0
    for aid, atom_type in list(nta.items()):
        if atom_type != 'hc':
            continue
        if any(nta.get(nb, '').lower().startswith('c')
               for nb in neighbours.get(aid, [])):
            nta[aid] = 'hpan'
            changed += 1
    return changed


def canonical_bond(types):
    a, b = types
    return tuple(sorted((a, b)))


def canonical_angle(types):
    a, b, c = types
    return min((a, b, c), (c, b, a))


def canonical_dihedral(types):
    return min(tuple(types), tuple(reversed(types)))


def canonical_improper(types):
    first, center, third, fourth = types
    ends = sorted((first, third, fourth))
    return (ends[0], center, ends[1], ends[2])


def atom_type_tuple(atom_ids, nta):
    return tuple(nta[i] for i in atom_ids)


def build_topology(bonds, neighbours, nta):
    bond_items = [tuple(sorted((b['a1'], b['a2']))) for b in bonds]
    bond_items = sorted(set(bond_items))

    angle_items = []
    for center, nbs in neighbours.items():
        for a, c in itertools.combinations(nbs, 2):
            angle_items.append((a, center, c))
    angle_items.sort()

    dihedral_set = set()
    for b, c in bond_items:
        for a in neighbours[b]:
            if a == c:
                continue
            for d in neighbours[c]:
                if d == b or d == a:
                    continue
                dihedral_set.add(min((a, b, c, d), (d, c, b, a)))
    dihedral_items = sorted(dihedral_set)

    improper_items = []
    for center, nbs in neighbours.items():
        if len(nbs) < 3:
            continue
        for trio in itertools.combinations(nbs, 3):
            ids = (trio[0], center, trio[1], trio[2])
            improper_items.append(ids)
    improper_items.sort()

    return {
        'Bonds': bond_items,
        'Angles': angle_items,
        'Dihedrals': dihedral_items,
        'Impropers': improper_items,
    }


def unique_type_ids(items, nta, canonical):
    ids = {}
    ordered = []
    for item in items:
        key = canonical(atom_type_tuple(item, nta))
        if key not in ids:
            ids[key] = len(ids) + 1
            ordered.append(key)
    return ids, ordered


def _equivalence_indices(section, ntypes, automatic=False,
                         has_bond_increment=False):
    """Return column indices after ``Type`` in an .frc equivalence row.

    Standard CVFF tables contain::

        NonB Bond Angle Torsion OOP

    CFF91-derived tables (including PCFF and COMPASS) contain an additional
    ``Bond Inct`` column between Bond and Angle.  Detecting that layout from
    each parsed row avoids force-field-name special cases.
    """
    bond_sections = {'quadratic_bond', 'quartic_bond'}
    angle_sections = {'quadratic_angle', 'quartic_angle', 'bond-bond', 'bond-angle'}
    torsion_sections = {'torsion_1', 'torsion_3', 'angle-angle-torsion_1',
                        'end_bond-torsion_3', 'middle_bond-torsion_3',
                        'bond-bond_1_3', 'angle-torsion_3'}
    oop_sections = {'out_of_plane', 'wilson_out_of_plane', 'angle-angle'}
    if not automatic:
        shift = 1 if has_bond_increment else 0
        idx = (1 if section in bond_sections
               else 2 + shift if section in angle_sections
               else 3 + shift if section in torsion_sections
               else 4 + shift if section in oop_sections
               else 0)
        return [idx] * ntypes
    if section in bond_sections:
        return [2] * ntypes
    if section in angle_sections:
        return [3, 4, 3][:ntypes]
    if section in torsion_sections:
        return [5, 6, 6, 5][:ntypes]
    if section in oop_sections:
        return [7, 8, 7, 7][:ntypes]
    return [0] * ntypes


def _map_equivalence(key, section, table, automatic=False):
    mapped = []
    for position, atom_type in enumerate(key):
        row = table.get(atom_type, [])
        # Six standard values means NonB/Bond/BondInct/Angle/Torsion/OOP.
        # Automatic tables have their own fixed nine-column schema.
        has_bond_increment = not automatic and len(row) >= 6
        indices = _equivalence_indices(
            section, len(key), automatic, has_bond_increment)
        idx = indices[position]
        mapped.append(row[idx] if idx < len(row) else atom_type)
    return tuple(mapped)


def equivalence_type(t, section, equiv, auto_equiv):
    """Compatibility helper for callers needing a single standard mapping."""
    return _map_equivalence((t,), section, equiv, False)[0]


def equivalent_key(key, section, equiv, auto_equiv):
    if section == 'nonbond(9-6)':
        return key
    return _map_equivalence(key, section, equiv, False)


def equivalence_keys(key, section, equiv, auto_equiv):
    """Return standard and automatic equivalence keys without collapsing them."""
    keys = []
    for table, automatic in ((equiv, False), (auto_equiv, True)):
        mapped = _map_equivalence(key, section, table, automatic)
        if mapped != tuple(key) and mapped not in keys:
            keys.append(mapped)
    return keys


def wildcard_match(pattern, key):
    # .frc files use both '*' and ranked forms such as '*7' as wildcards.
    return len(pattern) == len(key) and all(
        p.startswith('*') or p == k for p, k in zip(pattern, key))


def type_family(atom_type):
    m = re.match(r'[A-Za-z]+', atom_type.replace('_', ''))
    return m.group(0).lower() if m else atom_type.lower()


def general_element(atom_type):
    family = type_family(atom_type)
    if family.startswith('c'):
        return 'C'
    if family.startswith('h'):
        return 'H'
    if family.startswith('o'):
        return 'O'
    if family.startswith('n'):
        return 'N'
    if family.startswith('s'):
        return 'S'
    return family.upper()


def ordered_key_variants(key, canonical):
    variants = []

    def add(k):
        k = tuple(k)
        if k not in variants:
            variants.append(k)

    add(key)
    add(canonical(key))
    add(tuple(reversed(key)))
    add(canonical(tuple(reversed(key))))
    return variants


def similar_score(pattern, key):
    if len(pattern) != len(key) or '*' in pattern:
        return None
    score = 0
    for p, k in zip(pattern, key):
        if p == k:
            score += 4
        elif type_family(p) == type_family(k):
            score += 2
        else:
            return None
    return score


def lookup_similar_coeff(keys, table):
    best = None
    best_score = None
    for key in keys:
        for pattern, vals in table.items():
            score = similar_score(pattern, key)
            if score is None:
                continue
            if best_score is None or score > best_score:
                best = (vals, pattern)
                best_score = score
    return best


def lookup_exact_equiv_wildcard(raw_key, section, table, equiv, auto_equiv, canonical):
    exact_keys = ordered_key_variants(tuple(raw_key), canonical)
    for key in exact_keys:
        if key in table:
            return table[key], 'exact', key, exact_keys

    equiv_keys = []
    for equiv_key in equivalence_keys(raw_key, section, equiv, auto_equiv):
        for key in ordered_key_variants(equiv_key, canonical):
            if key not in equiv_keys:
                equiv_keys.append(key)
            if key in table:
                return table[key], 'equivalence', key, exact_keys + equiv_keys

    wildcard_best = None
    wildcard_specificity = -1
    for key in exact_keys + equiv_keys:
        for pattern, vals in table.items():
            if wildcard_match(pattern, key):
                # Prefer the pattern that fixes the most atom positions. This
                # reproduces .frc ranked wildcards such as '*7 c\'_ o_' taking
                # precedence over the less-specific '* c\'_ *'.
                specificity = sum(not p.startswith('*') for p in pattern)
                if specificity > wildcard_specificity:
                    wildcard_best = (vals, pattern)
                    wildcard_specificity = specificity
    if wildcard_best is not None:
        vals, pattern = wildcard_best
        return vals, 'wildcard', pattern, exact_keys + equiv_keys

    return None, 'zero', None, exact_keys + equiv_keys


def lookup_general_bond_equiv(raw_key, section, table, equiv, auto_equiv, canonical):
    if section != 'quadratic_bond' or len(raw_key) != 2:
        return None, None
    element_key = tuple(general_element(t) for t in raw_key)
    mapped = GENERAL_BOND_EQUIV.get(element_key)
    if mapped is None:
        mapped = GENERAL_BOND_EQUIV.get(tuple(reversed(element_key)))
        if mapped is not None:
            mapped = tuple(reversed(mapped))
    if mapped is None:
        return None, None
    vals, _method, matched_key, _searched = lookup_exact_equiv_wildcard(
        mapped, section, table, equiv, auto_equiv, canonical)
    if vals is not None:
        return vals, matched_key
    return None, mapped


def lookup_coeff(raw_key, section, coeffs, equiv, auto_equiv, canonical):
    table = coeffs.get(section, {})
    raw_key = tuple(raw_key)

    vals, method, matched_key, searched_keys = lookup_exact_equiv_wildcard(
        raw_key, section, table, equiv, auto_equiv, canonical)
    if vals is not None:
        return vals, method, matched_key

    if 'hpan' in raw_key:
        hc_key = tuple('hc' if t == 'hpan' else t for t in raw_key)
        vals, _method, matched_key, _hc_searched = lookup_exact_equiv_wildcard(
            hc_key, section, table, equiv, auto_equiv, canonical)
        if vals is not None:
            return vals, 'hpan_to_hc', matched_key

    vals, matched_key = lookup_general_bond_equiv(
        raw_key, section, table, equiv, auto_equiv, canonical)
    if vals is not None:
        return vals, 'general_equiv', matched_key

    return None, 'zero', None


def transform_coeff_values(title, ff_class, values, source_section=None):
    if ff_class == 'class2':
        count = COEFF_VALUE_COUNTS[title]
        if title == 'Pair Coeffs' and len(values) >= 2:
            # PCFF .frc nonbond(9-6) is r,epsilon; LAMMPS is epsilon,sigma.
            values = [values[1], values[0]]
        elif title == 'Dihedral Coeffs' and source_section == 'torsion_1':
            # K*[1+cos(n*phi-phi0)] -> class2 Vn/phase representation.
            converted = [0.0] * 6
            if len(values) >= 3:
                amplitude, periodicity, phase = values[:3]
                n = int(round(periodicity))
                if 1 <= n <= 3:
                    converted[2*(n-1)] = amplitude
                    # Match the established all2lmp class-II conversion:
                    # 180 degrees is represented as 0; other phases are kept.
                    converted[2*(n-1)+1] = 0.0 if phase == 180.0 else phase
            values = converted
        return (values[:count] + [0.0] * count)[:count]

    if ff_class == 'opls':
        count = OPLS_COEFF_VALUE_COUNTS[title]
        if title == 'Pair Coeffs' and len(values) >= 2:
            # Convert .frc A/r^12-B/r^6 to LAMMPS epsilon,sigma.
            a12, b6 = values[:2]
            if a12 > 0.0 and b6 > 0.0:
                epsilon = b6*b6/(4.0*a12)
                sigma = (a12/b6)**(1.0/6.0)
                values = [epsilon, sigma]
            else:
                values = [0.0, 0.0]
        elif title in ('Bond Coeffs', 'Angle Coeffs') and len(values) >= 2:
            values = [values[1], values[0]]
        return (values[:count] + [0.0] * count)[:count]

    count = CLASS1_COEFF_VALUE_COUNTS[title]
    if title == 'Pair Coeffs' and len(values) >= 2:
        a12, b6 = values[:2]
        if a12 > 0.0 and b6 > 0.0:
            values = [b6*b6/(4.0*a12), (a12/b6)**(1.0/6.0)]
        else:
            values = [0.0, 0.0]
    elif title in ('Bond Coeffs', 'Angle Coeffs') and len(values) >= 2:
        values = [values[1], values[0]]
    else:
        values = values[:count]
    return (values + [0.0] * count)[:count]


def lookup_any_coeff(raw_key, frc_sections, coeffs, equiv, auto_equiv, canonical):
    best_zero = (None, 'zero', None, None)
    for section in frc_sections:
        vals, method, matched_key = lookup_coeff(raw_key, section, coeffs, equiv, auto_equiv, canonical)
        if vals is not None:
            return vals, method, matched_key, section
    return best_zero


def fmt_values(values):
    return ''.join(f'{v:16.8f}' for v in values)


def equilibrium_value(key, section, coeffs, equiv, auto_equiv, canonical):
    """Look up R0 or theta0 from the primary class-II coefficient table."""
    vals, _method, _matched = lookup_coeff(
        key, section, coeffs, equiv, auto_equiv, canonical)
    return vals[0] if vals else 0.0


def complete_class2_cross_term(title, key, values, coeffs, equiv, auto_equiv):
    """Append the equilibrium geometry required by LAMMPS class2 styles."""
    k = tuple(key)
    bond = lambda pair: equilibrium_value(
        pair, 'quartic_bond', coeffs, equiv, auto_equiv, canonical_bond)
    angle = lambda trio: equilibrium_value(
        trio, 'quartic_angle', coeffs, equiv, auto_equiv, canonical_angle)

    if title == 'AngleAngle Coeffs':
        permutations = ((k[0], k[1], k[2], k[3]),
                        (k[2], k[1], k[0], k[3]),
                        (k[0], k[1], k[3], k[2]))
        couplings = []
        for pattern in permutations:
            vals, _method, _matched = lookup_coeff(
                pattern, 'angle-angle', coeffs, equiv, auto_equiv,
                canonical_improper)
            couplings.append(vals[0] if vals else 0.0)
        refs = [angle((k[0], k[1], k[2])), angle((k[0], k[1], k[3])),
                angle((k[2], k[1], k[3]))]
        return couplings + refs
    if title in ('BondBond Coeffs', 'BondAngle Coeffs'):
        refs = [bond((k[0], k[1])), bond((k[1], k[2]))]
    elif title in ('AngleAngleTorsion Coeffs', 'AngleTorsion Coeffs'):
        refs = [angle((k[0], k[1], k[2])), angle((k[1], k[2], k[3]))]
    elif title in ('EndBondTorsion Coeffs', 'BondBond13 Coeffs'):
        refs = [bond((k[0], k[1])), bond((k[2], k[3]))]
    elif title == 'MiddleBondTorsion Coeffs':
        refs = [bond((k[1], k[2]))]
    else:
        refs = []
    count = COEFF_VALUE_COUNTS[title]
    coupling_count = count - len(refs)
    if title == 'BondAngle Coeffs' and len(values) == 1:
        values = values * 2
    if title in ('EndBondTorsion Coeffs', 'AngleTorsion Coeffs') and len(values) == 3:
        values = values * 2
    return (values[:coupling_count] + [0.0] * coupling_count)[:coupling_count] + refs


def comment(types, extra=''):
    base = ''.join(f'{t:>10} ' for t in types)
    return f'# {base}{extra}'.rstrip()


def write_coeff_section(fh, title, style, ordered_keys, frc_sections,
                        coeffs, equiv, auto_equiv, canonical, ff_class):
    fh.write(f'{title}  # {style}\n\n')
    missing = []
    fallback = defaultdict(list)
    stats = {'total': len(ordered_keys), 'exact': 0, 'fallback': 0,
             'matched': 0, 'missing': 0}
    for i, key in enumerate(ordered_keys, 1):
        vals, method, matched_key, source_section = lookup_any_coeff(
            key, frc_sections, coeffs, equiv, auto_equiv, canonical)
        if vals is None:
            missing.append(key)
            stats['missing'] += 1
            count = (COEFF_VALUE_COUNTS[title] if ff_class == 'class2'
                     else OPLS_COEFF_VALUE_COUNTS[title] if ff_class == 'opls'
                     else CLASS1_COEFF_VALUE_COUNTS[title])
            if ff_class == 'class2' and title in (
                    'BondBond Coeffs', 'BondAngle Coeffs',
                    'AngleAngleTorsion Coeffs', 'EndBondTorsion Coeffs',
                    'MiddleBondTorsion Coeffs', 'BondBond13 Coeffs',
                    'AngleTorsion Coeffs', 'AngleAngle Coeffs'):
                vals = complete_class2_cross_term(
                    title, key, [], coeffs, equiv, auto_equiv)
            else:
                vals = [0.0] * count
        else:
            stats['matched'] += 1
            if method == 'exact':
                stats['exact'] += 1
            else:
                stats['fallback'] += 1
            if ff_class == 'class2' and title in (
                    'BondBond Coeffs', 'BondAngle Coeffs',
                    'AngleAngleTorsion Coeffs', 'EndBondTorsion Coeffs',
                    'MiddleBondTorsion Coeffs', 'BondBond13 Coeffs',
                    'AngleTorsion Coeffs', 'AngleAngle Coeffs'):
                vals = complete_class2_cross_term(
                    title, key, vals, coeffs, equiv, auto_equiv)
            else:
                vals = transform_coeff_values(
                    title, ff_class, vals, source_section)
            if method != 'exact':
                fallback[method].append((key, matched_key))
        extra = ''
        if method in ('equivalence', 'wildcard', 'hpan_to_hc', 'general_equiv', 'similar') and matched_key:
            extra = f' via {method}: {"-".join(matched_key)}'
        elif method == 'zero':
            extra = ' WARNING: zero fallback'
        fh.write(f'{i:<6d}{fmt_values(vals)} {comment(key, extra)}\n')
    fh.write('\n')
    return missing, fallback, stats


def write_data(out_path, source_data, frc_path, atoms, bonds, nta,
               masses, coeffs, equiv, auto_equiv, box, ff_class):
    neighbours = build_graph(bonds)
    topology = build_topology(bonds, neighbours, nta)
    if ff_class == 'opls':
        # OPLS harmonic impropers enforce planarity at three-connected
        # centers. Do not generate them for tetrahedral (four-connected)
        # centers as class-II angle-angle topology may do.
        topology['Impropers'] = [item for item in topology['Impropers']
                                 if len(neighbours.get(item[1], [])) == 3]

    atom_types = []
    for aid in sorted(atoms):
        t = nta[aid]
        if t not in atom_types:
            atom_types.append(t)
    atom_type_ids = {t: i + 1 for i, t in enumerate(atom_types)}

    bond_type_ids, bond_type_keys = unique_type_ids(
        topology['Bonds'], nta, canonical_bond)
    angle_type_ids, angle_type_keys = unique_type_ids(
        topology['Angles'], nta, canonical_angle)
    dihedral_type_ids, dihedral_type_keys = unique_type_ids(
        topology['Dihedrals'], nta, canonical_dihedral)
    improper_type_ids, improper_type_keys = unique_type_ids(
        topology['Impropers'], nta, canonical_improper)

    type_keys = {
        'Bond Coeffs': bond_type_keys,
        'Angle Coeffs': angle_type_keys,
        'Dihedral Coeffs': dihedral_type_keys,
        'Improper Coeffs': improper_type_keys,
        'BondBond Coeffs': angle_type_keys,
        'BondAngle Coeffs': angle_type_keys,
        'AngleAngleTorsion Coeffs': dihedral_type_keys,
        'EndBondTorsion Coeffs': dihedral_type_keys,
        'MiddleBondTorsion Coeffs': dihedral_type_keys,
        'BondBond13 Coeffs': dihedral_type_keys,
        'AngleTorsion Coeffs': dihedral_type_keys,
        'AngleAngle Coeffs': improper_type_keys,
    }
    canon_by_title = {
        'Bond Coeffs': canonical_bond,
        'Angle Coeffs': canonical_angle,
        'Dihedral Coeffs': canonical_dihedral,
        'Improper Coeffs': canonical_improper,
        'BondBond Coeffs': canonical_angle,
        'BondAngle Coeffs': canonical_angle,
        'AngleAngleTorsion Coeffs': canonical_dihedral,
        'EndBondTorsion Coeffs': canonical_dihedral,
        'MiddleBondTorsion Coeffs': canonical_dihedral,
        'BondBond13 Coeffs': canonical_dihedral,
        'AngleTorsion Coeffs': canonical_dihedral,
        'AngleAngle Coeffs': canonical_improper,
    }

    missing_report = defaultdict(list)
    fallback_report = defaultdict(lambda: defaultdict(list))
    stats_report = {}
    coeff_specs = (CLASS2_COEFF_SPECS if ff_class == 'class2'
                   else OPLS_COEFF_SPECS if ff_class == 'opls'
                   else CLASS1_COEFF_SPECS)

    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(
            f'HEADER,  {os.path.basename(source_data)}  converted with '
            f'build_input_lmp.py using {os.path.basename(frc_path)}  '
            f'Class: {2 if ff_class == "class2" else "OPLS" if ff_class == "opls" else 1}\n\n'
        )
        fh.write(f'{len(atoms)} atoms\n')
        fh.write(f'{len(topology["Bonds"])} bonds\n')
        fh.write(f'{len(topology["Angles"])} angles\n')
        fh.write(f'{len(topology["Dihedrals"])} dihedrals\n')
        fh.write(f'{len(topology["Impropers"])} impropers\n\n')
        fh.write(f'{len(atom_types)} atom types\n')
        fh.write(f'{len(bond_type_keys)} bond types\n')
        fh.write(f'{len(angle_type_keys)} angle types\n')
        fh.write(f'{len(dihedral_type_keys)} dihedral types\n')
        fh.write(f'{len(improper_type_keys)} improper types\n\n')

        if box:
            for line in box:
                fh.write(line + '\n')
        else:
            xs = [a['x'] for a in atoms.values()]
            ys = [a['y'] for a in atoms.values()]
            zs = [a['z'] for a in atoms.values()]
            fh.write(f'{min(xs)-0.5:.9f} {max(xs)+0.5:.9f} xlo xhi\n')
            fh.write(f'{min(ys)-0.5:.9f} {max(ys)+0.5:.9f} ylo yhi\n')
            fh.write(f'{min(zs)-0.5:.9f} {max(zs)+0.5:.9f} zlo zhi\n')
        fh.write('\n')

        fh.write('Masses\n\n')
        for t in atom_types:
            fh.write(f'{atom_type_ids[t]:<6d}{masses.get(t, 1.0):14.8f}  #  {t}\n')
        fh.write('\n')

        pair_keys = [(t,) for t in atom_types]
        for title, frc_sections, _ntypes, style in coeff_specs:
            if isinstance(frc_sections, str):
                frc_sections = (frc_sections,)
            keys = pair_keys if title == 'Pair Coeffs' else type_keys[title]
            canonical = (lambda x: x) if title == 'Pair Coeffs' else canon_by_title[title]
            missing, fallback, stats = write_coeff_section(
                fh, title, style, keys, frc_sections,
                coeffs, equiv, auto_equiv, canonical, ff_class)
            stats_report[title] = stats
            if missing:
                missing_report[title].extend(missing)
            for method, entries in fallback.items():
                fallback_report[method][title].extend(entries)

        fh.write('Atoms # full\n\n')
        for aid in sorted(atoms):
            a = atoms[aid]
            t = nta[aid]
            tail = ' '.join(a['tail'])
            tail = (' ' + tail) if tail else ''
            fh.write(
                f'{aid:6d} {a["mol"]:5d} {atom_type_ids[t]:5d} {a["charge"]:12.6f} '
                f'{a["x"]:14.9f} {a["y"]:14.9f} {a["z"]:14.9f}{tail} #  {t}\n'
            )
        fh.write('\n')

        fh.write('Bonds\n\n')
        for i, item in enumerate(topology['Bonds'], 1):
            tid = bond_type_ids[canonical_bond(atom_type_tuple(item, nta))]
            fh.write(f'{i:6d} {tid:5d} {item[0]:6d} {item[1]:6d}\n')
        fh.write('\n')

        fh.write('Angles\n\n')
        for i, item in enumerate(topology['Angles'], 1):
            tid = angle_type_ids[canonical_angle(atom_type_tuple(item, nta))]
            fh.write(f'{i:6d} {tid:5d} {item[0]:6d} {item[1]:6d} {item[2]:6d}\n')
        fh.write('\n')

        fh.write('Dihedrals\n\n')
        for i, item in enumerate(topology['Dihedrals'], 1):
            tid = dihedral_type_ids[canonical_dihedral(atom_type_tuple(item, nta))]
            fh.write(
                f'{i:6d} {tid:5d} {item[0]:6d} {item[1]:6d} '
                f'{item[2]:6d} {item[3]:6d}\n'
            )
        fh.write('\n')

        fh.write('Impropers\n\n')
        for i, item in enumerate(topology['Impropers'], 1):
            tid = improper_type_ids[canonical_improper(atom_type_tuple(item, nta))]
            fh.write(
                f'{i:6d} {tid:5d} {item[0]:6d} {item[1]:6d} '
                f'{item[2]:6d} {item[3]:6d}\n'
            )

    return missing_report, fallback_report, stats_report


def main():
    args = sys.argv[1:]
    if args:
        if len(args) != 4:
            print('Usage: python build_input_lmp.py input.data types.nta forcefield.frc output.data')
            sys.exit(1)
        data_path, nta_path, frc_path, out_path = args
    else:
        print()
        print('Build LAMMPS data from .data + .nta + .frc')
        print('-------------------------------------------')
        data_path = prompt_path('input .data file')
        nta_path = prompt_path('.nta file')
        frc_path = prompt_path('.frc force-field file')
        out_path = prompt_output()

    for path in (data_path, nta_path, frc_path):
        if not os.path.isfile(path):
            print(f'[ERROR] File not found: {path}')
            sys.exit(1)

    nta = parse_nta(nta_path)
    _header, box, atoms, bonds = parse_lammps_data(data_path)
    missing_atoms = sorted(set(atoms) - set(nta))
    if missing_atoms:
        print(f'[ERROR] Missing .nta atom types for atom IDs: {missing_atoms[:10]}')
        sys.exit(1)
    hpan_changed = apply_pcff_iff_hpan_override(frc_path, bonds, nta)
    if hpan_changed:
        print(f'\nApplied hpan override to {hpan_changed} carbon-bonded hc atom(s).')

    masses, coeffs, equiv, auto_equiv, ff_class = parse_frc(frc_path)
    print(f'\nDetected force-field class: {ff_class}')
    missing, fallback, stats = write_data(
        out_path, data_path, frc_path, atoms, bonds, nta,
        masses, coeffs, equiv, auto_equiv, box, ff_class)

    print(f'\nWrote: {out_path}')
    print('\nCoefficient match statistics:')
    print(f'  {"Section":<30} {"Total":>7} {"Exact":>7} '
          f'{"Fallback":>9} {"Missing":>9} {"Matched":>10}')
    print(f'  {"-"*30} {"-"*7} {"-"*7} {"-"*9} {"-"*9} {"-"*10}')
    grand = {'total': 0, 'exact': 0, 'fallback': 0,
             'matched': 0, 'missing': 0}
    for title, item in stats.items():
        pct = 100.0 * item['matched'] / item['total'] if item['total'] else 100.0
        print(f'  {title:<30} {item["total"]:>7d} {item["exact"]:>7d} '
              f'{item["fallback"]:>9d} {item["missing"]:>9d} {pct:>9.2f}%')
        for key in grand:
            grand[key] += item[key]
    pct = 100.0 * grand['matched'] / grand['total'] if grand['total'] else 100.0
    print(f'  {"TOTAL":<30} {grand["total"]:>7d} {grand["exact"]:>7d} '
          f'{grand["fallback"]:>9d} {grand["missing"]:>9d} {pct:>9.2f}%')
    if fallback:
        print('\nCoefficient fallback matches used:')
        for method in ('equivalence', 'wildcard', 'hpan_to_hc', 'general_equiv', 'similar'):
            sections = fallback.get(method)
            if not sections:
                continue
            total = sum(len(entries) for entries in sections.values())
            print(f'  {method}: {total} coefficient(s)')
    if missing:
        print('\nMissing coefficient entries were written as zeros:')
        for title, keys in missing.items():
            shown = ', '.join('-'.join(k) for k in keys[:8])
            more = '' if len(keys) <= 8 else f' ... ({len(keys)} total)'
            print(f'  {title}: {shown}{more}')
        print('\nCheck these against the .frc file before running LAMMPS.')


if __name__ == '__main__':
    main()

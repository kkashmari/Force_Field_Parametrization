# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0
# Required packages: Python standard library; optional Open Babel, AmberTools, or RDKit for charges
"""
Steps performed:
  1. Read a MOL2, LAMMPS data, or VASP POSCAR file
  2. Assign the FF atom type to every atom
  3. Ask which charge method to use (options depend on the FF chosen)
  4. Calculate partial charges
  5. Write both the .nta type file and the LAMMPS .data file

Charge methods from .frc file
    - Bond increments and Gasteiger : PCFF / CFF91 / COMPASS
    - Fixed charges     : CLAYFF
    - Gasteiger         : any FF
    - AM1-BCC           : GAFF / AMBER
"""

import sys
import os
import re
import glob
import math
import itertools
import subprocess
from collections import defaultdict

FF_DESCRIPTIONS = {
    'pcff':              'PCFF — Polymer Consistent Force Field',
    'pcff_iff':          'PCFF-IFF — PCFF + Interface Force Field extensions',
    'cff91':             'CFF91 — Consistent Force Field 1991',
    'compass':           'COMPASS — Condensed-phase Optimised Molecular Potentials',
    'compass_published': 'COMPASS (published parameters)',
    'clayff':            'CLAYFF — Clay minerals & interfaces',
    'oplsaa':            'OPLS-AA — Optimised Potentials for Liquid Simulations',
    'cvff':              'CVFF — Consistent Valence Force Field',
    'amber':             'AMBER — Assisted Model Building with Energy Refinement',
    'gaff':              'GAFF — General AMBER Force Field',
    'charmm':            'CHARMM — Chemistry at HARvard Macromolecular Mechanics',
    'uff':               'UFF — Universal Force Field',
}

FF_HAS_BOND_INCREMENTS = ('pcff', 'cff91', 'compass')
FF_HAS_FIXED_CHARGES = ('clayff',)

HPAN_HC_FF_FILE = 'pcff_iff_v1_5_CNT_poly_solv_hpan_hc_mod.frc'
HPAN_HC_FF_NAME = os.path.splitext(HPAN_HC_FF_FILE)[0]

# masses (LAMMPS data file)
ELEMENT_MASS = {
    'C': 12.00000, 'H':  1.00800, 'O': 15.99900, 'N': 14.00670,
    'S': 32.06400, 'F': 18.99840, 'Cl':35.45300, 'Br':79.90900,
    'Si':28.08600, 'P': 30.97380, 'Al':26.98200, 'Ca':40.08000,
    'Mg':24.30500, 'Na':22.99000, 'K': 39.10200,
}
ELEMENT_ORDER = ['C','H','O','N','S','F','Cl','Br','Si','P','Al','Ca','Mg','Na','K']

# Covalent radii (angstrom); for bonds missing from POSCAR.
COVALENT_RADIUS = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'Si': 1.11,
    'Al': 1.21, 'Ca': 1.76, 'Mg': 1.41, 'Na': 1.66, 'K': 2.03,
}

# hybridisation string (for .data comment)
SYBYL_HYB_LABEL = {
    'C.3':'Sp3','C.2':'Sp2','C.1':'Sp','C.ar':'Sp2',
    'N.3':'Sp3','N.2':'Sp2','N.1':'Sp','N.ar':'Sp2','N.pl3':'Sp2','N.4':'Sp3',
    'O.3':'Sp3','O.2':'Sp2','O.co2':'Sp2','O.ar':'Sp2',
    'S.3':'Sp3','S.2':'Sp2','S.ar':'Sp2',
    'H':'Terminal',
}

# hybridisation for scoring
SYBYL_HYB = {
    'C.3':'sp3','C.2':'sp2','C.1':'sp','C.ar':'ar',
    'N.3':'sp3','N.2':'sp2','N.1':'sp','N.ar':'ar','N.pl3':'sp2','N.4':'sp3',
    'O.3':'sp3','O.2':'sp2','O.co2':'sp2','O.ar':'ar',
    'S.3':'sp3','S.2':'sp2','S.ar':'ar','S.O':'sp2','S.O2':'sp2',
    'P.3':'sp3','H':'h',
}

def get_element(sybyl):
    return sybyl.split('.')[0]


def normalize_element(value):
    """Return a consistently capitalised element symbol."""
    m = re.match(r'\s*([A-Za-z]{1,2})', str(value))
    if not m:
        return 'X'
    symbol = m.group(1).capitalize()
    return symbol if symbol in ELEMENT_MASS else symbol[0]


def sybyl_from_element(element, coordination=0):
    """Make a conservative Sybyl-like label when the input has no such type."""
    el = normalize_element(element)
    if el == 'H': return 'H'
    if el == 'C': return 'C.3' if coordination >= 4 else ('C.2' if coordination == 3 else 'C.1')
    if el == 'N': return 'N.3' if coordination >= 3 else ('N.2' if coordination == 2 else 'N.1')
    if el == 'O': return 'O.3' if coordination >= 2 else 'O.2'
    if el == 'S': return 'S.3' if coordination >= 2 else 'S.2'
    if el == 'P': return 'P.3'
    return el


def describe_ff(stem):
    sl = stem.lower()
    for key, desc in FF_DESCRIPTIONS.items():
        if key in sl:
            return desc
    return f'{stem} force field'


def ff_family(frc_path):
    """Return the canonical FF family name for a .frc file path."""
    sl = os.path.basename(frc_path).lower()
    for key in FF_DESCRIPTIONS:
        if key in sl:
            return key
    return ''


def _word_in(text, word):
    """True if `word` appears in `text` NOT immediately after 'non'/'not'."""
    idx = 0
    while True:
        p = text.find(word, idx)
        if p == -1:
            return False
        pre = text[max(0, p - 5):p]
        if 'non' not in pre and 'not' not in pre:
            return True
        idx = p + 1


def _h_count_hint(desc):
    m = re.search(r'\b([123])\s*h(?:h|c)?(?:\'s)?(?!\w)', desc)
    return int(m.group(1)) if m else None

#  MOL2 PARSER


def parse_mol2(path):
    atoms, bonds = [], []
    section = None
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('@<TRIPOS>'):
                section = line[9:].strip(); continue
            if section == 'ATOM':
                parts = line.split()
                if len(parts) < 6: continue
                atoms.append({
                    'id':     int(parts[0]),
                    'name':   parts[1],
                    'x':      float(parts[2]),
                    'y':      float(parts[3]),
                    'z':      float(parts[4]),
                    'sybyl':  parts[5],
                    'charge': float(parts[8]) if len(parts) > 8 else 0.0,
                })
            elif section == 'BOND':
                parts = line.split()
                if len(parts) < 4: continue
                bonds.append({
                    'id':   int(parts[0]),
                    'a1':   int(parts[1]),
                    'a2':   int(parts[2]),
                    'type': parts[3],
                })
    return atoms, bonds


def _element_from_mass(mass):
    return min(ELEMENT_MASS, key=lambda el: abs(ELEMENT_MASS[el] - mass))


def parse_lammps_data(path):
    """Read common atomic/charge/molecular/full LAMMPS data-file variants."""
    with open(path, encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()
    bounds, masses, atoms, bonds = {}, {}, [], []
    section = None
    atom_style = None
    section_names = {'Masses', 'Atoms', 'Bonds', 'Angles', 'Dihedrals',
                     'Impropers', 'Velocities', 'Pair Coeffs', 'Bond Coeffs',
                     'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs'}
    for raw in lines:
        body, _, comment = raw.partition('#')
        s = body.strip()
        words = s.split()
        if len(words) >= 4 and words[-2:] in (['xlo','xhi'], ['ylo','yhi'], ['zlo','zhi']):
            bounds[words[-2][0]] = (float(words[0]), float(words[1])); continue
        if len(words) >= 6 and words[-3:] == ['xy', 'xz', 'yz']:
            bounds['tilts'] = tuple(float(v) for v in words[:3]); continue
        heading = raw.strip().split('#', 1)[0].strip()
        if heading in section_names:
            section = heading
            if heading == 'Atoms' and '#' in raw:
                atom_style = raw.split('#', 1)[1].strip().lower().split()[0]
            continue
        if not s:
            continue
        try:
            if section == 'Masses' and len(words) >= 2:
                tid, mass = int(words[0]), float(words[1])
                hint = comment.strip().split()[0] if comment.strip() else ''
                masses[tid] = normalize_element(hint) if hint else _element_from_mass(mass)
            elif section == 'Atoms':
                vals = words
                aid = int(vals[0])
                # LAMMPS styles: atomic=id type x y z; charge=id type q x y z;
                # molecular=id mol type x y z; full=id mol type q x y z.
                if atom_style == 'atomic' or (atom_style is None and len(vals) == 5):
                    tid, q, xyz = int(vals[1]), 0.0, vals[2:5]
                elif atom_style == 'charge' or (atom_style is None and len(vals) == 6):
                    tid, q, xyz = int(vals[1]), float(vals[2]), vals[3:6]
                elif atom_style == 'molecular':
                    tid, q, xyz = int(vals[2]), 0.0, vals[3:6]
                else:
                    tid, q, xyz = int(vals[2]), float(vals[3]), vals[4:7]
                hint = comment.strip().split()[0] if comment.strip() else ''
                el = normalize_element(hint) if hint else masses.get(tid, 'X')
                atoms.append({'id': aid, 'name': f'{el}{aid}', 'x': float(xyz[0]),
                              'y': float(xyz[1]), 'z': float(xyz[2]),
                              'sybyl': el, 'charge': q})
            elif section == 'Bonds' and len(words) >= 4:
                bonds.append({'id': int(words[0]), 'a1': int(words[2]),
                              'a2': int(words[3]), 'type': '1'})
        except (ValueError, IndexError):
            continue
    if not atoms:
        raise ValueError('no atoms found in LAMMPS data file')
    degree = defaultdict(int)
    for b in bonds:
        degree[b['a1']] += 1; degree[b['a2']] += 1
    for a in atoms:
        a['sybyl'] = sybyl_from_element(get_element(a['sybyl']), degree[a['id']])
    return atoms, bonds, {'bounds': bounds}


def _matvec(frac, lattice):
    return tuple(sum(frac[j] * lattice[j][i] for j in range(3)) for i in range(3))


def infer_poscar_bonds(atoms, fractional, lattice):
    bonds = []
    for i, a in enumerate(atoms):
        e1 = get_element(a['sybyl'])
        for j in range(i + 1, len(atoms)):
            e2 = get_element(atoms[j]['sybyl'])
            df = [fractional[j][k] - fractional[i][k] for k in range(3)]
            df = [v - round(v) for v in df]
            dv = _matvec(df, lattice)
            dist = math.sqrt(sum(v*v for v in dv))
            cutoff = 1.20 * (COVALENT_RADIUS.get(e1, 0.8) + COVALENT_RADIUS.get(e2, 0.8))
            if 0.1 < dist <= cutoff:
                bonds.append({'id': len(bonds)+1, 'a1': a['id'],
                              'a2': atoms[j]['id'], 'type': '1'})
    return bonds


def parse_poscar(path):
    """Read VASP 5 POSCAR/CONTCAR (Direct or Cartesian coordinates)."""
    with open(path, encoding='utf-8', errors='replace') as fh:
        rows = [line.strip() for line in fh if line.strip()]
    if len(rows) < 8:
        raise ValueError('POSCAR is incomplete')
    scale = float(rows[1].split()[0])
    lattice = [[float(x) * scale for x in rows[i].split()[:3]] for i in range(2, 5)]
    elements = [normalize_element(x) for x in rows[5].split()]
    try:
        counts = [int(x) for x in rows[6].split()]
    except ValueError as exc:
        raise ValueError('VASP 4 POSCAR needs element names on line 6') from exc
    if len(elements) != len(counts):
        raise ValueError('POSCAR element and atom-count columns differ')
    idx = 7
    if rows[idx].lower().startswith('s'):
        idx += 1
    direct = rows[idx].lower().startswith(('d', 'f'))
    idx += 1
    labels = [el for el, count in zip(elements, counts) for _ in range(count)]
    atoms, fractional = [], []
    # inverse lattice is avoided: Cartesian POSCARs are converted via a small 3x3 inverse.
    a, b, c = lattice
    det = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
           + a[2]*(b[0]*c[1]-b[1]*c[0]))
    inv = [[(b[1]*c[2]-b[2]*c[1])/det, (a[2]*c[1]-a[1]*c[2])/det, (a[1]*b[2]-a[2]*b[1])/det],
           [(b[2]*c[0]-b[0]*c[2])/det, (a[0]*c[2]-a[2]*c[0])/det, (a[2]*b[0]-a[0]*b[2])/det],
           [(b[0]*c[1]-b[1]*c[0])/det, (a[1]*c[0]-a[0]*c[1])/det, (a[0]*b[1]-a[1]*b[0])/det]]
    for n, el in enumerate(labels, 1):
        xyz = [float(x) for x in rows[idx+n-1].split()[:3]]
        if direct:
            frac, cart = xyz, _matvec(xyz, lattice)
        else:
            cart = [v * scale for v in xyz]
            frac = [sum(inv[r][k] * cart[k] for k in range(3)) for r in range(3)]
        fractional.append(frac)
        atoms.append({'id': n, 'name': f'{el}{n}', 'x': cart[0], 'y': cart[1],
                      'z': cart[2], 'sybyl': el, 'charge': 0.0})
    bonds = infer_poscar_bonds(atoms, fractional, lattice)
    # Convert the general VASP cell to LAMMPS' restricted triclinic form.
    norm = lambda v: math.sqrt(sum(x*x for x in v))
    dot = lambda u, v: sum(u[i]*v[i] for i in range(3))
    lx = norm(a); ahat = [x/lx for x in a]
    xy = dot(b, ahat); xz = dot(c, ahat)
    ly = math.sqrt(max(0.0, dot(b, b) - xy*xy))
    yz = (dot(b, c) - xy*xz) / ly
    lz = math.sqrt(max(0.0, dot(c, c) - xz*xz - yz*yz))
    restricted = [[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]]
    for atom, frac in zip(atoms, fractional):
        atom['x'], atom['y'], atom['z'] = _matvec(frac, restricted)
    degree = defaultdict(int)
    for bond in bonds:
        degree[bond['a1']] += 1; degree[bond['a2']] += 1
    for atom in atoms:
        atom['sybyl'] = sybyl_from_element(get_element(atom['sybyl']), degree[atom['id']])
    return atoms, bonds, {'lattice': lattice,
                          'bounds': {'x': (0.0, lx), 'y': (0.0, ly),
                                     'z': (0.0, lz), 'tilts': (xy, xz, yz)}}


def detect_input_format(path):
    name = os.path.basename(path).lower()
    ext = os.path.splitext(name)[1]
    if ext == '.mol2': return 'mol2'
    if ext in ('.data', '.lmp', '.lammps'): return 'lammps'
    if name in ('poscar', 'contcar') or ext in ('.vasp', '.poscar'): return 'poscar'
    raise ValueError('unsupported input format; use .mol2, .data/.lmp, POSCAR, CONTCAR, or .vasp')


def parse_structure(path):
    fmt = detect_input_format(path)
    if fmt == 'mol2':
        atoms, bonds = parse_mol2(path)
        return atoms, bonds, {}, fmt
    parser = parse_lammps_data if fmt == 'lammps' else parse_poscar
    atoms, bonds, metadata = parser(path)
    return atoms, bonds, metadata, fmt


def build_graph(bonds):
    nbrs = defaultdict(list)
    btm  = {}
    for b in bonds:
        a1, a2, bt = b['a1'], b['a2'], b['type']
        nbrs[a1].append(a2); nbrs[a2].append(a1)
        btm[(a1, a2)] = bt;  btm[(a2, a1)] = bt
    return dict(nbrs), btm


def find_small_rings(neighbours, max_size=6):
    rings, visited = [], set()
    def dfs(start, cur, path, depth):
        for nb in neighbours.get(cur, []):
            if depth > 0 and nb == start and len(path) >= 3:
                key = frozenset(path)
                if key not in visited:
                    visited.add(key); rings.append(key)
                continue
            if nb not in path and len(path) < max_size:
                dfs(start, nb, path + [nb], depth + 1)
    for aid in neighbours:
        dfs(aid, aid, [aid], 0)
    return rings


def atom_ring_sizes(neighbours):
    rings = find_small_rings(neighbours)
    rmap  = defaultdict(set)
    for ring in rings:
        for aid in ring:
            rmap[aid].add(len(ring))
    return dict(rmap)

#  FORCE-FIELD PARSER

def parse_frc(path):
    """Parse #atom_types section → list of type dicts."""
    types      = []
    seen       = set()
    in_section = False
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#atom_types'):
                in_section = True; continue
            if s.startswith('#') and in_section:
                in_section = False; continue
            if not in_section or not s or s.startswith('!') or s.startswith('>'):
                continue
            parts = s.split()
            if len(parts) < 5: continue
            try:
                float(parts[0]); int(float(parts[1]))
            except (ValueError, IndexError):
                continue
            type_name = parts[2]
            try:
                mass = float(parts[3])
            except ValueError:
                continue
            element     = parts[4]
            connections = None
            desc_start  = 5
            if len(parts) > 5:
                try:
                    connections = int(parts[5]); desc_start = 6
                except ValueError:
                    pass
            description = ' '.join(parts[desc_start:]) if len(parts) > desc_start else ''
            key = type_name.lower()
            if key in seen:
                continue
            seen.add(key)
            types.append({
                'type':        type_name,
                'element':     element,
                'connections': connections,
                'mass':        mass,
                'description': description.lower(),
                'version':     float(parts[0]),
                '_raw_desc':   description,   # keep original for charge parsing
            })
    return types


def parse_bond_increments(frc_path):
    """Return dict (type_i, type_j) -> (delta_ij, delta_ji)."""
    increments = {}
    in_bi = False
    with open(frc_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#bond_increments'):
                in_bi = True; continue
            if s.startswith('#') and in_bi:
                in_bi = False; continue
            if not in_bi or not s or s.startswith('!') or s.startswith('>'):
                continue
            parts = s.split()
            if len(parts) < 6: continue
            try:
                float(parts[0]); int(float(parts[1]))
            except (ValueError, IndexError):
                continue
            try:
                t1, t2 = parts[2], parts[3]
                d12, d21 = float(parts[4]), float(parts[5])
                increments[(t1, t2)] = (d12, d21)
                increments[(t2, t1)] = (d21, d12)
            except (ValueError, IndexError):
                continue
    return increments


def parse_auto_equivalence_bond_inct(frc_path):
    """
    Return dict: type_name -> bond_increment_key
    from the 'Bond Inct' column (col 4) of #auto_equivalence.
    """
    equiv = {}
    in_ae = False
    with open(frc_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#auto_equivalence'):
                in_ae = True; continue
            if s.startswith('#') and in_ae:
                in_ae = False; continue
            if not in_ae or not s or s.startswith('!') or s.startswith('>'):
                continue
            parts = s.split()
            if len(parts) < 5: continue
            try:
                float(parts[0]); int(float(parts[1]))
            except (ValueError, IndexError):
                continue
            equiv[parts[2]] = parts[4]   # col 4 = Bond Inct
    return equiv


def parse_fixed_charges(ff_types):
    """
    For CLAYFF-style FFs: extract the charge value from the end of
    each atom-type description string.
    Returns dict: type_name -> charge (float)
    """
    charges = {}
    for t in ff_types:
        desc = t['_raw_desc']
        m = re.search(r'([-+]?\d+\.\d+)\s*$', desc)
        if m:
            try:
                charges[t['type']] = float(m.group(1))
            except ValueError:
                pass
    return charges

#  FF TYPE FEATURE PRECOMPUTE  (for scoring)

def precompute_ff_features(ff_types):
    for t in ff_types:
        d  = t['description']
        n  = t['type'].lower()
        el = t['element']
        c  = t['connections']

        hyb = None; hyb_explicit = False
        if 'sp3' in d:            hyb = 'sp3'; hyb_explicit = True
        elif 'sp2' in d:          hyb = 'sp2'; hyb_explicit = True
        elif any(k in d for k in ('triple bond','triple-bond','sp carbon','sp nitrogen')):
            hyb = 'sp'
        elif _word_in(d, 'aromatic'):   hyb = 'ar'
        elif 'doubly bonded' in d or 'double bond' in d: hyb = 'sp2'
        elif 'carbonyl' in d:    hyb = 'sp2'
        if hyb is None and c is not None:
            if   el in ('C','N','S') and c == 4: hyb = 'sp3'
            elif el in ('C','N','S') and c == 3: hyb = 'sp2'
            elif el in ('C','N','S') and c == 2: hyb = 'sp2'
            elif el in ('C','N','S') and c == 1: hyb = 'sp'
            elif el == 'O' and c == 2:           hyb = 'sp3'
            elif el == 'O' and c == 1:           hyb = 'sp2'
            elif el == 'H':                      hyb = 'h'
        t['hyb'] = hyb; t['hyb_explicit'] = hyb_explicit

        t['ring3'] = ('3-membered' in d or '3- membered' in d or
                      'three membered' in d or 'three-membered' in d or
                      bool(re.search(r'\b(c3m|c3h|n3m|s3e|o3e)\b', n)))
        t['ring4'] = ('4-membered' in d or '4- membered' in d or
                      'four membered' in d or 'four-membered' in d or
                      bool(re.search(r'\b(c4m|c4h|n4m|s4e|o4e)\b', n)))
        t['ring5'] = '5-membered' in d or '5 membered' in d
        t['ring6'] = '6-membered' in d or '6 membered' in d
        t['has_H_in_ring'] = (t['ring3'] or t['ring4']) and \
                              ('with h' in d or 'with hydrogen' in d or n.endswith('h'))

        t['is_carbonate']         = 'carbonate' in d
        t['is_inorganic_mineral'] = el == 'O' and (
            n.startswith('oc') or
            any(k in d for k in ('silicate','ettringite','tobermorite',
                                  'aluminate','clay','gypsum','apatite',
                                  'monosulfate','sulfate unit')))

        t['is_carbonyl_O']         = el == 'O' and ('carbonyl' in d or (c is not None and c == 1))
        t['is_organic_carbonyl_O'] = el == 'O' and 'carbonyl group' in d
        t['is_ester_carbonyl_O']   = (el == 'O' and 'ester' in d and
                                       ('c=o' in d or 'carbonyl' in d))
        t['is_ester_O']            = (el == 'O' and 'ester' in d and
                                      'carbonate' not in d and
                                      not t['is_ester_carbonyl_O'])
        t['is_hydroxyl_O']         = el == 'O' and any(
            k in d for k in ('hydroxyl', 'bonded to hydrogen', '-oh'))
        t['is_ether_O']            = (el == 'O' and ('ether' in d or 'acetal' in d)
                                      and not t['is_hydroxyl_O'])
        t['is_water_O']            = el == 'O' and 'water' in d

        t['is_carbonyl_C']  = el == 'C' and any(
            k in d for k in ('carbonyl', 'c=o', '-coor'))
        t['is_ester_C']     = (el == 'C' and 'ester' in d and
                                'carbonate' not in d and 'och' not in d)
        t['is_carbamate_C'] = el == 'C' and ('carbamate' in d or 'urea' in d or 'urethane' in d)
        t['is_amide_C']     = el == 'C' and ('amide' in d or 'amid' in d)
        t['is_aromatic_C']  = el == 'C' and _word_in(d, 'aromatic')
        t['is_db_terminal'] = el == 'C' and ('end doubly bonded' in d or
                                              ('terminal' in d and 'doubly' in d))
        t['is_db_internal'] = (el == 'C' and 'doubly bonded carbon' in d
                                and 'end' not in d and 'terminal' not in d)
        t['is_db_next']     = el == 'C' and 'next to' in d and 'doubly bonded' in d

        t['bonded_to_c'] = el == 'H' and 'bonded to c' in d
        t['bonded_to_n'] = el == 'H' and ('bonded to n' in d or 'bonded to nitrogen' in d)
        t['bonded_to_o'] = el == 'H' and ('bonded to o' in d or 'bonded to oxygen' in d)
        t['bonded_to_s'] = el == 'H' and ('bonded to s' in d or 'bonded to sulfur' in d)
        t['bonded_to_p'] = el == 'H' and 'bonded to p' in d
        t['is_polar_H']  = el == 'H' and any(k in d for k in ('nitrogen','oxygen'))
        t['is_water_H']  = el == 'H' and 'water' in d

        t['h_count_hint']     = _h_count_hint(d)
        t['is_charged']       = ('charged' in d or 'ion' in d
                                  or n.endswith('+') or n.endswith('-'))
        t['is_generic']       = ('generic' in d
                                  or n in ('c','h','n','o','s','p','f','cl','br'))
        t['next_to_carbonyl'] = 'next to' in d and ('carbonyl' in d or 'c=o' in d)
        t['is_united_atom'] = ('(ua)' in d or 'united atom' in d)

        # Tinker OPLS descriptions encode chemical context in compact strings.
        t['is_alkane_C'] = el == 'C' and 'alkane' in d
        t['is_alkene_C'] = el == 'C' and ('alkene' in d or 'vinyl' in d)
        t['is_ether_C'] = el == 'C' and 'ether' in d
        t['is_ester_alkyl_C'] = (el == 'C' and 'ester' in d and
                                  any(k in d for k in ('och3', 'och2', 'ochr', 'ocr3')))
        if el == 'C' and t['h_count_hint'] is None:
            if 'ch3' in d: t['h_count_hint'] = 3
            elif 'ch2' in d: t['h_count_hint'] = 2
            elif re.search(r'(?<![a-z])h2-c=', d): t['h_count_hint'] = 2
            elif re.search(r'(?<![a-z])rh-c=', d): t['h_count_hint'] = 1
            elif re.search(r'(?<![a-z])r2-c=', d): t['h_count_hint'] = 0

        t['is_alkane_H'] = el == 'H' and ('alkane h-c' in d or 'aliphatic h' in d)
        t['is_alkene_H'] = el == 'H' and ('alkene h-c=' in d or 'vinyl h' in d)
        t['is_ether_H'] = el == 'H' and 'ether h-c' in d
        t['is_ester_alkyl_H'] = el == 'H' and 'ester' in d and 'och' in d
        t['is_hetero_H'] = el == 'H' and any(
            k in d for k in ('ammonia', 'amine', 'amide', 'thiol', 'hydroxyl',
                             'alcohol', 'water', 'phenol'))
        t['is_aromatic_ether_O'] = el == 'O' and any(
            k in d for k in ('diphenyl ether', 'phenyl ether', 'aryl ether'))


#  ATOM FEATURE EXTRACTION

def extract_atom_features(atom, atoms_by_id, neighbours, btm, ring_sizes):
    aid    = atom['id']
    sybyl  = atom['sybyl']
    el     = get_element(sybyl)
    nbs    = neighbours.get(aid, [])
    rsizes = ring_sizes.get(aid, set())
    nb_els = [get_element(atoms_by_id[n]['sybyl']) for n in nbs]
    h_cnt  = nb_els.count('H')
    heavy  = [n for n in nbs if get_element(atoms_by_id[n]['sybyl']) != 'H']
    hyb    = SYBYL_HYB.get(sybyl, 'sp3' if len(nbs) >= 4 else 'sp2')

    def bt(nb):
        return btm.get((aid, nb), '1')

    double_to = {get_element(atoms_by_id[nb]['sybyl'])
                 for nb in nbs if bt(nb) in ('2','ar','am')}
    cc_double_nbs = [nb for nb in nbs
                     if get_element(atoms_by_id[nb]['sybyl']) == 'C'
                     and bt(nb) in ('2','ar','am')]

    def is_carbonyl_c(cid):
        return any(
            get_element(atoms_by_id[nb2]['sybyl']) == 'O'
            and btm.get((cid, nb2), '1') in ('2','ar','am')
            for nb2 in neighbours.get(cid, [])
        )

    is_ester_O = False
    if el == 'O' and hyb == 'sp3':
        for nb in nbs:
            if get_element(atoms_by_id[nb]['sybyl']) == 'C' and is_carbonyl_c(nb):
                is_ester_O = True; break

    o_on_parent_c = 0
    if el == 'O' and nbs:
        pc = [n for n in nbs if get_element(atoms_by_id[n]['sybyl']) == 'C']
        if pc:
            o_on_parent_c = sum(1 for nb2 in neighbours.get(pc[0],[])
                                if get_element(atoms_by_id[nb2]['sybyl']) == 'O')

    is_carbonyl_O = (sybyl in ('O.2','O.co2') or
                     (el == 'O' and 'C' in double_to and len(nbs) == 1))
    is_ester_carbonyl_O = False
    if is_carbonyl_O and nbs:
        parent = nbs[0]
        is_ester_carbonyl_O = sum(
            1 for n in neighbours.get(parent, [])
            if get_element(atoms_by_id[n]['sybyl']) == 'O') >= 2
    is_ether_O    = el == 'O' and hyb == 'sp3' and 'H' not in nb_els and not is_ester_O
    ether_has_aromatic_neighbor = el == 'O' and any(
        atoms_by_id[n]['sybyl'] in ('C.ar', 'N.ar') for n in nbs)
    is_hydroxyl_O = el == 'O' and 'H' in nb_els
    is_water      = el == 'O' and set(nb_els) == {'H'} and h_cnt == 2

    has_cc_double = el == 'C' and 'C' in double_to
    has_co_double = el == 'C' and 'O' in double_to

    o_neighbor_count = sum(1 for e in nb_els if e == 'O') if el == 'C' else 0
    is_carbonate_C   = el == 'C' and o_neighbor_count >= 3

    heavy_outside_cc_double = [n for n in heavy if n not in cc_double_nbs]
    is_db_terminal = (el=='C' and hyb=='sp2' and has_cc_double
                      and not heavy_outside_cc_double and not has_co_double)
    is_db_internal = (el=='C' and hyb=='sp2' and has_cc_double
                      and bool(heavy_outside_cc_double) and not has_co_double)
    is_db_next     = (el=='C' and hyb=='sp2' and not has_cc_double and not has_co_double
                      and any(atoms_by_id[n]['sybyl']=='C.2' for n in nbs))
    is_amide_C     = el=='C' and has_co_double and 'N' in set(nb_els)
    is_carbamate_C = (el=='C' and has_co_double
                      and sum(1 for e in nb_els if e=='N') >= 1
                      and sum(1 for e in nb_els if e=='O') >= 2)
    next_to_carbonyl = el=='C' and any(
        is_carbonyl_c(n) for n in nbs
        if get_element(atoms_by_id[n]['sybyl']) == 'C'
    )
    is_ester_alkyl_C = False
    if el == 'C' and hyb == 'sp3':
        for oid in nbs:
            if get_element(atoms_by_id[oid]['sybyl']) != 'O':
                continue
            if any(other != aid and get_element(atoms_by_id[other]['sybyl']) == 'C'
                   and is_carbonyl_c(other) for other in neighbours.get(oid, [])):
                is_ester_alkyl_C = True
                break
    H_bonded_to = (get_element(atoms_by_id[nbs[0]]['sybyl'])
                   if el == 'H' and nbs else None)
    h_parent_hyb = h_parent_is_ether = h_parent_is_ester_alkyl = False
    if el == 'H' and nbs:
        parent = nbs[0]
        parent_atom = atoms_by_id[parent]
        h_parent_hyb = SYBYL_HYB.get(parent_atom['sybyl'])
        parent_nbs = neighbours.get(parent, [])
        oxygen_nbs = [n for n in parent_nbs
                      if get_element(atoms_by_id[n]['sybyl']) == 'O']
        h_parent_is_ether = bool(oxygen_nbs)
        for oid in oxygen_nbs:
            for other in neighbours.get(oid, []):
                if other != parent and get_element(atoms_by_id[other]['sybyl']) == 'C' \
                        and is_carbonyl_c(other):
                    h_parent_is_ester_alkyl = True

    return {
        'element':          el, 'sybyl': sybyl, 'hyb': hyb,
        'valence':          len(nbs), 'h_count': h_cnt, 'heavy_count': len(heavy),
        'double_to':        double_to, 'has_cc_double': has_cc_double,
        'has_co_double':    has_co_double,
        'in_ring3':         3 in rsizes, 'in_ring4': 4 in rsizes,
        'in_ring5':         5 in rsizes, 'in_ring6': 6 in rsizes,
        'is_ester_O':       is_ester_O, 'is_ether_O': is_ether_O,
        'is_ester_carbonyl_O': is_ester_carbonyl_O,
        'ether_has_aromatic_neighbor': ether_has_aromatic_neighbor,
        'is_carbonyl_O':    is_carbonyl_O, 'is_hydroxyl_O': is_hydroxyl_O,
        'is_water':         is_water, 'o_on_parent_c': o_on_parent_c,
        'o_neighbor_count': o_neighbor_count, 'is_carbonate_C': is_carbonate_C,
        'is_db_terminal':   is_db_terminal, 'is_db_internal': is_db_internal,
        'is_db_next':       is_db_next, 'is_amide_C': is_amide_C,
        'is_carbamate_C':   is_carbamate_C, 'next_to_carbonyl': next_to_carbonyl,
        'is_ester_alkyl_C': is_ester_alkyl_C,
        'H_bonded_to':      H_bonded_to,
        'h_parent_hyb': h_parent_hyb,
        'h_parent_is_ether': h_parent_is_ether,
        'h_parent_is_ester_alkyl': h_parent_is_ester_alkyl,
    }

#  SCORING ENGINE

def score_type(ff_type, feat):
    t  = ff_type
    el = feat['element']
    if t['element'].upper() != el.upper():
        return None
    if t['connections'] is not None and t['connections'] > 0:
        if abs(t['connections'] - feat['valence']) > 1:
            return None
    # Explicit-H MOL2/LAMMPS structures must not receive united-atom types.
    if t.get('is_united_atom'):
        return None

    hyb   = feat['hyb']
    t_hyb = t['hyb']
    score = 0; spec = 0

    if t_hyb and hyb:
        if t_hyb == hyb:
            score += 30; spec += 1
            if t['hyb_explicit']: score += 10; spec += 1
        else:
            if   {t_hyb,hyb} >= {'sp3','sp2'}: score -= 25
            elif {t_hyb,hyb} >= {'sp3','sp'}:  score -= 20
            elif {t_hyb,hyb} >= {'sp2','sp'}:  score -= 10
            elif 'ar' in (t_hyb,hyb):          score -= 20

    for rn, fkey, pts in [(3,'in_ring3',25),(4,'in_ring4',25),
                           (5,'in_ring5',10),(6,'in_ring6',5)]:
        if t[f'ring{rn}']:
            if feat[fkey]:  score += pts; spec += (3 if rn<=4 else 1)
            else:           score -= pts
    if t['has_H_in_ring'] and (feat['in_ring3'] or feat['in_ring4']):
        score += (10 if feat['h_count']>0 else -10)
        if feat['h_count']>0: spec += 1

    if t['h_count_hint'] is not None:
        if t['h_count_hint'] == feat['h_count']: score += 25; spec += 3
        else:                                     score -= 15

    if t['is_carbonate']:
        hit = feat['is_carbonate_C'] if el=='C' else feat['o_on_parent_c']>=3
        score += (20 if hit else -20); spec += (2 if hit else 0)
    if t['is_inorganic_mineral']:
        score -= 15

    if el == 'O':
        if t.get('is_aromatic_ether_O'):
            hit = feat['ether_has_aromatic_neighbor']
            score += (20 if hit else -20); spec += (2 if hit else 0)
        if t.get('is_ester_carbonyl_O'):
            hit = feat['is_ester_carbonyl_O']
            score += (30 if hit else -25); spec += (4 if hit else 0)
        if t['is_carbonyl_O']:
            hit = feat['is_carbonyl_O']
            score += (30 if hit else -30); spec += (3 if hit else 0)
        if t['is_organic_carbonyl_O']:
            hit = feat['is_carbonyl_O']
            score += (15 if hit else -5);  spec += (2 if hit else 0)
        if t['is_ester_O']:
            hit = feat['is_ester_O']
            score += (25 if hit else -20); spec += (3 if hit else 0)
        if t['is_ether_O']:
            hit = feat['is_ether_O']
            score += (20 if hit else -15); spec += (2 if hit else 0)
        if t['is_hydroxyl_O']:
            hit = feat['is_hydroxyl_O']
            score += (20 if hit else -10); spec += (2 if hit else 0)
        if t['is_water_O']:
            hit = feat['is_water']
            score += (20 if hit else -10); spec += (2 if hit else 0)

    if el == 'C':
        if t.get('is_alkane_C'):
            score += (15 if hyb == 'sp3' else -20)
        if t.get('is_alkene_C'):
            score += (20 if feat['has_cc_double'] else -25)
            spec += (2 if feat['has_cc_double'] else 0)
        if t.get('is_ether_C'):
            hit = feat['hyb'] == 'sp3' and feat['o_neighbor_count'] > 0
            score += (15 if hit else -15); spec += (2 if hit else 0)
        if t.get('is_ester_alkyl_C'):
            hit = feat['is_ester_alkyl_C']
            score += (20 if hit else -20); spec += (2 if hit else 0)
        if t['is_db_terminal']:
            hit = feat['is_db_terminal']
            score += (25 if hit else -20); spec += (3 if hit else 0)
        if t['is_db_internal']:
            hit = feat['is_db_internal']
            score += (25 if hit else -20); spec += (3 if hit else 0)
        if t['is_db_next']:
            hit = feat['is_db_next']
            score += (15 if hit else -10); spec += (2 if hit else 0)
        if t['is_carbonyl_C']:
            hit = feat['has_co_double']
            score += (20 if hit else -20); spec += (2 if hit else 0)
        if t['is_ester_C']:
            hit = feat['has_co_double']
            score += (30 if hit else -20); spec += (3 if hit else 0)
        if t['is_carbamate_C']:
            hit = feat['is_carbamate_C']
            score += (20 if hit else -10); spec += (2 if hit else 0)
        if t['is_amide_C']:
            hit = feat['is_amide_C']
            score += (10 if hit else -5);  spec += (1 if hit else 0)
        if t['is_aromatic_C']:
            hit = hyb == 'ar'
            score += (15 if hit else -15); spec += (1 if hit else 0)
        if t['next_to_carbonyl']:
            hit = feat['next_to_carbonyl']
            score += (10 if hit else -5);  spec += (1 if hit else 0)

    if el == 'H':
        nb = feat['H_bonded_to']
        if t.get('is_hetero_H') and nb == 'C':
            score -= 30
        if t.get('is_alkene_H'):
            hit = nb == 'C' and feat['h_parent_hyb'] == 'sp2'
            score += (25 if hit else -20); spec += (3 if hit else 0)
        if t.get('is_alkane_H'):
            hit = nb == 'C' and feat['h_parent_hyb'] == 'sp3' \
                and not feat['h_parent_is_ether']
            score += (20 if hit else -15); spec += (2 if hit else 0)
        if t.get('is_ether_H'):
            hit = nb == 'C' and feat['h_parent_is_ether']
            score += (20 if hit else -15); spec += (2 if hit else 0)
        if t.get('is_ester_alkyl_H'):
            hit = nb == 'C' and feat['h_parent_is_ester_alkyl']
            score += (25 if hit else -20); spec += (3 if hit else 0)
        for flag, target in [('bonded_to_c','C'),('bonded_to_n','N'),
                              ('bonded_to_o','O'),('bonded_to_s','S'),
                              ('bonded_to_p','P')]:
            if t[flag]:
                score += (25 if nb==target else -20)
                spec  += (2  if nb==target else 0)
        if t['is_polar_H']:
            score += (10 if nb in ('N','O','F') else -10)
        if t['is_water_H']:
            hit = feat['is_water']
            score += (15 if hit else -10); spec += (1 if hit else 0)

    if t['is_charged']: score -= 5
    if t['is_generic']: score -= 3
    return score, spec


def assign_types(atoms, atoms_by_id, neighbours, btm, ring_sizes, ff_types):
    results = {}
    for atom in atoms:
        feat = extract_atom_features(atom, atoms_by_id, neighbours, btm, ring_sizes)
        best_type = None; best_score = None; best_spec = -1
        for t in ff_types:
            r = score_type(t, feat)
            if r is None: continue
            sc, sp = r
            if (best_score is None or sc > best_score
                    or (sc == best_score and sp > best_spec)):
                best_score = sc; best_spec = sp; best_type = t['type']
        results[atom['id']] = best_type or feat['element'].lower()
    return results


def _load_lammps_parameter_matcher():
    """Load the colocated builder so atom typing and coefficient lookup agree."""
    import importlib.util

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'build_input_lmp.py')
    if not os.path.isfile(path):
        return None
    spec = importlib.util.spec_from_file_location('_atom_typing_lmp_builder', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def optimize_types_for_parameter_coverage(
        atoms, bonds, atoms_by_id, neighbours, btm, ring_sizes, ff_types,
        type_assignments, frc_path, chemistry_window=6, report_path=None):
    """Prefer chemically valid atom types with the best FF parameter coverage.

    Chemical perception remains the gatekeeper: only candidates within
    ``chemistry_window`` points of the best chemical score are considered.
    Among those candidates, the routine uses the exact same equivalence and
    wildcard matching as build_input_lmp.py and maximises available primary
    bonded parameters first, followed by class-II cross terms.
    """
    try:
        matcher = _load_lammps_parameter_matcher()
        if matcher is None:
            return 0
        _masses, coeffs, equiv, auto_equiv, ff_class = matcher.parse_frc(frc_path)
    except (OSError, ImportError, AttributeError, ValueError) as exc:
        print(f'  [Typing] Parameter-coverage ranking skipped: {exc}')
        return 0

    topology = matcher.build_topology(bonds, neighbours, type_assignments)
    interactions_by_atom = defaultdict(list)
    for kind, items in topology.items():
        for item in items:
            for aid in item:
                interactions_by_atom[aid].append((kind, item))

    primary_title = {
        'Bonds': 'Bond Coeffs', 'Angles': 'Angle Coeffs',
        'Dihedrals': 'Dihedral Coeffs', 'Impropers': 'Improper Coeffs',
    }
    cross_titles = {
        'Bonds': (),
        'Angles': ('BondBond Coeffs', 'BondAngle Coeffs'),
        'Dihedrals': ('AngleAngleTorsion Coeffs', 'EndBondTorsion Coeffs',
                      'MiddleBondTorsion Coeffs', 'BondBond13 Coeffs',
                      'AngleTorsion Coeffs'),
        'Impropers': ('AngleAngle Coeffs',),
    }
    canonical = {
        'Bonds': matcher.canonical_bond,
        'Angles': matcher.canonical_angle,
        'Dihedrals': matcher.canonical_dihedral,
        'Impropers': matcher.canonical_improper,
    }
    specs = matcher.CLASS2_COEFF_SPECS if ff_class == 'class2' \
        else matcher.CLASS1_COEFF_SPECS
    sections_by_title = {}
    for title, sections, _ntypes, _style in specs:
        sections_by_title[title] = (sections,) if isinstance(sections, str) else sections

    lookup_cache = {}

    def has_coeff(title, key, canon):
        cache_key = (title, tuple(key))
        if cache_key not in lookup_cache:
            sections = sections_by_title.get(title, ())
            if not sections:
                lookup_cache[cache_key] = False
            else:
                vals, _method, _matched, _source = matcher.lookup_any_coeff(
                    key, sections, coeffs, equiv, auto_equiv, canon)
                lookup_cache[cache_key] = vals is not None
        return lookup_cache[cache_key]

    def local_coverage(aid, candidate):
        old = type_assignments[aid]
        type_assignments[aid] = candidate
        primary = cross = 0
        for kind, item in interactions_by_atom.get(aid, []):
            key = canonical[kind](tuple(type_assignments[i] for i in item))
            if has_coeff(primary_title[kind], key, canonical[kind]):
                primary += 1
            for title in cross_titles[kind]:
                if has_coeff(title, key, canonical[kind]):
                    cross += 1
        type_assignments[aid] = old
        return primary, cross

    report = open(report_path, 'w', encoding='utf-8') if report_path else None
    if report:
        report.write(
            '# Atom-type scoring and force-field parameter coverage\n'
            f'# Force field: {os.path.basename(frc_path)}\n'
            f'# Chemical eligibility window: best score - {chemistry_window}\n'
            '# Only chemically eligible candidates are listed.\n'
            '# primary/cross are counts for local interactions containing the atom\n'
            'atom_id\tatom_name\tsybyl\tcandidate\tchemical_score\t'
            'specificity\tprimary_coeffs\tcross_coeffs\teligible\tselected\n')

    changes = 0
    for atom in atoms:
        aid = atom['id']
        # Preserve explicit project-specific overrides applied before this pass.
        locked_override = type_assignments.get(aid) == 'hpan'
        feat = extract_atom_features(atom, atoms_by_id, neighbours, btm, ring_sizes)
        candidates = []
        for ff_type in ff_types:
            result = score_type(ff_type, feat)
            if result is not None:
                candidates.append((ff_type['type'], result[0], result[1]))
        if not candidates:
            continue
        best_chemistry = max(score for _name, score, _spec in candidates)
        eligible_names = {
            name for name, score, _spec in candidates
            if score >= best_chemistry - chemistry_window
        }
        candidates = [item for item in candidates if item[0] in eligible_names]

        ranked_all = []
        for name, chemical_score, specificity in candidates:
            primary, cross = local_coverage(aid, name)
            ranked_all.append(
                (primary, cross, chemical_score, specificity, name))

        old = type_assignments[aid]
        eligible_ranked = ranked_all
        # Prefer the existing selection when every meaningful score is tied.
        winner = old if locked_override else max(
            eligible_ranked,
            key=lambda item: (item[0], item[1], item[2], item[3],
                              item[-1] == old))[-1]
        if winner != old:
            old_cov = next((r[:2] for r in ranked_all if r[-1] == old), (0, 0))
            new_cov = next(r[:2] for r in ranked_all if r[-1] == winner)
            type_assignments[aid] = winner
            changes += 1
            print(f'  [Typing] Atom {aid}: {old} -> {winner} '
                  f'(local primary/cross coverage {old_cov[0]}/{old_cov[1]} '
                  f'-> {new_cov[0]}/{new_cov[1]})')

        if report:
            for primary, cross, chemical_score, specificity, name in sorted(
                    ranked_all, key=lambda item: (-item[2], -item[3], item[-1])):
                report.write(
                    f'{aid}\t{atom["name"]}\t{atom["sybyl"]}\t{name}\t'
                    f'{chemical_score}\t{specificity}\t{primary}\t{cross}\t'
                    f'{"yes" if name in eligible_names else "no"}\t'
                    f'{"yes" if name == winner else "no"}\n')
            report.write('\n')
    if report:
        report.close()
    return changes


def apply_pcff_iff_hpan_override(frc_path, atoms, atoms_by_id, neighbours, type_assignments):
    if os.path.basename(frc_path) != HPAN_HC_FF_FILE:
        return
    for atom in atoms:
        aid = atom['id']
        if get_element(atom['sybyl']) != 'H':
            continue
        if any(get_element(atoms_by_id[nb]['sybyl']) == 'C'
               for nb in neighbours.get(aid, [])):
            type_assignments[aid] = 'hpan'


#  CHARGE CALCULATION METHODS

def _rdkit_gasteiger(atoms, bonds):
    """Calculate Gasteiger charges from the parsed molecular graph with RDKit."""
    from rdkit import Chem
    from rdkit.Chem import rdPartialCharges

    mol = Chem.RWMol()
    id_to_idx = {}
    for atom in atoms:
        rd_atom = Chem.Atom(get_element(atom['sybyl']))
        # These are the two common charged Sybyl types for which the formal
        # charge is unambiguous.  Other types are left neutral and inferred.
        if atom['sybyl'] == 'N.4':
            rd_atom.SetFormalCharge(1)
        id_to_idx[atom['id']] = mol.AddAtom(rd_atom)

    bond_types = {
        '1': Chem.BondType.SINGLE, 'single': Chem.BondType.SINGLE,
        '2': Chem.BondType.DOUBLE, 'double': Chem.BondType.DOUBLE,
        '3': Chem.BondType.TRIPLE, 'triple': Chem.BondType.TRIPLE,
        'ar': Chem.BondType.AROMATIC, 'aro': Chem.BondType.AROMATIC,
    }
    for bond in bonds:
        bt = bond_types.get(str(bond.get('type', '1')).lower(),
                            Chem.BondType.SINGLE)
        mol.AddBond(id_to_idx[bond['a1']], id_to_idx[bond['a2']], bt)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    rdPartialCharges.ComputeGasteigerCharges(mol, nIter=12,
                                             throwOnParamFailure=True)
    charges = {}
    for atom in atoms:
        value = float(mol.GetAtomWithIdx(id_to_idx[atom['id']])
                      .GetProp('_GasteigerCharge'))
        if not math.isfinite(value):
            raise ValueError(f'non-finite charge for atom {atom["id"]}')
        charges[atom['id']] = value
    return charges


def _write_temporary_mol2(path, atoms, bonds):
    """Write the parsed graph as MOL2 for the Open Babel fallback."""
    atom_ids = {a['id'] for a in atoms}
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('@<TRIPOS>MOLECULE\nstructure\n')
        fh.write(f'{len(atoms)} {len(bonds)} 0 0 0\nSMALL\nNO_CHARGES\n\n')
        fh.write('@<TRIPOS>ATOM\n')
        for atom in atoms:
            fh.write(f'{atom["id"]:7d} {atom["name"]:<8} '
                     f'{atom["x"]:12.6f} {atom["y"]:12.6f} {atom["z"]:12.6f} '
                     f'{atom["sybyl"]:<7} 1 MOL 0.000000\n')
        fh.write('@<TRIPOS>BOND\n')
        for i, bond in enumerate(bonds, 1):
            if bond['a1'] in atom_ids and bond['a2'] in atom_ids:
                fh.write(f'{i:6d} {bond["a1"]:6d} {bond["a2"]:6d} '
                         f'{bond.get("type", "1")}\n')


def _openbabel_gasteiger(atoms, bonds):
    """Calculate Gasteiger charges with the Open Babel command-line tool."""
    import tempfile
    import shutil

    obabel = shutil.which('obabel')
    if not obabel:
        raise RuntimeError('Open Babel executable (obabel) is not installed')
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, 'input.mol2')
        out_path = os.path.join(tmpdir, 'charged.mol2')
        _write_temporary_mol2(in_path, atoms, bonds)
        result = subprocess.run(
            [obabel, in_path, '-O', out_path, '--partialcharge', 'gasteiger'],
            capture_output=True, text=True)
        if result.returncode != 0 or not os.path.isfile(out_path):
            raise RuntimeError(result.stderr.strip() or 'Open Babel failed')
        charged_atoms, _ = parse_mol2(out_path)
        if len(charged_atoms) != len(atoms):
            raise RuntimeError('Open Babel output has a different atom count')
        return {source['id']: charged['charge']
                for source, charged in zip(atoms, charged_atoms)}


def _antechamber_gasteiger(atoms, bonds):
    """Calculate Gasteiger charges with AmberTools/antechamber."""
    import tempfile
    import shutil

    antechamber = shutil.which('antechamber')
    if not antechamber:
        raise RuntimeError('AmberTools executable (antechamber) is not installed')
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, 'input.mol2')
        out_path = os.path.join(tmpdir, 'charged.mol2')
        _write_temporary_mol2(in_path, atoms, bonds)
        result = subprocess.run(
            [antechamber, '-i', in_path, '-fi', 'mol2', '-o', out_path,
             '-fo', 'mol2', '-c', 'gas', '-s', '0'],
            capture_output=True, text=True, cwd=tmpdir)
        if result.returncode != 0 or not os.path.isfile(out_path):
            raise RuntimeError(result.stderr.strip() or 'antechamber failed')
        charged_atoms, _ = parse_mol2(out_path)
        if len(charged_atoms) != len(atoms):
            raise RuntimeError('antechamber output has a different atom count')
        return {source['id']: charged['charge']
                for source, charged in zip(atoms, charged_atoms)}


def charges_gasteiger(atoms, bonds, **_):
    """Calculate Gasteiger-Marsili charges, independent of force-field type."""
    try:
        charges = _rdkit_gasteiger(atoms, bonds)
        print('  [Charges] Gasteiger charges calculated with RDKit.')
        return charges
    except Exception as rdkit_error:
        try:
            charges = _antechamber_gasteiger(atoms, bonds)
            print('  [Charges] Gasteiger charges calculated with AmberTools.')
            return charges
        except Exception as amber_error:
            # Open Babel also handles some unusual valence patterns leniently.
            try:
                charges = _openbabel_gasteiger(atoms, bonds)
                print('  [Charges] Gasteiger charges calculated with Open Babel.')
                return charges
            except Exception as obabel_error:
                raise RuntimeError(
                    'Could not calculate Gasteiger charges. Install RDKit, '
                    'AmberTools, or Open Babel from conda-forge. '
                    f'RDKit: {rdkit_error}; AmberTools: {amber_error}; '
                    f'Open Babel: {obabel_error}')


def charges_input(atoms, **_):
    """Keep charges supplied by the input (zero for POSCAR)."""
    return {a['id']: a.get('charge', 0.0) for a in atoms}


def charges_bond_increments(atoms, bonds, type_assignments, frc_path, **_):
    """
    Compute PCFF/CFF91/COMPASS bond-increment charges.
    charge(i) = Σ_bonds  Δ(bond_equiv(type_i) → bond_equiv(type_j))
    """
    increments = parse_bond_increments(frc_path)
    bond_equiv = parse_auto_equivalence_bond_inct(frc_path)

    def beq(t):
        return bond_equiv.get(t, t)

    charges = {a['id']: 0.0 for a in atoms}
    missing = set()

    for b in bonds:
        a1, a2 = b['a1'], b['a2']
        t1 = beq(type_assignments[a1])
        t2 = beq(type_assignments[a2])
        key = (t1, t2)
        if key in increments:
            d12, d21 = increments[key]
            charges[a1] += d12
            charges[a2] += d21
        else:
            missing.add((t1, t2))

    if missing:
        print(f'\n  [Charges] {len(missing)} bond-increment pair(s) not found in .frc '
              f'(set to 0): {", ".join(f"{a}-{b}" for a,b in sorted(missing)[:5])}')

    return charges


def charges_fixed(atoms, type_assignments, ff_types, **_):
    """
    Return fixed charges from the FF atom_types description (CLAYFF style).
    """
    fixed = parse_fixed_charges(ff_types)
    charges = {}
    missing = []
    for a in atoms:
        t = type_assignments[a['id']]
        if t in fixed:
            charges[a['id']] = fixed[t]
        else:
            charges[a['id']] = 0.0
            missing.append(t)
    if missing:
        print(f'\n  [Charges] No fixed charge found for types: {set(missing)} — set to 0')
    return charges


def charges_am1bcc(atoms, bonds, mol2_path, **_):
    """
    Compute AM1-BCC charges via antechamber (AmberTools).
    Requires antechamber to be installed and in PATH.
    """
    import tempfile, shutil

    antechamber = shutil.which('antechamber')
    if not antechamber:
        print('\n  [Charges] antechamber not found in PATH.')
        print('           Install AmberTools: conda install -c conda-forge ambertools')
        print('           Falling back to Gasteiger charges.')
        return charges_gasteiger(atoms, bonds)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_mol2 = os.path.join(tmpdir, 'out_bcc.mol2')
        cmd = [
            antechamber,
            '-i', mol2_path, '-fi', 'mol2',
            '-o', out_mol2, '-fo', 'mol2',
            '-c', 'bcc', '-s', '2', '-nc', '0',
        ]
        print('\n  [Charges] Running antechamber AM1-BCC …')
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmpdir)
        if result.returncode != 0 or not os.path.isfile(out_mol2):
            print('  [Charges] antechamber failed — falling back to Gasteiger.')
            print(f'            {result.stderr[:200]}')
            return charges_gasteiger(atoms, bonds)

        bcc_atoms, _ = parse_mol2(out_mol2)
        bcc_map = {a['id']: a['charge'] for a in bcc_atoms}
        return {a['id']: bcc_map.get(a['id'], a['charge']) for a in atoms}


#  CHARGE METHOD MENU
def build_charge_menu(frc_path, input_format='mol2'):
    """
    Return list of (label, key, description) tuples for the charge method menu,
    based on what the chosen force field supports.
    """
    family = ff_family(frc_path)
    options = []

    has_bi  = any(k in family for k in FF_HAS_BOND_INCREMENTS)
    has_fix = any(k in family for k in FF_HAS_FIXED_CHARGES)

    if has_bi:
        options.append(('Bond increments', 'bond_increments',
                        'Native PCFF/CFF91/COMPASS method — reads Δ values from .frc  [Recommended]'))
    if has_fix:
        options.append(('Fixed charges', 'fixed',
                        'Published charges embedded in .frc atom_types  [Recommended for CLAYFF]'))

    # Gasteiger charges depend on the molecular graph and element types, not
    # on the selected force field, so this option is valid for every FF.
    options.append(('Gasteiger', 'gasteiger',
                    'Gasteiger-Marsili via RDKit/AmberTools/Open Babel  [Any FF]'))

    if input_format == 'mol2':
        antechamber_avail = bool(subprocess.run(['which','antechamber'],
                                 capture_output=True).stdout)
        suffix = '[~seconds]' if antechamber_avail else '[needs: conda install ambertools]'
        options.append(('AM1-BCC', 'am1bcc',
                        f'Semi-empirical QM via antechamber — best for GAFF/AMBER  {suffix}'))
        options.append(('MOL2 charges', 'input', 'Keep the partial charges stored in the MOL2 file'))
    elif input_format == 'lammps':
        options.append(('Input charges', 'input', 'Keep charges from the LAMMPS Atoms section'))
    else:
        options.append(('Zero charges', 'input', 'POSCAR contains no partial charges'))
    return options


def select_charge_method(frc_path, input_format='mol2'):
    """Present charge method menu and return chosen key."""
    options = build_charge_menu(frc_path, input_format)
    family  = ff_family(frc_path)

    print()
    print('╔══════════════════════════════════════════════════════════════╗')
    print('║              Partial Charge Method                          ║')
    print('╠══════════════════════════════════════════════════════════════╣')
    for i, (label, key, desc) in enumerate(options, 1):
        print(f'║  {i}.  {label:<20}                              ║')
        print(f'║      {desc[:56]:<56}  ║')
        if i < len(options):
            print('║                                                              ║')
    print('╚══════════════════════════════════════════════════════════════╝')
    print()

    while True:
        try:
            raw = input(f'  Select charge method [1-{len(options)}]: ').strip()
            idx = int(raw)
            if 1 <= idx <= len(options):
                label, key, _ = options[idx - 1]
                return key, label
            print(f'  Please enter a number between 1 and {len(options)}.')
        except (ValueError, EOFError):
            print('  Invalid input.')


def compute_charges(method_key, atoms, bonds, type_assignments,
                    ff_types, frc_path, mol2_path):
    """Dispatch to the correct charge function."""
    kwargs = dict(
        atoms=atoms, bonds=bonds, type_assignments=type_assignments,
        ff_types=ff_types, frc_path=frc_path, mol2_path=mol2_path,
    )
    if method_key == 'bond_increments':
        return charges_bond_increments(**kwargs)
    elif method_key == 'fixed':
        return charges_fixed(**kwargs)
    elif method_key == 'am1bcc':
        return charges_am1bcc(**kwargs)
    elif method_key == 'input':
        return charges_input(**kwargs)
    else:  # gasteiger
        return charges_gasteiger(**kwargs)

#  Topology

def avg_bond_angle(atom, atoms_by_id, neighbours):
    aid  = atom['id']
    nbs  = neighbours.get(aid, [])
    if len(nbs) < 2:
        return 0.0
    p1 = (atom['x'], atom['y'], atom['z'])
    angles = []
    for a, b in itertools.combinations(nbs, 2):
        p0 = (atoms_by_id[a]['x'], atoms_by_id[a]['y'], atoms_by_id[a]['z'])
        p2 = (atoms_by_id[b]['x'], atoms_by_id[b]['y'], atoms_by_id[b]['z'])
        v1 = tuple(p0[i]-p1[i] for i in range(3))
        v2 = tuple(p2[i]-p1[i] for i in range(3))
        dot = sum(v1[i]*v2[i] for i in range(3))
        m1  = math.sqrt(sum(x**2 for x in v1))
        m2  = math.sqrt(sum(x**2 for x in v2))
        if m1 > 0 and m2 > 0:
            angles.append(math.degrees(math.acos(max(-1.0, min(1.0, dot/(m1*m2))))))
    return sum(angles)/len(angles) if angles else 0.0


#  OUTPUTS

def write_nta(atoms, type_assignments, ff_name, charge_method, out_path):
    stem = os.path.splitext(os.path.basename(out_path))[0]
    with open(out_path, 'w') as fh:
        fh.write(f'New type assignment file for {stem}.data '
                 f'> atom_typing v3.0 / {ff_name} / charges: {charge_method}\n')
        fh.write('style id\n')
        for atom in atoms:
            aid  = atom['id']
            atyp = type_assignments[aid]
            if ff_name == HPAN_HC_FF_NAME and atyp == 'hpan':
                atyp = 'hc'
            fh.write(f'{aid:<6} {atyp:<15} #  Correctly found\n')
    print(f'  [NTA]  {out_path}')


def write_data(atoms, bonds, type_assignments, charges,
               ff_name, charge_method, atoms_by_id, neighbours, out_path,
               metadata=None):
    """Write a LAMMPS  .data file."""
    mol_stem = os.path.splitext(os.path.basename(out_path))[0]

    # Element type numbering (stable order)
    elements_seen = []
    for a in atoms:
        el = get_element(a['sybyl'])
        if el not in elements_seen:
            elements_seen.append(el)
    elements_seen.sort(key=lambda e: ELEMENT_ORDER.index(e)
                       if e in ELEMENT_ORDER else 99)
    el_to_type = {el: i+1 for i, el in enumerate(elements_seen)}

    metadata = metadata or {}
    # Centroid-centre coordinates, except when an existing LAMMPS box is known.
    n  = len(atoms)
    preserve_box = all(k in metadata.get('bounds', {}) for k in 'xyz')
    cx = 0.0 if preserve_box else sum(a['x'] for a in atoms) / n
    cy = 0.0 if preserve_box else sum(a['y'] for a in atoms) / n
    cz = 0.0 if preserve_box else sum(a['z'] for a in atoms) / n
    for a in atoms:
        a['xc'] = a['x'] - cx
        a['yc'] = a['y'] - cy
        a['zc'] = a['z'] - cz

    if preserve_box:
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = (metadata['bounds'][k] for k in 'xyz')
    else:
        PAD = 0.5
        xlo, xhi = min(a['xc'] for a in atoms)-PAD, max(a['xc'] for a in atoms)+PAD
        ylo, yhi = min(a['yc'] for a in atoms)-PAD, max(a['yc'] for a in atoms)+PAD
        zlo, zhi = min(a['zc'] for a in atoms)-PAD, max(a['zc'] for a in atoms)+PAD

    lines = []

    # Header
    lines.append(
        f'HEADER,  {os.path.basename(out_path).replace(".data",".mol2")}  '
        f'read w/ atom_typing v3.0 > {ff_name} / charges: {charge_method}\n'
    )
    lines.append('\n')
    lines.append(f'{len(atoms)} atoms\n')
    lines.append(f'{len(bonds)} bonds\n')
    lines.append('\n')
    lines.append(f'{len(elements_seen)} atom types\n')
    lines.append(f'1 bond types\n')
    lines.append('\n')
    lines.append(f' {xlo:.10f}   {xhi:.10f}    xlo  xhi  \n')
    lines.append(f' {ylo:.10f}   {yhi:.10f}    ylo  yhi  \n')
    lines.append(f' {zlo:.10f}   {zhi:.10f}    zlo  zhi  \n')
    if preserve_box and 'tilts' in metadata['bounds']:
        xy, xz, yz = metadata['bounds']['tilts']
        lines.append(f' {xy:.10f}   {xz:.10f}   {yz:.10f}    xy xz yz\n')
    lines.append('\n')

    # Masses
    lines.append('Masses\n\n')
    for el, tid in el_to_type.items():
        mass = ELEMENT_MASS.get(el, 1.0)
        lines.append(f' {tid}   {mass:<9}  #  {el:<5}\n')
    lines.append('\n')

    # Bond Coeffs
    lines.append('Bond Coeffs\n\n')
    lines.append(' 1  #      1\n')
    lines.append('\n')

    # Atoms # full
    lines.append('Atoms # full\n\n')
    mol_id = 1
    for a in atoms:
        aid    = a['id']
        el     = get_element(a['sybyl'])
        atype  = el_to_type[el]
        q      = charges.get(aid, 0.0)
        x, y, z = a['xc'], a['yc'], a['zc']
        pcff   = type_assignments[aid]
        hybr   = SYBYL_HYB_LABEL.get(a['sybyl'], 'Sp3')
        avg_a  = avg_bond_angle(a, atoms_by_id, neighbours)
        lines.append(
            f'  {aid:<5}  {mol_id}   {atype}  {q:+.6f}   '
            f'{x:.9f}    {y:.9f}     {z:.9f}   0  0  0'
            f'    #     {el:<5}  {hybr:>8}     avg-angle: {avg_a:>11.6f}'
            f'      type: {pcff}\n'
        )
    lines.append('\n')

    # Bonds
    lines.append('Bonds\n\n')
    for b in bonds:
        lines.append(f'{b["id"]:<2} 1  {b["a1"]:<2} {b["a2"]:<2}\n')

    with open(out_path, 'w') as fh:
        fh.writelines(lines)
    print(f'  [DATA] {out_path}')

#  INTERACTIVE MENUS

def find_frc_files(search_dirs):
    found = []
    for d in search_dirs:
        found.extend(glob.glob(os.path.join(d, '*.frc')))
    seen, unique = set(), []
    for p in found:
        canon = os.path.realpath(p)
        if canon not in seen:
            seen.add(canon); unique.append(p)
    return sorted(unique, key=lambda p: os.path.basename(p).lower())


def select_ff_interactive(frc_files):
    print()
    print('╔══════════════════════════════════════════════════════════════╗')
    print('║              Available Force Fields                         ║')
    print('╠══════════════════════════════════════════════════════════════╣')
    for i, path in enumerate(frc_files, 1):
        stem = os.path.splitext(os.path.basename(path))[0]
        desc = describe_ff(stem)
        print(f'║  {i:>2}.  {stem:<30}  ║')
        print(f'║        {desc:<52}  ║')
        if i < len(frc_files):
            print('║                                                              ║')
    print('╚══════════════════════════════════════════════════════════════╝')
    print()
    while True:
        try:
            raw = input(f'  Select force field [1-{len(frc_files)}]: ').strip()
            idx = int(raw)
            if 1 <= idx <= len(frc_files):
                return frc_files[idx - 1]
            print(f'  Please enter a number between 1 and {len(frc_files)}.')
        except (ValueError, EOFError):
            print('  Invalid input.')


def select_input_interactive():
    print()
    while True:
        raw = input('  Path to MOL2, LAMMPS data, or POSCAR file: ').strip()
        if not raw: continue
        if not os.path.isfile(raw):
            print(f'  File not found: {raw}'); continue
        return raw


def select_frc_path_interactive():
    print()
    while True:
        raw = input('  Path to force-field (.frc) file: ').strip()
        if not raw: continue
        if not os.path.isfile(raw):
            print(f'  File not found: {raw}'); continue
        return raw


def select_output_name_interactive():
    print()
    while True:
        raw = input('  Output name: ').strip()
        if not raw: continue
        return os.path.splitext(raw)[0]

def main():
    args = sys.argv[1:]
    input_path = frc_path = ffdir = None
    i = 0
    while i < len(args):
        if args[i] == '--ffdir' and i+1 < len(args):
            ffdir = args[i+1]; i += 2
        elif input_path is None:
            input_path = args[i]; i += 1
        elif frc_path is None:
            frc_path = args[i]; i += 1
        else:
            i += 1

    print()
    print('══════════════════════════════════════════════════════════')
    print('       ATOM TYPING + CHARGE + DATA FILE TOOL  v3.0       ')
    print('══════════════════════════════════════════════════════════')

    # ── structure input ───────────────────────────────────────────────────
    if input_path is None:
        input_path = select_input_interactive()
    if not os.path.isfile(input_path):
        print(f'\n[ERROR] Input file not found: {input_path}'); sys.exit(1)
    try:
        input_format = detect_input_format(input_path)
    except ValueError as exc:
        print(f'\n[ERROR] {exc}'); sys.exit(1)

    # ── force field ───────────────────────────────────────────────────────
    if frc_path is None:
        frc_path = select_frc_path_interactive()
    elif not os.path.isfile(frc_path):
        print(f'\n[ERROR] Force-field file not found: {frc_path}'); sys.exit(1)

    ff_name = os.path.splitext(os.path.basename(frc_path))[0]

    # ── charge method ─────────────────────────────────────────────────────
    charge_key, charge_label = select_charge_method(frc_path, input_format)

    # ── output paths ──────────────────────────────────────────────────────
    out_name = select_output_name_interactive()
    out_dir  = os.getcwd() if os.access(os.getcwd(), os.W_OK) else \
               os.path.dirname(os.path.realpath(input_path))
    base     = os.path.join(out_dir, out_name)
    nta_path  = base + '.nta'
    data_path = base + '.data'
    score_path = base + '_typing_scores.txt'

    # ── load FF ───────────────────────────────────────────────────────────
    print(f'\n[FF]   {ff_name}')
    print(f'       {describe_ff(ff_name)}')
    ff_types = parse_frc(frc_path)
    precompute_ff_features(ff_types)
    print(f'       {len(ff_types)} atom types loaded')

    # ── load structure ────────────────────────────────────────────────────
    print(f'\n[{input_format}] {os.path.basename(input_path)}')
    try:
        atoms, bonds, metadata, input_format = parse_structure(input_path)
    except (ValueError, OSError) as exc:
        print(f'\n[ERROR] Could not parse input: {exc}'); sys.exit(1)
    neighbours, btm = build_graph(bonds)
    ring_sizes      = atom_ring_sizes(neighbours)
    atoms_by_id     = {a['id']: a for a in atoms}
    print(f'       {len(atoms)} atoms, {len(bonds)} bonds')

    # ── assign types ──────────────────────────────────────────────────────
    print('\n  Assigning atom types …')
    type_assignments = assign_types(
        atoms, atoms_by_id, neighbours, btm, ring_sizes, ff_types)
    apply_pcff_iff_hpan_override(
        frc_path, atoms, atoms_by_id, neighbours, type_assignments)
    changed = optimize_types_for_parameter_coverage(
        atoms, bonds, atoms_by_id, neighbours, btm, ring_sizes, ff_types,
        type_assignments, frc_path, report_path=score_path)
    print(f'  Parameter-coverage optimization changed {changed} atom assignment(s).')
    print(f'  Typing score report: {score_path}')

    # ── compute charges ───────────────────────────────────────────────────
    print(f'  Computing charges ({charge_label}) …')
    try:
        charges = compute_charges(
            charge_key, atoms, bonds, type_assignments,
            ff_types, frc_path, input_path)
    except RuntimeError as exc:
        print(f'\n[ERROR] {exc}')
        sys.exit(1)

    total_q = sum(charges.values())

    # ── summary table ─────────────────────────────────────────────────────
    print()
    print(f"  {'ID':>4}  {'Name':<5}  {'Sybyl':<7}  {'FF type':<8}  {'Charge':>8}")
    print(f"  {'─'*4}  {'─'*5}  {'─'*7}  {'─'*8}  {'─'*8}")
    for atom in atoms:
        aid = atom['id']
        q   = charges.get(aid, 0.0)
        print(f"  {aid:>4}  {atom['name']:<5}  {atom['sybyl']:<7}  "
              f"{type_assignments[aid]:<8}  {q:>+8.4f}")
    print(f"\n  Total charge: {total_q:+.6f} e")

    # ── write outputs ──────────────────────────────────────────────────────
    print()
    write_nta(atoms, type_assignments, ff_name, charge_label, nta_path)
    write_data(atoms, bonds, type_assignments, charges,
               ff_name, charge_label, atoms_by_id, neighbours, data_path, metadata)

    print(f'\n  Done.\n')


if __name__ == '__main__':
    main()

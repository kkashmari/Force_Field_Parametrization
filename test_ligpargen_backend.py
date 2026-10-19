from __future__ import annotations

import tempfile
import unittest
import zipfile
from pathlib import Path

from ligpargen_backend import (
    LigParGenError, import_ligpargen, repair_known_ligpargen_header_issues,
    validate_ligpargen_text,
)


SAMPLE = """LAMMPS data file generated with LigParGen 2.1 - Jorgensen Lab

2 atoms
1 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
1 bond types
0 angle types
0 dihedral types
0 improper types

-5 5 xlo xhi
-5 5 ylo yhi
-5 5 zlo zhi

Masses

1 12.011
2 1.008

Pair Coeffs

1 0.066 3.500
2 0.030 2.500

Bond Coeffs

1 340.0 1.09

Angle Coeffs

Dihedral Coeffs

Improper Coeffs

Atoms

1 1 1 -0.1 0 0 0
2 1 2  0.1 1 0 0

Bonds

1 1 1 2

Angles

Dihedrals

Impropers
"""


class LigParGenBackendTests(unittest.TestCase):
    def test_validation(self):
        audit = validate_ligpargen_text(SAMPLE)
        self.assertEqual(audit["counts"]["atoms"], 2)
        self.assertAlmostEqual(audit["net_charge"], 0.0)

    def test_rejects_out_of_range_type(self):
        broken = SAMPLE.replace("2 1 2  0.1", "2 1 3  0.1")
        with self.assertRaises(LigParGenError):
            validate_ligpargen_text(broken)

    def test_repairs_only_legacy_improper_count(self):
        broken = SAMPLE.replace("0 impropers", "3 impropers", 1)
        repaired, notes = repair_known_ligpargen_header_issues(broken)
        self.assertIn("0 impropers", repaired)
        self.assertEqual(len(notes), 1)
        validate_ligpargen_text(repaired)

    def test_zip_and_verbatim_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "download.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr("MOL.lmp", SAMPLE)
            out = root / "out"
            audit = import_ligpargen(archive, out, "MOL")
            copied = Path(audit["data_file"])
            self.assertEqual(copied.read_text(), SAMPLE)
            self.assertFalse(audit["tinker_parameters_used"])
            self.assertIn(
                "improper_style  cvff", Path(audit["lammps_input"]).read_text())


if __name__ == "__main__":
    unittest.main()

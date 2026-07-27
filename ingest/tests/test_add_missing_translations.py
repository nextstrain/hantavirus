"""Tests for ingest/scripts/add_missing_translations.py.

Runs in any Nextstrain env (unittest is stdlib):

    python -m unittest discover -s ingest/tests

Also runnable directly:

    python ingest/tests/test_add_missing_translations.py
"""

import sys
import unittest
from pathlib import Path

from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from add_missing_translations import find_largest_orf  # noqa: E402


class FindLargestOrf(unittest.TestCase):
    def test_start_snaps_to_met_not_stretch_boundary(self):
        # Regression: within a stop-delimited stretch the ORF must start at the
        # first Met, and the *coordinate* must advance to that Met as well. The
        # old code advanced only the translation, leaving the start at the stretch
        # boundary, so seq[start:end] no longer matched the stored translation.
        #
        #   frame 0:  CAA  CAA  ATG  GGT  GCT  AAG  ... (ORF body) ...  TAA
        #              Q    Q    M    G    A    K                        *
        #   nt:        0    3    6    9   12   15                       135
        #                        ^ largest ORF starts at the Met -> start = 6
        #                          (the leading "QQ" at nt 0-5 is NOT part of it)
        seq = Seq(
            "CAACAAATGGGTGCTAAGCATTCTGGTAAGTCTCCAGTATTTCATTGGAAGTGGGATGTAA"
            "AGCATTATGTATGGAAGGATTATGATACTCCATCTCATGATCCAGGTCCAGTAACTTTTCAT"
            "TCTTCTCATAAGTAA"
        )
        start, end, protein = find_largest_orf(seq)

        self.assertEqual(start, 6)  # the Met, not the "QQ" boundary at 0
        self.assertTrue(protein.startswith("MGAKHSGK"))  # ORF reads M G A K H S G K ...
        self.assertEqual(len(protein), 43)
        # The coordinates reproduce the returned translation exactly (the bug):
        self.assertEqual(str(seq[start:end].translate()), protein)

    def test_largest_orf_is_not_the_first(self):
        # The largest ORF may come *after* an earlier, shorter ORF and its stop.
        #
        #   frame 0:  ATG  AAA  TAA | ATG  TTT  TGG  ... TAT  TAA
        #              M    K    *  |  M    F    W        Y    *
        #   nt:        0    3    6     9   12   15       51   54
        #             \--- ORF1 ---/  \------ largest ORF ------/
        #                                ^ start = 9 (the second Met)
        seq = Seq("ATGAAATAAATGTTTTGGAAGGCTGTAGTATATGCTACTAAGGGTACTTTTTATTAA")
        start, end, protein = find_largest_orf(seq)

        self.assertEqual((start, end, protein), (9, 54, "MFWKAVVYATKGTFY"))
        self.assertEqual(str(seq[start:end].translate()), protein)

    def test_stop_before_end_excludes_trailing_sequence(self):
        # The ORF ends at its stop codon; the stop itself and any sequence after
        # it are excluded from the coordinates (end < len(seq)).
        #
        #   frame 0:  ATG ... GTA ACT | TAA | GGG TAA CCC TAG TTT
        #              M  ...  V   T  |  *  |  (trailing, excluded)
        #   nt:        0          48    51    54                  69
        #             \--- ORF (17 aa) --/     ^ end = 51, not 69
        seq = Seq(
            "ATGTTTTTTACTAAGGATGATAAGGTACCACATTTTCCACATGATGTAACTTAAGGGTAACCCTAGTTT"
        )
        start, end, protein = find_largest_orf(seq)

        self.assertEqual((start, end), (0, 51))
        self.assertEqual(protein, "MFFTKDDKVPHFPHDVT")
        self.assertLess(end, len(seq))  # trailing sequence excluded
        self.assertEqual(str(seq[start:end].translate()), protein)

    def test_no_met_uses_whole_stretch(self):
        # A stop-delimited stretch with no Met is used in full.
        #   frame 0:  AAA  AAA  TAA  ->  K K *  -> "KK" over nt 0..6
        self.assertEqual(find_largest_orf(Seq("AAAAAATAA")), (0, 6, "KK"))

    def test_no_codons_returns_empty(self):
        self.assertEqual(find_largest_orf(Seq("AT")), (0, 0, ""))


if __name__ == "__main__":
    unittest.main()

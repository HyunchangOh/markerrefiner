"""Behaviour the filter depends on."""

import os
import tempfile
import unittest

from mmprefilter.qc import (
    annotate_seq_id,
    gene_id,
    longest_isoform_records,
    n50,
    parse_busco_short_summary,
    parse_species_from_filename,
)
from mmprefilter.structure import (
    genomic_intron_lengths,
    msa_segments_from_ref,
    row_matches,
    strip_intron_n,
    summarise_lengths,
)


THRESHOLDS = {
    "n_exons_gt": 0,
    "n_introns_lt": 5,
    "max_exon_bp_gt": 200,
    "max_intron_bp_lt": 100,
    "avg_exon_bp_gt": 100,
    "avg_intron_bp_lt": 50,
}


class StructureTests(unittest.TestCase):
    def test_strip_intron_mask(self):
        self.assertEqual(strip_intron_n("ACGTNNNAC"), "ACGTAC")

    def test_segments_skip_gaps(self):
        segments = msa_segments_from_ref("AC--GTN-NA")
        kinds = [kind for kind, _ in segments]
        self.assertEqual(kinds, ["coding", "intron", "coding"])
        coding = [columns for kind, columns in segments if kind == "coding"]
        # Gaps are removed before ACGT/N runs are split, so AC and GT are one exon.
        self.assertEqual(coding[0], [0, 1, 4, 5])
        self.assertEqual(coding[1], [9])

    def test_blank_metric_does_not_fail(self):
        row = {
            "n_exons": "2",
            "n_introns": "1",
            "max_exon_bp": "300",
            "max_intron_bp": "",
            "avg_exon_bp": "250",
            "avg_intron_bp": "",
        }
        self.assertTrue(row_matches(row, THRESHOLDS))

    def test_intron_cutoff_fails(self):
        row = {
            "n_exons": "2",
            "n_introns": "8",
            "max_exon_bp": "300",
            "max_intron_bp": "10",
            "avg_exon_bp": "250",
            "avg_intron_bp": "10",
        }
        self.assertFalse(row_matches(row, THRESHOLDS))

    def test_genomic_intron_on_plus_strand(self):
        blocks = [
            {"chr": "Chr1", "start": 10, "end": 20, "strand": "+"},
            {"chr": "Chr1", "start": 31, "end": 40, "strand": "+"},
        ]
        self.assertEqual(genomic_intron_lengths(blocks), [10])

    def test_summary_leaves_missing_intron_blank(self):
        summary = summarise_lengths([100, 300], [])
        self.assertEqual(summary["n_exons"], 2)
        self.assertEqual(summary["max_exon_bp"], 300)
        self.assertEqual(summary["max_intron_bp"], "")


class QcTests(unittest.TestCase):
    def test_trinity_gene_and_longest_isoform(self):
        records = [
            ("TRINITY_DN1_c0_g1_i1", "A" * 10),
            ("TRINITY_DN1_c0_g1_i2", "A" * 25),
            ("TRINITY_DN2_c0_g1_i1", "C" * 4),
        ]
        self.assertEqual(gene_id(records[0][0]), "TRINITY_DN1_c0_g1")
        kept = longest_isoform_records(records)
        self.assertEqual([seq_id for seq_id, _ in kept], [
            "TRINITY_DN1_c0_g1_i2",
            "TRINITY_DN2_c0_g1_i1",
        ])

    def test_n50(self):
        self.assertEqual(n50([10, 10, 10, 50]), 50)

    def test_species_name_and_header_typo(self):
        code, species = parse_species_from_filename(
            "DAT5-echitum_plantagineum_localscratch_trinity.Trinity.fasta"
        )
        self.assertEqual((code, species), ("DAT5", "echitum_plantagineum"))
        renamed = annotate_seq_id(
            "DAT5_TRINITY_DN1_c0_g1_i1",
            {"DAT5": species},
            [["echitum_plantagineum", "echium_plantagineum"]],
        )
        self.assertEqual(renamed, "DAT5_echium_plantagineum_TRINITY_DN1_c0_g1_i1")

    def test_busco_short_summary(self):
        text = (
            "# Summarized benchmarking in BUSCO notation for file "
            "echium_vulgare.Trinity.fasta\n"
            "\tC:41.4%[S:25.9%,D:15.5%],F:15.9%,M:42.7%,n:2805\n"
            "\t1160\tComplete BUSCOs (C)\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "short_summary.txt")
            with open(path, "w") as handle:
                handle.write(text)
            parsed = parse_busco_short_summary(path, tmp)
        self.assertEqual(parsed["sample"], "echium_vulgare")
        self.assertEqual(parsed["Cpct"], "41.4")
        self.assertEqual(parsed["C"], "1160")
        self.assertEqual(parsed["n"], "2805")


if __name__ == "__main__":
    unittest.main()

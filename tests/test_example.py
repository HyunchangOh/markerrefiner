"""The shipped example must reproduce the documented keep and drop decisions."""

import os
import tempfile
import unittest

from mmprefilter.config import load_config
from mmprefilter.pipeline import run_config
from mmprefilter.structure import read_table


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class ExamplePipelineTests(unittest.TestCase):
    def test_example_config(self):
        cfg = load_config(os.path.join(REPO, "example", "config.json"))
        with tempfile.TemporaryDirectory() as tmp:
            cfg["output_dir"] = tmp
            run_config(cfg)
            cds = read_table(os.path.join(tmp, "reference_cds_check.tsv"))
            self.assertEqual([row["status"] for row in cds], ["MATCH", "MATCH", "MATCH"])
            with open(os.path.join(tmp, "filtered", "athaliana", "namelist.txt")) as handle:
                kept = handle.read().split()
            self.assertEqual(kept, ["AT1G01180", "AT1G01280"])
            with open(os.path.join(tmp, "filtered", "vitis", "namelist.txt")) as handle:
                vitis = handle.read().strip()
            self.assertEqual(vitis, "")
            alignment = os.path.join(
                tmp,
                "filtered",
                "athaliana",
                "AT1G01280.mafft.added.reference.align.fna",
            )
            with open(alignment) as handle:
                text = handle.read()
            self.assertIn(">DAT5_echium_plantagineum_TRINITY_DN15068_c0_g1_i3", text)
            self.assertNotIn("echitum", text)
            report = os.path.join(tmp, "report.html")
            with open(report) as handle:
                page = handle.read()
            self.assertIn("Download thresholds.json", page)
            self.assertIn("AT1G01180", page)
            self.assertIn("VV00G01720", page)
            comparison = read_table(os.path.join(tmp, "gff_vs_msa.tsv"))
            self.assertTrue(all(row["coding_bp_match"] == "yes" for row in comparison))


if __name__ == "__main__":
    unittest.main()

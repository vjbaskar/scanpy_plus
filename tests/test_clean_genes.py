from pathlib import Path
import unittest


class CleanGenesTest(unittest.TestCase):
    def test_download_cc_uses_bundled_resource(self):
        gene_file = (
            Path(__file__).resolve().parents[1]
            / "scanpy_plus"
            / "data"
            / "regev_lab_cell_cycle_genes.txt"
        )
        self.assertTrue(gene_file.is_file())
        genes = gene_file.read_text().splitlines()
        self.assertIn("MCM5", genes)
        self.assertIn("CENPA", genes)


if __name__ == "__main__":
    unittest.main()

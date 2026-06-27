import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("probeat_cli", ROOT / "workflow" / "probeat_cli.py")
CLI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLI)


class CliTests(unittest.TestCase):
    def test_profiles_are_defined(self):
        for name in ("core", "full", "probiotic", "ont"):
            self.assertIn("tools", CLI.profile_config(name))
        self.assertIn("fastani", CLI.profile_config("ont")["tools"])

    def test_validate_sample_rejects_missing_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            read = Path(directory) / "r1.fastq"
            read.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
            with self.assertRaises(SystemExit):
                CLI.validate_sample("sample", "PE", [str(read)])

    def test_validate_sample_accepts_ont(self):
        with tempfile.TemporaryDirectory() as directory:
            read = Path(directory) / "reads.fastq"
            read.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
            self.assertEqual([read.resolve()], CLI.validate_sample("sample_1", "ONT", [str(read)]))

    def test_manifests_have_no_local_workspace_path(self):
        for path in (ROOT / "manifests").glob("*.yaml"):
            self.assertNotIn("/home/vebaev/codex", path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()

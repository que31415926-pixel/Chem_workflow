"""toolchain.lock.json, the BYOL Dockerfile, environment.yml and cloud_setup must agree."""
import json
import pathlib
import re
import unittest

from chemflow2 import cloud_setup

V2 = pathlib.Path(__file__).resolve().parents[1]
REPO = V2.parent
LOCK = json.loads((REPO / "toolchain.lock.json").read_text())


class ToolchainLock(unittest.TestCase):
    def dockerfiles(self):
        paths = [V2 / "deploy" / "byol.Dockerfile", V2 / "chemflow2" / "_deploy" / "byol.Dockerfile"]
        return [p for p in paths if p.exists()]

    def test_dockerfile_conda_pins_match_lock(self):
        self.assertTrue(self.dockerfiles())
        for path in self.dockerfiles():
            text = path.read_text()
            line = next(l for l in text.splitlines() if "micromamba install" in l)
            pins = dict(re.findall(r"(?<![\w-])([a-z][\w.-]*)=([\w.]+)", line.replace('"', "")))
            self.assertEqual(pins, LOCK["conda_forge"], path)

    def test_dockerfile_images_and_openmpi_match_lock(self):
        for path in self.dockerfiles():
            text = path.read_text()
            self.assertIn(f"FROM {LOCK['images']['final_stage']}", text, path)
            self.assertIn(f"FROM {LOCK['images']['orca_and_openmpi_build_stages']}", text, path)
            self.assertIn(f"ARG OMPI_VERSION={LOCK['openmpi']['version']}", text, path)
            self.assertIn("OMPI_MCA_plm=isolated", text, path)

    def test_environment_yml_matches_lock(self):
        text = (REPO / "environment.yml").read_text()
        for name in ("python", "crest", "xtb", "tblite", "rdkit", "numpy", "scipy"):
            self.assertIn(f"  - {name}={LOCK['conda_forge'][name]}\n", text, name)
        self.assertIn(f"python-libsbml=={LOCK['conda_forge']['python-libsbml']}", text)

    def test_cloud_setup_accepts_only_locked_installer(self):
        locked = LOCK["orca"]["platforms"]["linux-x86_64"]["installer"]
        self.assertTrue(cloud_setup.INSTALLER_NAME.fullmatch(locked))
        for other in ("orca_6_1_1_linux_x86-64_shared_openmpi418.run",
                      "orca_6_1_0_linux_x86-64_shared_openmpi416.run",
                      "orca_6_0_1_linux_x86-64_shared_openmpi416.run"):
            self.assertIsNone(cloud_setup.INSTALLER_NAME.fullmatch(other), other)

    def test_orca_version_consistent(self):
        self.assertEqual(LOCK["orca"]["version"], "6.1.0")
        name = LOCK["orca"]["platforms"]["linux-x86_64"]["installer"]
        self.assertTrue(name.startswith("orca_" + LOCK["orca"]["version"].replace(".", "_") + "_"))
        self.assertIn(LOCK["openmpi"]["version"].replace(".", ""), name)


if __name__ == "__main__":
    unittest.main()

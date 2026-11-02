"""Meaningful offline tests for the new distribution interface."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('prepare_data', ROOT/'tools/prepare_data.py')
prepare=importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)

class ReleaseTests(unittest.TestCase):
    def test_archive_traversal_and_absolute_paths_rejected(self):
        for name in ['../escape.txt','/escape.txt','C:/escape.txt','nested/../../escape.txt',r'..\escape.txt',r'nested\..\..\escape.txt']:
            with self.subTest(name=name), self.assertRaises(ValueError):
                prepare.safe_member(zipfile.ZipInfo(name))

    def test_normal_archive_member_allowed(self):
        self.assertEqual(prepare.safe_member(zipfile.ZipInfo('tables/S01_genes.csv')), 'tables/S01_genes.csv')

    def test_launcher_lists_real_stages_without_optional_packages(self):
        p=subprocess.run([sys.executable,str(ROOT/'run.py'),'list'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn('nested-original',p.stdout)
        self.assertIn('virtual-ko',p.stdout)

if __name__=='__main__':
    unittest.main()

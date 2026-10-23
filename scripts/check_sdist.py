"""Verify test inputs and run offline tests from an extracted source archive."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
OFFLINE_TESTS = (
    "tests/test_sequence_data.py",
    "tests/test_ucsc_gtf.py",
    "tests/test_assembly_selection.py",
    "tests/test_genetic_code_metadata.py",
    "tests/test_annotation_genetic_code.py",
    "tests/test_legacy_cache_compatibility.py",
    "tests/test_release_cache_combinations.py",
    "tests/test_release_catalogue.py",
    "tests/test_release_comparisons.py",
    "tests/test_annotation_integrity.py",
    "tests/test_cli_errors.py",
)


def check_sdist(archive_path):
    with TemporaryDirectory(prefix="pyensembl-sdist-") as temporary:
        destination = Path(temporary)
        with tarfile.open(archive_path) as archive:
            archive.extractall(destination, filter="data")
        roots = list(destination.iterdir())
        if len(roots) != 1 or not roots[0].is_dir():
            raise ValueError("Expected one source-distribution root directory")
        extracted = roots[0]
        # A checkout can independently prove that the manifest included every
        # tracked test input and omitted developer-local/generated test files.
        # Extracted sdists do not need Git to run their offline smoke tests.
        if (ROOT / ".git").exists():
            tracked = subprocess.check_output(
                ["git", "ls-files", "-z", "tests/"], cwd=ROOT,
            ).decode().split("\0")
            expected = {name for name in tracked if name}
            actual = {path.relative_to(extracted).as_posix()
                      for path in (extracted / "tests").rglob("*") if path.is_file()}
            if expected != actual:
                raise ValueError("Source archive test inputs differ: missing=%s, extra=%s"
                                 % (sorted(expected - actual), sorted(actual - expected)))
            for name in expected:
                if (ROOT / name).read_bytes() != (extracted / name).read_bytes():
                    raise ValueError("Source archive changed test input: " + name)
            print("Verified %d tracked test files and fixtures" % len(expected), flush=True)
        environment = dict(os.environ)
        for name in ("PYTHONPATH", "PYENSEMBL_NETWORK_TESTS", "PYENSEMBL_BENCHMARKS"):
            environment.pop(name, None)
        environment["PYENSEMBL_CACHE_DIR"] = str(destination / "cache")
        environment["XDG_CACHE_HOME"] = str(destination / "xdg-cache")
        # -I excludes the invoking checkout and PYTHONPATH. Explicitly import
        # the extracted package, so an editable install cannot mask omissions.
        runner = """
import pathlib
import sys
sys.path.insert(0, sys.argv[1])
import pyensembl
assert pathlib.Path(pyensembl.__file__).resolve().is_relative_to(pathlib.Path(sys.argv[1]).resolve()), pyensembl.__file__
import pytest
raise SystemExit(pytest.main(['-q', '--import-mode=importlib', *sys.argv[2:]]))
"""
        subprocess.run(
            [sys.executable, "-I", "-c", runner, str(extracted), *OFFLINE_TESTS],
            cwd=extracted, env=environment, check=True,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    check_sdist(parser.parse_args().archive.resolve())

"""Check that downstream callers get useful types and argument diagnostics."""

from pathlib import Path
import subprocess
import sys


def test_typed_consumer_and_invalid_arguments(tmp_path):
    root = Path(__file__).resolve().parents[1]
    consumer = tmp_path / "consumer.py"
    consumer.write_text('''from typing import Any, Literal, assert_type
from pyensembl.database import Database
from pyensembl.normalization import normalize_chromosome, normalize_strand

assert_type(normalize_chromosome(1), str)
assert_type(normalize_strand("+"), Literal["+", "-"])
db = Database("example.gtf")
assert_type(db.column_values_at_locus("gene_id", "gene", "1", 10), list[Any])
assert_type(db.distinct_column_values_at_locus("start", "gene", 1, 10), list[Any])
''')
    command = [sys.executable, "-m", "mypy", "--config-file", str(root / "pyproject.toml"),
               "--cache-dir", str(tmp_path / "mypy-cache"), str(consumer)]
    checked = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=60)
    assert checked.returncode == 0, checked.stdout + checked.stderr

    with consumer.open("a") as handle:
        handle.write('normalize_chromosome(1.5)\n'
                     'db.column_values_at_locus("gene_id", "gene", "1", "ten")\n')
    rejected = subprocess.run(command, cwd=root, text=True, capture_output=True, timeout=60)
    assert rejected.returncode == 1, rejected.stdout + rejected.stderr
    assert rejected.stdout.count("[arg-type]") == 2, rejected.stdout

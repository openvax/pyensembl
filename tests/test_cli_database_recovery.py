"""Explicit CLI recovery uses datacache without migrating existing caches."""

import os
from pathlib import Path
import sqlite3
import stat

import pytest

from pyensembl import Genome
from .test_cli_errors import run_cli
from .test_legacy_cache_compatibility import snapshot


@pytest.fixture
def installed(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache))
    source = tmp_path / "genes.gtf"
    source.write_text('1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g"; gene_name "G";\n')
    with Genome(reference_name="test", annotation_name="custom", gtf_path_or_url=str(source)) as genome:
        genome.index()
        database = Path(genome.db.local_db_path)
    arguments = ["install", "--reference-name", "test", "--annotation-name", "custom", "--gtf", str(source)]
    return cache, source, database, arguments


@pytest.mark.parametrize("damage", ["non-sqlite", "truncated"])
def test_explicit_overwrite_recovers_and_preserves_database_identity(installed, tmp_path, damage):
    cache, source, database, arguments = installed
    original_source = source.read_bytes()
    original_source_stat = source.stat()
    database.write_bytes(b"not a sqlite database" if damage == "non-sqlite" else database.read_bytes()[:200])
    database.chmod(0o640)
    alias = tmp_path / "existing-link.db"
    os.link(database, alias)
    original_database = database.stat()
    before = snapshot(cache)

    ordinary = run_cli(arguments, cache)
    assert ordinary.returncode == 1
    assert "DatabaseError" in ordinary.stderr
    assert snapshot(cache) == before

    recovered = run_cli([*arguments, "--overwrite"], cache)
    assert recovered.returncode == 0, recovered.stderr
    current = database.stat()
    assert (current.st_dev, current.st_ino) == (original_database.st_dev, original_database.st_ino)
    assert stat.S_IMODE(current.st_mode) == stat.S_IMODE(original_database.st_mode)
    assert (current.st_uid, current.st_gid, current.st_nlink) == (
        original_database.st_uid, original_database.st_gid, original_database.st_nlink,
    )
    with sqlite3.connect(alias) as connection:
        assert connection.execute("PRAGMA integrity_check").fetchall() == [("ok",)]
        assert connection.execute("SELECT gene_id FROM gene").fetchall() == [("g",)]
    with Genome(reference_name="test", annotation_name="custom", gtf_path_or_url=str(source)) as genome:
        assert genome.gene_by_id("g").name == "G"
    assert source.read_bytes() == original_source
    current_source_stat = source.stat()
    for attribute in ("st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_size", "st_mtime_ns", "st_ctime_ns"):
        assert getattr(current_source_stat, attribute) == getattr(original_source_stat, attribute)

    recovered_snapshot = snapshot(cache)
    reused = run_cli(arguments, cache)
    assert reused.returncode == 0, reused.stderr
    assert snapshot(cache) == recovered_snapshot


def test_invalid_replacement_gtf_preserves_corrupt_database(installed):
    cache, source, database, arguments = installed
    database.write_bytes(b"corrupt saved database")
    source.write_text("not\ta\tvalid\tGTF\n")
    before = snapshot(cache)
    original_source = source.read_bytes()
    result = run_cli([*arguments, "--overwrite"], cache)
    assert result.returncode == 1
    assert "ParsingError" in result.stderr
    assert snapshot(cache) == before
    assert source.read_bytes() == original_source

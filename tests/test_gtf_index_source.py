"""SQLite annotation indexes belong to the GTF that produced them."""

import os
from pathlib import Path
import sqlite3

import pytest

from pyensembl import Genome
from pyensembl.database import DATABASE_SCHEMA_VERSION


def write_gtf(path, name):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "g1"; gene_name "%s";\n' % name
    )


def genome_for(source, cache, **kwargs):
    return Genome(
        "test", "source-identity", gtf_path_or_url=str(source),
        cache_directory_path=str(cache), **kwargs,
    )


def test_changed_gtf_rebuilds_without_explicit_overwrite(tmp_path):
    source = tmp_path / "annotation.gtf"
    cache = tmp_path / "cache"
    write_gtf(source, "CTA")
    with genome_for(source, cache) as first:
        first.index()
        assert first.gene_name_of_gene_id("g1") == "CTA"
        before = source.stat()
        write_gtf(source, "NEW")  # Same length and path.
        os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
        with genome_for(source, cache) as second:
            assert not second.installed()
            with pytest.raises(ValueError, match="GTF database needs to be created"):
                second.gene_name_of_gene_id("g1")
            second.index()
            assert second.gene_name_of_gene_id("g1") == "NEW"
        first.index()
        assert first.gene_name_of_gene_id("g1") == "NEW"


def test_same_basename_sources_keep_separate_indexes(tmp_path):
    cache = tmp_path / "cache"
    first_source = tmp_path / "a" / "annotation.gtf"
    second_source = tmp_path / "b" / "annotation.gtf"
    write_gtf(first_source, "FIRST")
    write_gtf(second_source, "SECOND")
    with genome_for(first_source, cache) as first, genome_for(second_source, cache) as second:
        first.index()
        second.index()
        assert first.db.local_db_path != second.db.local_db_path
        assert first.gene_name_of_gene_id("g1") == "FIRST"
        assert second.gene_name_of_gene_id("g1") == "SECOND"
        assert first.installed() and second.installed()


def test_changed_copied_gtf_refreshes_copy_and_index(tmp_path):
    source = tmp_path / "annotation.gtf"
    cache = tmp_path / "cache"
    write_gtf(source, "OLD")
    options = {"copy_local_files_to_cache": True}
    with genome_for(source, cache, **options) as first:
        first.index()
        copied = Path(first.gtf_path)
        write_gtf(source, "NEW")
        with genome_for(source, cache, **options) as second:
            assert not second.installed()
            second.index()
            assert second.gene_name_of_gene_id("g1") == "NEW"
            assert 'gene_name "NEW"' in copied.read_text()
        write_gtf(source, "LAST")
        first.index()
        assert first.gene_name_of_gene_id("g1") == "LAST"


def test_copied_gtf_still_works_after_original_removed(tmp_path):
    source = tmp_path / "annotation.gtf"
    cache = tmp_path / "cache"
    write_gtf(source, "SAVED")
    options = {"copy_local_files_to_cache": True}
    with genome_for(source, cache, **options) as first:
        first.index()
    source.unlink()
    with genome_for(source, cache, **options) as second:
        assert second.installed()
        second.index()
        assert second.gene_name_of_gene_id("g1") == "SAVED"


def test_legacy_index_without_source_identity_is_rebuilt(tmp_path):
    source = tmp_path / "annotation.gtf"
    write_gtf(source, "FRESH")
    with genome_for(source, tmp_path / "cache") as genome:
        index = Path(genome._annotation_index_paths()[0])
        index.parent.mkdir()
        with sqlite3.connect(index) as connection:
            connection.execute('CREATE TABLE "_datacache_metadata" (version INT)')
            connection.execute(
                'INSERT INTO "_datacache_metadata" VALUES (?)', (DATABASE_SCHEMA_VERSION,)
            )
        assert not genome.installed()
        genome.index()
        assert genome.gene_name_of_gene_id("g1") == "FRESH"


def test_read_only_gtf_and_explicit_overwrite(tmp_path):
    source = tmp_path / "annotation.gtf"
    write_gtf(source, "ORIGINAL")
    source.chmod(0o444)
    try:
        with genome_for(source, tmp_path / "cache") as genome:
            genome.index()
            with sqlite3.connect(genome.db.local_db_path) as connection:
                connection.execute('UPDATE gene SET gene_name = "TAMPERED"')
            genome.index(overwrite=True)
            assert genome.gene_name_of_gene_id("g1") == "ORIGINAL"
            assert source.stat().st_mode & 0o777 == 0o444
    finally:
        source.chmod(0o600)

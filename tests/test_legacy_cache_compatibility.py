"""Existing installations remain usable without changing any cache files."""

import gzip
from pathlib import Path
from urllib.parse import urlsplit

import datacache
import pytest

from pyensembl import Genome
from pyensembl.common import dump_pickle
from pyensembl.database import Database
from pyensembl.download_cache import DownloadCache
from .test_gtf_index_source import write_gtf


def snapshot(directory):
    result = {}
    for path in directory.rglob("*"):
        info = path.lstat()
        result[str(path.relative_to(directory))] = (
            info.st_mode, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns,
            path.read_bytes() if path.is_file() else None,
        )
    return result


def legacy_genome(directory, origin, compressed=False):
    """Build the pre-2.27.5 layout directly, without download path selection."""
    directory.mkdir()
    gtf = directory / "annotation.gtf"
    transcript = directory / "transcripts.fa"
    protein = directory / "proteins.fa"
    write_gtf(gtf, "SAVED")
    transcript.write_text(">t\nACGT\n")
    protein.write_text(">p\nMAV\n")
    if compressed:
        paths = []
        for path in (gtf, transcript, protein):
            zipped = path.with_name(path.name + ".gz")
            zipped.write_bytes(gzip.compress(path.read_bytes()))
            path.unlink()
            paths.append(zipped)
        gtf, transcript, protein = paths
    settings = dict(
        reference_name="test", annotation_name="legacy",
        gtf_path_or_url=origin + gtf.name,
        transcript_fasta_paths_or_urls=[origin + transcript.name],
        protein_fasta_paths_or_urls=[origin + protein.name],
        cache_directory_path=str(directory),
    )
    with Database(str(gtf), source=settings["gtf_path_or_url"],
                  cache_directory_path=str(directory)) as database:
        database.create()
    # Build historical dictionaries directly, independent of current naming.
    for path, sequences in ((transcript, {"t": "ACGT"}), (protein, {"p": "MAV"})):
        dump_pickle(sequences, str(path) + ".pickle")
    return settings


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("symlink", [False, True])
def test_complete_legacy_cache_is_unchanged_by_offline_use(tmp_path, monkeypatch, compressed, symlink):
    directory = tmp_path / "cache"
    settings = legacy_genome(directory, "https://example.invalid/", compressed)
    if symlink:
        gtf = directory / Path(urlsplit(settings["gtf_path_or_url"]).path).name
        target = tmp_path / gtf.name
        gtf.rename(target)
        gtf.symlink_to(target)
    for path in directory.iterdir():
        path.chmod(0o444)
    directory.chmod(0o555)
    before = snapshot(directory)

    def forbidden(*args, **kwargs):
        pytest.fail("using an installed legacy cache must not acquire or rebuild data")

    monkeypatch.setattr(datacache, "fetch_file", forbidden)
    monkeypatch.setattr(Database, "create", forbidden)
    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)
    try:
        for _ in range(2):
            with Genome(**settings) as genome:
                assert genome.installed()
                assert genome.inspect_data()["installed"]
                genome.download()
                genome.index()
                assert genome.gene_name_of_gene_id("g1") == "SAVED"
                assert genome.transcript_sequence("t") == "ACGT"
                assert genome.protein_sequence("p") == "MAV"
        assert snapshot(directory) == before
    finally:
        directory.chmod(0o755)
        for path in directory.iterdir():
            path.chmod(0o644)


def test_existing_mirror_receipt_does_not_invalidate_legacy_setup(tmp_path, monkeypatch):
    source = tmp_path / "mirror" / "source.fa"
    source.parent.mkdir()
    source.write_bytes(b">t\nACGT\n")
    directory = tmp_path / "cache"
    legacy = directory / source.name
    datacache.fetch_file(source.as_uri(), destination=str(legacy), record_provenance=True)
    before = snapshot(directory)
    source.unlink()
    cache = DownloadCache("test", "legacy", cache_directory_path=str(directory))

    def forbidden(*args, **kwargs):
        pytest.fail("legacy mirror compatibility must not fetch or rewrite provenance")

    monkeypatch.setattr(datacache, "fetch_file", forbidden)
    assert cache.download_or_copy_if_necessary("https://example.invalid/source.fa", True) == str(legacy)
    assert snapshot(directory) == before


@pytest.mark.parametrize("url_name,legacy_name,decompress", [
    ("data.fa.gz", "data.fa.gz", False),
    ("data.fa.gz", "data.fa", True),
    ("data.fa.gz?download=1", "data.fa.gz_download_1", True),
])
def test_legacy_compression_and_query_filenames_still_resolve(
    tmp_path, url_name, legacy_name, decompress,
):
    legacy = tmp_path / legacy_name
    legacy.write_bytes(b"existing bytes")
    before = snapshot(tmp_path)
    cache = DownloadCache("test", "legacy", cache_directory_path=str(tmp_path),
                          decompress_on_download=decompress)
    result = cache.download_or_copy_if_necessary("https://example.invalid/" + url_name)
    assert result == str(legacy)
    assert snapshot(tmp_path) == before


def test_explicit_download_keeps_legacy_files_and_switches_loaded_genome(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    directory = tmp_path / "cache"
    settings = legacy_genome(directory, source.as_uri() + "/")
    before = snapshot(directory)
    write_gtf(source / "annotation.gtf", "NEW")
    (source / "transcripts.fa").write_text(">t\nGGCC\n")
    (source / "proteins.fa").write_text(">p\nMGM\n")
    with Genome(**settings) as genome:
        assert genome.gene_name_of_gene_id("g1") == "SAVED"
        assert genome.transcript_sequence("t") == "ACGT"
        assert genome.protein_sequence("p") == "MAV"
        genome.download(overwrite=True)
        genome.index()
        assert genome.gene_name_of_gene_id("g1") == "NEW"
        assert genome.transcript_sequence("t") == "GGCC"
        assert genome.protein_sequence("p") == "MGM"
        assert genome.gtf_path != str(directory / "annotation.gtf")
    after = snapshot(directory)
    assert {path: after[path] for path in before} == before
    with Genome(**settings) as reopened:
        assert reopened.installed()
        assert reopened.gene_name_of_gene_id("g1") == "NEW"

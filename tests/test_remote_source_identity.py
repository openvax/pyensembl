"""Downloaded annotations and sequences belong to their complete source URL."""

import gzip
from pathlib import Path

import datacache
import pytest

from pyensembl import Genome
from pyensembl.download_cache import DownloadCache, MissingRemoteFile
from .test_gtf_index_source import write_gtf


def make_cache(directory, decompress=False):
    return DownloadCache(
        "test", "remote-identity", cache_directory_path=str(directory),
        decompress_on_download=decompress,
    )


@pytest.mark.parametrize("decompress", [False, True])
def test_same_basename_downloads_remain_independent_offline(tmp_path, decompress):
    sources = [tmp_path / name / "annotation.gtf.gz" for name in ("a", "b")]
    for source, contents in zip(sources, (b"FIRST", b"SECOND")):
        source.parent.mkdir()
        source.write_bytes(gzip.compress(contents))
    cache = make_cache(tmp_path / "cache", decompress)
    first = Path(cache.download_or_copy_if_necessary(sources[0].as_uri(), True))
    assert cache.inspect(sources[1].as_uri()).status == "missing"
    with pytest.raises(MissingRemoteFile):
        cache.download_or_copy_if_necessary(sources[1].as_uri())
    second = Path(cache.download_or_copy_if_necessary(sources[1].as_uri(), True))
    assert first != second
    for path, expected in ((first, b"FIRST"), (second, b"SECOND")):
        contents = path.read_bytes()
        assert (contents if decompress else gzip.decompress(contents)) == expected

    # Overwriting one source cannot replace the other's downloaded bytes.
    sources[1].write_bytes(gzip.compress(b"REPLACED"))
    cache.download_or_copy_if_necessary(sources[1].as_uri(), True, overwrite=True)
    assert (first.read_bytes() if decompress else gzip.decompress(first.read_bytes())) == b"FIRST"
    for source in sources:
        source.unlink()
    reopened = make_cache(tmp_path / "cache", decompress)
    for source, path in zip(sources, (first, second)):
        assert reopened.download_or_copy_if_necessary(source.as_uri()) == str(path)
        assert reopened.download_or_copy_if_necessary(source.as_uri(), True) == str(path)


@pytest.mark.parametrize("other", [
    "https://other.test/a/source.fa.gz?version=one",
    "https://example.test/b/source.fa.gz?version=one",
    "http://example.test/a/source.fa.gz?version=one",
    "https://example.test/a/source.fa.gz?version=two",
    "https://user:secret@example.test/a/source.fa.gz?version=one",
    "https://example.test/a/source.fa.gz?version=one#selection",
])
def test_full_url_selects_cache_and_keeps_query_out_of_filename(tmp_path, other):
    cache = make_cache(tmp_path / "cache", decompress=True)
    source = "https://example.test/a/source.fa.gz?version=one"
    first, second = cache.cached_path(source), cache.cached_path(other)
    assert first != second
    assert first.endswith("source.fa") and second.endswith("source.fa")
    assert "version=" not in first + second
    assert "secret" not in first + second
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("compressed", [False, True])
def test_download_modes_do_not_share_a_file(tmp_path, compressed):
    source = tmp_path / ("source.fa.gz" if compressed else "source.fa")
    contents = b">t\nACGT\n"
    source.write_bytes(gzip.compress(contents) if compressed else contents)
    raw = make_cache(tmp_path / "cache")
    expanded = make_cache(tmp_path / "cache", decompress=True)
    raw_path = raw.download_or_copy_if_necessary(source.as_uri(), True)
    expanded_path = expanded.download_or_copy_if_necessary(source.as_uri(), True)
    assert raw_path != expanded_path
    assert Path(raw_path).read_bytes() == source.read_bytes()
    assert Path(expanded_path).read_bytes() == b">t\nACGT\n"


@pytest.mark.parametrize("filename", ["a" * 240 + ".fa", "é" * 110 + ".fa"])
def test_long_remote_filename_leaves_room_for_provenance_and_indexes(tmp_path, filename):
    source = tmp_path / filename
    source.write_bytes(b">t\nACGT\n")
    cache = make_cache(tmp_path / "cache")
    path = cache.download_or_copy_if_necessary(source.as_uri(), True)
    assert Path(path).read_bytes() == source.read_bytes()
    assert cache.inspect(source.as_uri()).source_url == source.as_uri()
    Path(path + ".pickle").write_bytes(b"index")


@pytest.mark.parametrize("receipt", [False, True])
def test_existing_legacy_download_is_reused_in_place(tmp_path, receipt):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nOLD\n")
    directory = tmp_path / "cache"
    directory.mkdir()
    legacy = directory / source.name
    if receipt:
        datacache.fetch_file(source.as_uri(), destination=str(legacy), record_provenance=True)
    else:
        legacy.write_bytes(source.read_bytes())
    source.write_bytes(b">t\nNEW\n")
    cache = make_cache(directory)
    assert cache.inspect(source.as_uri()).status == "available"
    assert cache.download_or_copy_if_necessary(source.as_uri()) == str(legacy)
    assert cache.download_or_copy_if_necessary(source.as_uri(), True) == str(legacy)
    assert legacy.read_bytes() == b">t\nOLD\n"


def test_genomes_isolate_remote_gtf_and_fasta_indexes(tmp_path, monkeypatch):
    options = []
    for label, sequence in (("FIRST", "ACGT"), ("SECOND", "GGCC")):
        directory = tmp_path / label
        gtf = directory / "annotation.gtf"
        write_gtf(gtf, label)
        transcript = directory / "transcript.fa"
        transcript.write_text(">t\n%s\n" % sequence)
        protein = directory / "protein.fa"
        protein.write_text(">p\n%s\n" % label)
        options.append(dict(
            reference_name="test", annotation_name="remote-identity",
            gtf_path_or_url=gtf.as_uri(),
            transcript_fasta_paths_or_urls=[transcript.as_uri()],
            protein_fasta_paths_or_urls=[protein.as_uri()],
            cache_directory_path=str(tmp_path / "cache"),
        ))
    for settings in options:
        with Genome(**settings) as genome:
            genome.download()
            genome.index()
    # Reopen both genomes: in-memory values must not mask an index collision.
    def offline(*args, **kwargs):
        pytest.fail("a downloaded genome must be usable offline")

    monkeypatch.setattr(datacache, "fetch_file", offline)
    with Genome(**options[0]) as first, Genome(**options[1]) as second:
        for genome, label, sequence in ((first, "FIRST", "ACGT"), (second, "SECOND", "GGCC")):
            assert genome.installed()
            assert genome.gene_name_of_gene_id("g1") == label
            assert genome.transcript_sequence("t") == sequence
            assert genome.protein_sequence("p") == label
        assert set(first._annotation_index_paths()).isdisjoint(second._annotation_index_paths())
        first.delete_index_files()
        assert not first.installed()
        assert second.installed()

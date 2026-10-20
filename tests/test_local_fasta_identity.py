"""Local FASTA queries never share newly created indexes across sources."""

import gzip
import os
from pathlib import Path

import pytest

from pyensembl import Genome, SequenceData
from pyensembl.common import dump_pickle
from .test_legacy_cache_compatibility import snapshot


def write_fasta(path, sequence, identifier="t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (">%s\n%s\n" % (identifier, sequence)).encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)


def genome_for(source, cache, copied=False, role="transcript"):
    return Genome(
        "test", "local-fasta", cache_directory_path=str(cache),
        copy_local_files_to_cache=copied,
        **{role + "_fasta_paths_or_urls": [str(source)]},
    )


@pytest.mark.parametrize("role", ["transcript", "protein"])
@pytest.mark.parametrize("copied", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
def test_same_basename_sources_remain_independent(tmp_path, monkeypatch, role, copied, compressed):
    name = "sequences.fa" + (".gz" if compressed else "")
    sources = [tmp_path / label / name for label in ("a", "b")]
    sequences = ("ACGT", "GGCC") if role == "transcript" else ("MAV", "MGM")
    cache = tmp_path / "cache"
    for source, sequence in zip(sources, sequences):
        write_fasta(source, sequence)
        with genome_for(source, cache, copied, role) as genome:
            genome.index()
            assert getattr(genome, role + "_sequence")("t") == sequence
    if copied:
        for source in sources:
            source.unlink()
    before = snapshot(cache)

    def forbidden(*args, **kwargs):
        pytest.fail("an indexed FASTA must reopen without acquisition or parsing")

    monkeypatch.setattr("datacache.fetch_file", forbidden)
    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)
    with genome_for(sources[0], cache, copied, role) as first, \
         genome_for(sources[1], cache, copied, role) as second:
        for genome, sequence in zip((first, second), sequences):
            assert genome.installed()
            assert genome.inspect_data()["installed"]
            genome.index()
            assert getattr(genome, role + "_sequence")("t") == sequence
        assert set(first._annotation_index_paths()).isdisjoint(second._annotation_index_paths())
        assert snapshot(cache) == before
        first.delete_index_files()
        assert not first.installed()
        assert second.installed()
        assert getattr(second, role + "_sequence")("t") == sequences[1]


@pytest.mark.parametrize("copied", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
def test_changed_bytes_with_same_size_and_mtime_select_new_index(tmp_path, copied, compressed):
    source = tmp_path / ("source.fa.gz" if compressed else "source.fa")
    cache = tmp_path / "cache"
    write_fasta(source, "ACGT")
    with genome_for(source, cache, copied) as first:
        first.index()
        assert first.transcript_sequence("t") == "ACGT"
        old_index = Path(first._annotation_index_paths()[0])
        saved = old_index.read_bytes()
        info = source.stat()
        write_fasta(source, "GGCC")
        assert source.stat().st_size == info.st_size
        os.utime(source, ns=(info.st_atime_ns, info.st_mtime_ns))
        with genome_for(source, cache, copied) as second:
            assert not second.installed()
            second.index()
            assert second.transcript_sequence("t") == "GGCC"
            assert Path(second._annotation_index_paths()[0]) != old_index
        first.index()
        assert first.transcript_sequence("t") == "GGCC"
        assert old_index.read_bytes() == saved


def test_loaded_dictionary_notices_changed_attached_source(tmp_path):
    source = tmp_path / "source.fa"
    write_fasta(source, "ACGT")
    data = SequenceData(str(source), str(tmp_path / "cache"))
    assert data.get("t") == "ACGT"
    assert "t" in data
    write_fasta(source, "GGCC", identifier="u")
    assert data.get("u") == "GGCC"
    assert data.get("t") is None
    assert "u" in data and "t" not in data


@pytest.mark.parametrize("copied", [False, True])
@pytest.mark.parametrize("legacy", [False, True])
def test_read_only_existing_setup_is_preserved(tmp_path, monkeypatch, copied, legacy):
    source = tmp_path / "source" / "source.fa"
    cache = tmp_path / "cache"
    write_fasta(source, "ACGT")
    cache.mkdir()
    if legacy:
        if copied:
            (cache / source.name).write_bytes(source.read_bytes())
        dump_pickle({"t": "ACGT"}, cache / (source.name + ".pickle"))
    else:
        with genome_for(source, cache, copied) as genome:
            genome.index()
    source.chmod(0o444)
    for path in cache.iterdir():
        path.chmod(0o444)
    cache.chmod(0o555)
    before = snapshot(tmp_path)

    def forbidden(*args, **kwargs):
        pytest.fail("read-only existing setups must not acquire or rebuild")

    monkeypatch.setattr("datacache.fetch_file", forbidden)
    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)
    try:
        for _ in range(2):
            with genome_for(source, cache, copied) as genome:
                assert genome.installed()
                assert genome.inspect_data()["installed"]
                genome.download()
                genome.index()
                assert genome.transcript_sequence("t") == "ACGT"
        assert snapshot(tmp_path) == before
    finally:
        cache.chmod(0o755)
        source.chmod(0o644)
        for path in cache.iterdir():
            path.chmod(0o644)


@pytest.mark.parametrize("refresh", ["index", "download"])
def test_refresh_of_legacy_copy_leaves_saved_files_intact(tmp_path, refresh):
    source = tmp_path / "source" / "source.fa"
    cache = tmp_path / "cache"
    write_fasta(source, "ACGT")
    cache.mkdir()
    (cache / source.name).write_bytes(source.read_bytes())
    dump_pickle({"t": "ACGT"}, cache / (source.name + ".pickle"))
    before = snapshot(cache)
    with genome_for(source, cache, copied=True) as genome:
        assert genome.transcript_sequence("t") == "ACGT"
        write_fasta(source, "GGCC")
        getattr(genome, refresh)()
        assert genome.transcript_sequence("t") == "GGCC"
        after = snapshot(cache)
        assert {path: after[path] for path in before} == before


def test_deletion_finds_only_this_source_after_it_disappears(tmp_path):
    cache = tmp_path / "cache"
    sources = [tmp_path / label / "same.fa" for label in ("a", "b")]
    for source in sources:
        write_fasta(source, "ACGT")
        with genome_for(source, cache) as genome:
            genome.index()
    sources[0].unlink()
    with genome_for(sources[0], cache) as first, genome_for(sources[1], cache) as second:
        removed = first.delete_index_files()
        assert len(removed) == 1
        assert second.installed()
        assert sources[1].read_text() == ">t\nACGT\n"


def test_long_local_filename_can_publish_index_atomically(tmp_path):
    source = tmp_path / ("é" * 110 + ".fa")
    write_fasta(source, "ACGT")
    with genome_for(source, tmp_path / "cache", copied=True) as genome:
        genome.index()
        assert genome.installed()
        assert genome.transcript_sequence("t") == "ACGT"


def test_source_changed_during_parse_is_not_published(tmp_path, monkeypatch):
    source = tmp_path / "source.fa"
    cache = tmp_path / "cache"
    write_fasta(source, "ACGT")

    def changing_parser(path, progress):
        write_fasta(source, "GGCC")
        return {"t": "ACGT"}

    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", changing_parser)
    with pytest.raises(ValueError, match="FASTA changed while indexing"):
        SequenceData(str(source), str(cache)).index()
    assert not list(cache.glob("*.pickle"))


@pytest.mark.parametrize("legacy", [False, True])
def test_decompressed_copy_reuses_index_until_original_changes(tmp_path, monkeypatch, legacy):
    source = tmp_path / "source.fa.gz"
    write_fasta(source, "ACGT")
    cache = tmp_path / "cache"
    if legacy:
        cache.mkdir()
        (cache / "source.fa").write_bytes(gzip.decompress(source.read_bytes()))
        dump_pickle({"t": "ACGT"}, cache / "source.fa.pickle")
    def make():
        return Genome("test", "local-fasta", cache_directory_path=str(cache),
                      transcript_fasta_paths_or_urls=[str(source)],
                      copy_local_files_to_cache=True, decompress_on_download=True)
    with make() as genome:
        genome.index()
        assert genome.installed()
        assert genome.transcript_sequence("t") == "ACGT"
    before = snapshot(cache)
    with monkeypatch.context() as patch:
        def forbidden(*args, **kwargs):
            pytest.fail("unchanged decompressed copies must not be copied or reindexed")
        patch.setattr("datacache.fetch_file", forbidden)
        patch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)
        with make() as restored:
            assert restored.installed()
            restored.download()
            restored.index()
            assert restored.transcript_sequence("t") == "ACGT"
        assert snapshot(cache) == before
    write_fasta(source, "GGCC")
    with make() as changed:
        changed.index()
        assert changed.installed()
        assert changed.transcript_sequence("t") == "GGCC"

"""Tests for datacache-backed safe gzipped-JSON FASTA dictionary caches (#506)."""

import gzip
import json
from pathlib import Path
import pickle

import datacache
import pytest

from pyensembl import Genome, SequenceData
from pyensembl.common import dump_pickle
from pyensembl.fasta_cache import (
    FASTA_INDEX_FORMAT,
    FASTA_INDEX_SUFFIX,
    FASTA_INDEX_VERSION,
    FastaDictionaryCache,
    cached_dictionary_lost_versions,
    file_content_digest,
    is_fasta_index_filename,
    legacy_local_index_prefix,
    load_fasta_dictionary,
    safe_filename_tail,
    save_fasta_dictionary,
)


def test_newly_written_cache_uses_safe_gzipped_json_and_never_pickles(tmp_path, monkeypatch):
    source = tmp_path / "transcripts.fa"
    source.write_text(">ENST00000001.2\nACGT\n>ENST00000002.1\nGGCC\n")
    cache_dir = tmp_path / "cache"

    def forbid_pickle(*args, **kwargs):
        pytest.fail("New FASTA dictionary caches must never use pickle")

    monkeypatch.setattr(pickle, "dump", forbid_pickle)
    monkeypatch.setattr(pickle, "load", forbid_pickle)

    sd = SequenceData([str(source)], cache_directory_path=str(cache_dir))
    sd.index()

    index_path = Path(sd.fasta_dictionary_index_paths[0])
    assert index_path.name.endswith(FASTA_INDEX_SUFFIX)
    assert index_path.parent == cache_dir
    assert not list(cache_dir.glob("*.pickle"))

    with gzip.open(index_path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload == {
        "format": FASTA_INDEX_FORMAT,
        "version": FASTA_INDEX_VERSION,
        "sequences": {
            "ENST00000001.2": "ACGT",
            "ENST00000002.1": "GGCC",
        },
    }

    # Reopen in a fresh SequenceData: reloads from .json.gz without re-parsing FASTA or unpickling.
    monkeypatch.setattr(
        "pyensembl.sequence_data.parse_fasta_dictionary",
        lambda *a, **kw: pytest.fail("Should have reused .json.gz cache"),
    )
    reopened = SequenceData([str(source)], cache_directory_path=str(cache_dir))
    assert reopened.get("ENST00000001.2") == "ACGT"
    assert reopened.get("ENST00000002.1") == "GGCC"
    assert reopened.fasta_version("ENST00000001") == 2


def test_legacy_hashed_pickle_is_reused_and_deleted_cleanly(tmp_path, monkeypatch):
    """Pre-2.41 hashed local .pickle files stay readable offline and are deleted on cleanup."""
    source = tmp_path / "proteins.fa"
    source.write_text(">ENSP00000001.1\nMSEQ\n")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    digest = file_content_digest(str(source))
    legacy_hashed = cache_dir / (
        legacy_local_index_prefix(str(source))
        + digest
        + "-"
        + safe_filename_tail(str(source), 72)
        + ".pickle"
    )
    dump_pickle({"ENSP00000001.1": "MFROMLEGACY"}, legacy_hashed)

    cache = FastaDictionaryCache(cache_dir)
    assert cache.index_path(str(source)) == str(legacy_hashed)
    assert cache.inspect(str(source)).status == "available"

    monkeypatch.setattr(
        "pyensembl.sequence_data.parse_fasta_dictionary",
        lambda *a, **kw: pytest.fail("Should have reused existing hashed legacy pickle"),
    )
    sd = SequenceData([str(source)], cache_directory_path=str(cache_dir))
    assert sd.get("ENSP00000001.1") == "MFROMLEGACY"

    deleted = cache.delete(str(source))
    assert len(deleted) == 1
    assert deleted[0][0] == str(legacy_hashed)
    assert not legacy_hashed.exists()


@pytest.mark.parametrize(
    "bad_payload",
    [
        b"not gzip bytes",
        gzip.compress(b"{not valid json"),
        gzip.compress(json.dumps({"format": "wrong", "version": 1, "sequences": {}}).encode()),
        gzip.compress(json.dumps({"format": FASTA_INDEX_FORMAT, "version": 99, "sequences": {}}).encode()),
        gzip.compress(json.dumps({"format": FASTA_INDEX_FORMAT, "version": 1, "sequences": ["not-dict"]}).encode()),
        gzip.compress(json.dumps({"format": FASTA_INDEX_FORMAT, "version": 1, "sequences": {"k": 123}}).encode()),
    ],
)
def test_corrupt_json_cache_is_rejected_and_rebuilt(tmp_path, bad_payload):
    source = tmp_path / "transcripts.fa"
    source.write_text(">ENST00000001.3\nTTTT\n")
    cache_dir = tmp_path / "cache"
    cache = FastaDictionaryCache(cache_dir)
    target = Path(cache.target_path(str(source)))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(bad_payload)

    with pytest.raises(datacache.FileValidationError):
        load_fasta_dictionary(target)

    sd = SequenceData([str(source)], cache_directory_path=str(cache_dir))
    assert sd.get("ENST00000001.3") == "TTTT"
    assert load_fasta_dictionary(target) == {"ENST00000001.3": "TTTT"}


def test_corrupt_legacy_pickle_is_rebuilt_as_safe_json(tmp_path):
    source = tmp_path / "transcripts.fa"
    source.write_text(">ENST00000001.3\nACGT\n")
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    bad_pickle = cache_dir / "transcripts.fa.pickle"
    bad_pickle.write_bytes(b"corrupt pickle bytes")

    with pytest.raises(datacache.FileValidationError):
        load_fasta_dictionary(bad_pickle)

    sd = SequenceData([str(source)], cache_directory_path=str(cache_dir))
    sd.index()
    assert sd.get("ENST00000001.3") == "ACGT"
    rebuilt = Path(sd.fasta_dictionary_index_paths[0])
    assert rebuilt.name.endswith(FASTA_INDEX_SUFFIX)
    assert load_fasta_dictionary(rebuilt) == {"ENST00000001.3": "ACGT"}


def test_remote_long_filename_cache_publishes_atomically(tmp_path):
    cache_dir = tmp_path / "cache"
    fasta = cache_dir / ("a" * 200 + ".fa")
    cache_dir.mkdir()
    fasta.write_text(">t\nACGT\n")
    url = "https://example.test/" + fasta.name

    cache = FastaDictionaryCache(cache_dir)
    saved = cache.save({"t": "ACGT"}, str(fasta), source=url)
    assert len(Path(saved).name.encode("utf-8")) <= 180
    assert load_fasta_dictionary(saved) == {"t": "ACGT"}
    assert cache.index_path(str(fasta), source=url) == saved
    assert cache.index_files(str(fasta), source=url) == [saved]


def test_save_and_helpers_contract(tmp_path):
    assert is_fasta_index_filename("genes.fa.json.gz")
    assert is_fasta_index_filename("genes.fa.pickle")
    assert not is_fasta_index_filename("genes.fa")

    target = tmp_path / "nested" / "out.json.gz"
    save_fasta_dictionary({"a": "ACGT"}, target)
    assert load_fasta_dictionary(target) == {"a": "ACGT"}

    with pytest.raises(datacache.FileValidationError):
        save_fasta_dictionary({"a": 123}, target)

    fasta = tmp_path / "v.fa"
    fasta.write_text(">ENST00000001.2\nACGT\n")
    assert cached_dictionary_lost_versions(str(fasta), {"ENST00000001": "ACGT"}) is True
    assert cached_dictionary_lost_versions(str(fasta), {"ENST00000001.2": "ACGT"}) is False


def test_genome_inspect_and_delete_use_json_fasta_cache(tmp_path):
    source = tmp_path / "transcripts.fa"
    source.write_text(">t1\nACGT\n")
    cache_dir = tmp_path / "cache"
    genome = Genome(
        "synthetic",
        "json-cache",
        transcript_fasta_paths_or_urls=[str(source)],
        cache_directory_path=str(cache_dir),
    )
    assert genome.annotation_status() == "not indexed"
    genome.index()
    assert genome.annotation_status() == "indexed"
    report = genome.inspect_data()
    assert report["files"]["transcript_fasta_1_index"].path.endswith(FASTA_INDEX_SUFFIX)
    assert report["files"]["transcript_fasta_1_index"].status == "available"
    assert genome.transcript_sequence("t1") == "ACGT"

    deleted = genome.delete_index_files()
    assert len(deleted) == 1
    assert deleted[0][0].endswith(FASTA_INDEX_SUFFIX)
    assert genome.annotation_status() == "not indexed"

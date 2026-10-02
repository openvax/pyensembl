"""Real file acquisition and offline inspection through datacache."""

import gzip
import json
from pathlib import Path

import datacache
import pytest

from pyensembl import EnsemblRelease, Genome
from pyensembl.download_cache import DownloadCache
from .test_cli_output import complete_database, ensembl_files
from .test_genome_fasta import DNA, run_cli


def make_cache(tmp_path, decompress=False):
    return DownloadCache(
        "test", "imports", copy_local_files_to_cache=True,
        decompress_on_download=decompress,
        cache_directory_path=str(tmp_path / "cache"),
    )


def test_imported_file_survives_removal_of_original(tmp_path):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nACGT\n")
    cache = make_cache(tmp_path)
    imported = cache.download_or_copy_if_necessary(str(source))
    source.unlink()
    assert cache.download_or_copy_if_necessary(str(source)) == imported
    assert Path(imported).read_bytes() == b">t\nACGT\n"


def test_compressed_local_import_honors_decompression(tmp_path):
    source = tmp_path / "source.fa.gz"
    source.write_bytes(gzip.compress(b">t\nACGT\n"))
    cache = make_cache(tmp_path, decompress=True)
    imported = cache.download_or_copy_if_necessary(str(source))
    assert imported.endswith("source.fa")
    assert Path(imported).read_bytes() == b">t\nACGT\n"
    assert gzip.decompress(source.read_bytes()) == b">t\nACGT\n"


@pytest.mark.parametrize("bad_source", ["empty", "directory"])
def test_readiness_rejects_invalid_source_even_with_an_index(tmp_path, bad_source):
    source = tmp_path / "source.fa"
    if bad_source == "empty":
        source.touch()
    else:
        source.mkdir()
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "source.fa.pickle").write_bytes(b"old index")
    genome = Genome(
        "test", "imports", transcript_fasta_paths_or_urls=[str(source)],
        cache_directory_path=str(cache),
    )
    assert not genome.installed()
    assert not genome.required_local_files_exist()
    assert genome._annotation_status() == "invalid"


def test_empty_files_ok_does_not_accept_directories(tmp_path):
    source = tmp_path / "source.fa"
    source.touch()
    genome = Genome("test", "imports", transcript_fasta_paths_or_urls=[str(source)],
                    cache_directory_path=str(tmp_path / "cache"))
    assert genome.required_local_files_exist(empty_files_ok=True)
    source.unlink()
    source.mkdir()
    assert not genome.required_local_files_exist(empty_files_ok=True)


def test_failed_import_overwrite_preserves_destination_and_receipt(tmp_path):
    source = tmp_path / "source.fa.gz"
    source.write_bytes(gzip.compress(b">t\nACGT\n"))
    cache = make_cache(tmp_path, decompress=True)
    imported = Path(cache.download_or_copy_if_necessary(str(source)))
    before = cache.inspect(str(source))
    source.write_bytes(b"not gzip")
    with pytest.raises((OSError, ValueError)):
        cache.download_or_copy_if_necessary(str(source), overwrite=True)
    assert imported.read_bytes() == b">t\nACGT\n"
    assert cache.inspect(str(source)) == before


@pytest.mark.parametrize("decompress", [False, True])
def test_remote_acquisition_and_provenance(tmp_path, decompress):
    source = tmp_path / "remote.fa.gz"
    source.write_bytes(gzip.compress(b">t\nACGT\n"))
    cache = make_cache(tmp_path, decompress)
    imported = cache.download_or_copy_if_necessary(source.as_uri(), download_if_missing=True)
    assert Path(imported).read_bytes() == (b">t\nACGT\n" if decompress else source.read_bytes())
    inspection = cache.inspect(source.as_uri())
    assert inspection.status == "available"
    assert inspection.source_url == source.as_uri()
    assert inspection.fetched_at
    # Without a trusted expected digest, datacache avoids a checksum pass.
    assert inspection.recorded_sha256 is None
    assert not inspection.verified
    source.unlink()
    # Both an offline lookup and a normal download call must reuse valid data.
    assert cache.download_or_copy_if_necessary(source.as_uri()) == imported
    assert cache.download_or_copy_if_necessary(source.as_uri(), download_if_missing=True) == imported


def test_invalid_cached_download_requires_explicit_overwrite(tmp_path):
    source = tmp_path / "remote.fa"
    source.write_bytes(b">t\nACGT\n")
    cache = make_cache(tmp_path)
    destination = Path(cache.cached_path(source.as_uri()))
    destination.parent.mkdir()
    destination.touch()
    for download in (False, True):
        with pytest.raises(datacache.FileValidationError):
            cache.download_or_copy_if_necessary(source.as_uri(), download_if_missing=download)
    cache.download_or_copy_if_necessary(source.as_uri(), download_if_missing=True, overwrite=True)
    assert destination.read_bytes() == source.read_bytes()


@pytest.mark.parametrize("receipt", [None, "{broken", '{"url": "old"}'])
def test_legacy_or_bad_provenance_does_not_break_inspection(tmp_path, receipt):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nACGT\n")
    if receipt is not None:
        (tmp_path / ".source.fa.datacache.json").write_text(receipt)
    cache = DownloadCache("test", "imports", cache_directory_path=str(tmp_path / "absent"))
    inspection = cache.inspect(str(source))
    assert inspection.status == "available"
    assert inspection.source_url is None
    assert not inspection.verified
    assert not Path(cache.cache_directory_path).exists()


def test_stale_receipt_is_not_reported_as_current_provenance(tmp_path):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nACGT\n")
    cache = make_cache(tmp_path)
    imported = cache.download_or_copy_if_necessary(str(source))
    Path(imported).write_bytes(b">changed\nGGCC\n")
    inspection = cache.inspect(str(source))
    assert inspection.status == "available"
    assert inspection.source_url is None
    assert not inspection.verified


def test_inventory_is_read_only_and_reports_database_completeness(tmp_path, monkeypatch):
    source = tmp_path / "source.gtf"
    source.write_bytes(b"data")
    cache = tmp_path / "cache"
    database = cache / "source.db"
    complete_database(database)
    genome = Genome("test", "imports", gtf_path_or_url=str(source),
                    transcript_fasta_paths_or_urls=["https://example.org/transcripts.fa"],
                    cache_directory_path=str(cache))
    before = {str(path): path.stat().st_mtime_ns for path in tmp_path.rglob("*")}

    def forbidden(*args, **kwargs):
        pytest.fail("inspection must not acquire, create or index")

    monkeypatch.setattr(datacache, "fetch_file", forbidden)
    monkeypatch.setattr(datacache, "ensure_dir", forbidden)
    monkeypatch.setattr(genome, "index", forbidden)
    report = genome.inspect_data()
    assert report["annotation"] == "incomplete"
    assert not report["installed"]
    assert report["reference_dna"] is None
    assert report["files"]["gtf"].status == "available"
    assert report["files"]["gtf_index"].status == "available"
    assert report["files"]["transcript_fasta_1"].status == "missing"
    assert {str(path): path.stat().st_mtime_ns for path in tmp_path.rglob("*")} == before
    database.write_bytes(b"broken SQLite")
    report = genome.inspect_data()
    assert report["files"]["gtf_index"].status == "corrupt"
    assert "rebuild" in str(report["files"]["gtf_index"].error)


def test_unreadable_source_is_reported_and_cannot_be_used(tmp_path):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nACGT\n")
    genome = Genome("test", "imports", transcript_fasta_paths_or_urls=[str(source)],
                    cache_directory_path=str(tmp_path / "cache"))
    source.chmod(0)
    try:
        assert genome._annotation_status() == "inaccessible"
        assert not genome.installed()
        inspection = genome.inspect_data()["files"]["transcript_fasta_1"]
        assert isinstance(inspection.error, PermissionError)
        with pytest.raises(PermissionError):
            genome.download_cache.download_or_copy_if_necessary(str(source))
    finally:
        source.chmod(0o600)


def test_cli_inspects_multiple_releases_offline_as_json(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    location = ensembl_files(93)
    (location / Path(EnsemblRelease(93).gtf_url).name).write_bytes(b"")
    run_cli(monkeypatch, "inspect", "--release", "93", "94", "--json")
    reports = json.loads(capsys.readouterr().out)
    assert [report["annotation_version"] for report in reports] == [93, 94]
    assert reports[0]["cache_directory"] == str(location)
    assert [report["annotation"] for report in reports] == ["invalid", "missing"]
    assert reports[0]["files"]["gtf"]["status"] == "corrupt"
    assert "empty" in reports[0]["files"]["gtf"]["error"]
    assert not (tmp_path / "pyensembl" / "GRCh38" / "ensembl94").exists()


def test_cli_custom_inspection_and_human_output(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    source = tmp_path / "attached.fa"
    source.write_bytes(b">t\nACGT\n")
    run_cli(monkeypatch, "inspect", "--reference-name", "custom",
            "--annotation-name", "mine", "--transcript-fasta", str(source))
    output = capsys.readouterr().out
    assert "annotation=not indexed" in output
    assert "transcript_fasta_1" in output
    assert "available" in output
    assert str(source) in output
    assert not (tmp_path / "pyensembl").exists()


def test_json_option_only_applies_to_inspection(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, "list", "--json")
    assert "--json requires inspect" in capsys.readouterr().err


def test_missing_inspection_does_not_create_cache(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("inspection must not acquire or create files")

    monkeypatch.setattr(datacache, "fetch_file", forbidden)
    monkeypatch.setattr(datacache, "ensure_dir", forbidden)
    directory = tmp_path / "absent"
    genome = Genome("test", "imports", gtf_path_or_url="https://example.org/source.gtf",
                    cache_directory_path=str(directory))
    assert genome.inspect_data()["annotation"] == "missing"
    assert not directory.exists()


def test_import_can_reuse_a_source_already_in_its_cache(tmp_path):
    cache = make_cache(tmp_path)
    source = tmp_path / "cache" / "source.fa"
    source.parent.mkdir()
    source.write_bytes(b">t\nACGT\n")
    assert cache.download_or_copy_if_necessary(str(source), overwrite=True) == str(source)
    assert source.read_bytes() == b">t\nACGT\n"


def test_inspection_never_deserializes_pickle_contents(tmp_path):
    source = tmp_path / "source.fa"
    source.write_bytes(b">t\nACGT\n")
    index = tmp_path / "cache" / "source.fa.pickle"
    index.parent.mkdir()
    index.write_bytes(b"not a pickle, but an available file")
    genome = Genome("test", "imports", transcript_fasta_paths_or_urls=[str(source)],
                    cache_directory_path=str(index.parent))
    assert genome.inspect_data()["annotation"] == "indexed"
    index.unlink()
    index.mkdir()
    assert genome.inspect_data()["annotation"] == "not indexed"
    assert genome.inspect_data()["files"]["transcript_fasta_1_index"].status == "corrupt"


def test_dna_inspection_preserves_fingerprint_readiness_and_never_rebuilds(tmp_path):
    source = tmp_path / "dna.fa"
    source.write_bytes(DNA)
    with Genome("test", "imports", genome_fasta_path_or_url=str(source),
                cache_directory_path=str(tmp_path / "cache")) as genome:
        report = genome.inspect_data()
        assert report["reference_dna"] == "needs index"
        assert not report["installed"]
        genome.index_genome_fasta()
        report = genome.inspect_data(check_genome_fasta=True)
        assert report["reference_dna"] == "indexed"
        assert report["installed"]
        source.write_bytes(DNA + b"\n")
        before = {str(path): path.stat().st_mtime_ns for path in tmp_path.rglob("*")}
        report = genome.inspect_data()
        assert report["reference_dna"] == "needs index"
        assert report["files"]["genome_fasta_index"].status == "corrupt"
        assert not genome.installed()
        assert {str(path): path.stat().st_mtime_ns for path in tmp_path.rglob("*")} == before


def test_unreadable_dna_is_reported_without_crashing(tmp_path):
    source = tmp_path / "dna.fa"
    source.write_bytes(DNA)
    genome = Genome("test", "imports", genome_fasta_path_or_url=str(source),
                    cache_directory_path=str(tmp_path / "cache"))
    source.chmod(0)
    try:
        report = genome.inspect_data()
        assert report["reference_dna"] == "inaccessible"
        assert report["files"]["genome_fasta"].status == "inaccessible"
        assert not genome.installed()
        assert not genome.required_local_files_exist()
    finally:
        source.chmod(0o600)

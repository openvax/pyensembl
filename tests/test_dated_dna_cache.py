"""Dated DNA shares only verified artifacts and preserves private installs."""

from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil

import pytest
from datacache import FileValidationError

from pyensembl import EnsemblAnnotation, EnsemblRelease, prune_genome_fastas
from pyensembl import genome_fasta_cache as cache
from pyensembl.genome_fasta import GenomeFasta
from .test_genome_fasta import DNA, forbid_downloads, serve_downloads
from .test_legacy_cache_compatibility import snapshot


def dated(date="2026_04", explicit=False, **kwargs):
    if explicit:
        return EnsemblAnnotation("GCA_000001405.29", date, genome_fasta=True, **kwargs)
    return EnsemblRelease(date, genome_fasta=True, **kwargs)


@pytest.fixture
def platform(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "data"))
    state = {"payload": gzip.compress(DNA, mtime=0), "manifest": None, "size": None}
    calls = serve_downloads(monkeypatch, tmp_path, lambda url: state["payload"])

    def metadata(request, **kwargs):
        url = request if isinstance(request, str) else request.full_url
        if url.endswith("/md5sum.txt"):
            checksum = hashlib.md5(state["payload"], usedforsecurity=False).hexdigest()
            manifest = state["manifest"]
            if manifest is None:
                manifest = "\n".join(checksum + "  " + mask + ".fa.bgz"
                                     for mask in ("unmasked", "softmasked", "hardmasked"))
            return io.BytesIO(manifest.encode())
        assert request.get_method() == "HEAD"
        response = io.BytesIO()
        response.headers = {"Content-Length": str(state["size"] or len(state["payload"]))}
        return response

    monkeypatch.setattr(cache, "urlopen", metadata)
    return state, calls


@pytest.mark.parametrize("explicit", [False, True])
def test_dates_share_download_and_index_and_reopen_read_only(platform, monkeypatch, explicit):
    _, calls = platform
    first, second = dated("2025_12", explicit), dated("2026_04", explicit)
    for genome in (first, second):
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        assert genome.sequence("MT", 1, 4) == "GCTA"
        genome.close()
    assert len(calls) == 1
    assert first.genome_fasta_path == second.genome_fasta_path
    assert first._genome_fasta.index_path == second._genome_fasta.index_path
    assert cache.canonical_source_options(first._genome_fasta.source) == ("toplevel", "none")

    root = cache.dna_cache_root().parent
    paths = list(root.rglob("*")) + [root]
    for path in paths:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(root)

    def offline(*args, **kwargs):
        pytest.fail("installed dated DNA must reuse files and indexes without writes or network")

    monkeypatch.setattr(cache, "urlopen", offline)
    monkeypatch.setattr("pyensembl.genome_fasta._write_json", offline)
    monkeypatch.setattr(cache, "_write_json", offline)
    forbid_downloads(monkeypatch, "installed dated DNA must work offline")
    try:
        for date in ("2025_12", "2026_04"):
            with dated(date, explicit) as genome:
                genome.download_genome_fasta()
                genome.index_genome_fasta()
                assert genome.sequence("1", 3, 10) == "GTNNTTAA"
        assert snapshot(root) == before
    finally:
        for path in reversed(paths):
            path.chmod(0o755 if path.is_dir() else 0o644)


def test_release_and_explicit_accession_share_despite_labels(platform):
    _, calls = platform
    with dated() as first, dated(explicit=True, reference_name="custom label") as second:
        first.download_genome_fasta()
        second.download_genome_fasta()
        assert first.genome_fasta_path == second.genome_fasta_path
        assert len(calls) == 1


@pytest.mark.parametrize("change", ["bytes", "compression", "mask", "provider", "version", "gcf"])
def test_distinct_artifacts_never_alias(platform, change):
    state, calls = platform
    first = dated()
    first.download_genome_fasta()
    if change == "bytes":
        state["payload"] = gzip.compress(DNA.replace(b"GCTA", b"CCAA"), mtime=0)
    elif change == "compression":
        state["payload"] = gzip.compress(DNA, mtime=1)
    second = dated("2025_12", genome_fasta_mask="soft" if change == "mask" else "none")
    if change in ("provider", "version", "gcf"):
        accession = {"version": "GCA_000001405.28", "gcf": "GCF_000001405.29"}.get(
            change, "GCA_000001405.29"
        )
        second = EnsemblAnnotation(accession, "2026_04", genome_fasta=True,
                                   provider="community" if change == "provider" else "ensembl")
    try:
        second.download_genome_fasta()
        assert first.genome_fasta_path != second.genome_fasta_path
        assert len(calls) == 2
        assert first.sequence("MT", 1, 4) == "GCTA"
        assert second.sequence("MT", 1, 4) == ("CCAA" if change == "bytes" else "GCTA")
    finally:
        first.close()
        second.close()


@pytest.mark.parametrize("manifest", ["", "bad  unmasked.fa.bgz", "a" * 32 + "  other.fa.bgz",
                                       ("a" * 32 + "  unmasked.fa.bgz\n") * 2])
def test_incomplete_or_ambiguous_metadata_isolates_by_url(platform, manifest):
    state, calls = platform
    state["manifest"] = manifest
    with dated("2025_12") as first, dated() as second:
        for genome in (first, second):
            genome.download_genome_fasta()
            assert genome._genome_fasta.identity == {"source": genome.genome_fasta_urls[0]}
        assert first.genome_fasta_path != second.genome_fasta_path
        assert len(calls) == 2


def test_unreachable_metadata_isolates_by_url(platform, monkeypatch):
    def unavailable(*args, **kwargs):
        raise OSError("offline metadata")

    monkeypatch.setattr(cache, "urlopen", unavailable)
    with dated() as genome:
        genome.download_genome_fasta()
        assert "source" in genome._genome_fasta.identity


@pytest.mark.parametrize("bad", ["md5", "size"])
def test_mismatched_download_is_not_published_or_registered(platform, bad):
    state, _ = platform
    if bad == "md5":
        state["manifest"] = "0" * 32 + "  unmasked.fa.bgz\n"
    else:
        state["size"] = len(state["payload"]) + 1
    with dated() as genome:
        with pytest.raises(FileValidationError):
            genome.download_genome_fasta()
        assert genome.genome_fasta_path is None
        assert not genome._genome_fasta.reference_path.exists()
        assert not genome._genome_fasta.index_path.exists()


def test_failed_refresh_preserves_installed_reference_and_bytes(platform):
    state, _ = platform
    with dated() as genome:
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        directory = genome._genome_fasta.directory
        before = snapshot(directory)
        reference = genome._genome_fasta.reference_path.read_bytes()
        state["manifest"] = "0" * 32 + "  unmasked.fa.bgz\n"
        with pytest.raises(FileValidationError):
            genome.download_genome_fasta(overwrite=True)
        assert snapshot(directory) == before
        assert genome._genome_fasta.reference_path.read_bytes() == reference
    with dated() as reopened:
        assert reopened.sequence("MT", 1, 4) == "GCTA"


@pytest.mark.parametrize("explicit", [False, True])
def test_existing_private_install_stays_unchanged(platform, monkeypatch, explicit):
    selection = dated(explicit=explicit)
    old = GenomeFasta(selection._genome_fasta.source, selection.download_cache.cache_directory_path)
    old.prepare(download=True)
    old.open()
    old.remember()
    old.close()
    root = Path(selection.download_cache.cache_directory_path)
    paths = list(root.rglob("*")) + [root]
    for path in paths:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(root)

    def forbidden(*args, **kwargs):
        pytest.fail("legacy private DNA must not migrate, reindex, acquire or rewrite metadata")

    monkeypatch.setattr(cache, "urlopen", forbidden)
    monkeypatch.setattr(cache, "_write_json", forbidden)
    monkeypatch.setattr("pyensembl.genome_fasta._write_json", forbidden)
    forbid_downloads(monkeypatch, "preserve legacy private installs")
    try:
        with dated(explicit=explicit) as restored:
            assert type(restored._genome_fasta) is GenomeFasta
            restored.download_genome_fasta()
            restored.index_genome_fasta()
            assert restored.sequence("MT", 1, 4) == "GCTA"
            assert restored.genome_fasta_path == old.installed_path
        assert snapshot(root) == before
        assert not cache.dna_cache_root().exists()
    finally:
        for path in reversed(paths):
            path.chmod(0o755 if path.is_dir() else 0o644)


def test_concurrent_dates_share_one_download_and_prune_only_when_unused(platform):
    _, calls = platform
    genomes = [dated("2025_12"), dated()]

    def install(genome):
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        genome.close()

    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(install, genomes))
    assert len(calls) == 1
    directory = genomes[0]._genome_fasta.directory
    shutil.rmtree(genomes[0].download_cache.cache_directory_path)
    assert prune_genome_fastas() == []
    shutil.rmtree(genomes[1].download_cache.cache_directory_path)
    assert [path for path, _ in prune_genome_fastas(dry_run=True)] == [str(directory)]
    assert directory.exists()
    assert [path for path, _ in prune_genome_fastas()] == [str(directory)]
    assert not directory.exists()


def test_malformed_dated_identity_fails_closed(platform):
    with dated() as genome:
        genome.download_genome_fasta()
        ref = genome._genome_fasta.reference_path
        state = json.loads(ref.read_text())
        state["identity"]["annotation_provider"] = "../../outside"
        state["shared_key"] = cache._identity_key(state["identity"])
        ref.write_text(json.dumps(state))
        with pytest.raises(ValueError, match="Invalid dated"):
            dated()
        with pytest.raises(ValueError, match="Cannot safely prune"):
            prune_genome_fastas()


@pytest.mark.parametrize("replacement", ["http://ftp.ebi.ac.uk", "https://example.test"])
def test_mirrors_do_not_share(platform, replacement):
    with dated(server=replacement + "/pub/ensemblorganisms") as genome:
        assert type(genome._genome_fasta) is GenomeFasta
        genome.download_genome_fasta()
        assert not cache.dna_cache_root().exists()

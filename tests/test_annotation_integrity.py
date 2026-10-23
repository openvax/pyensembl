"""Real HTTP transfers, checksum pinning, and read-only legacy cache reuse."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from functools import partial
import gzip
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import pickle
import sys
from threading import Thread
from urllib.parse import urlsplit

import datacache
import pytest

from pyensembl import EnsemblRelease, annotation_integrity, shell
from pyensembl.download_cache import DownloadCache
from .test_legacy_cache_compatibility import snapshot
from .test_release_cache_combinations import TinySources


@pytest.fixture
def origin(tmp_path, monkeypatch):
    root = tmp_path / "origin"
    root.mkdir()
    calls = []

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            calls.append(self.path)
            super().do_GET()

    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(Handler, directory=str(root)))
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = "http://127.0.0.1:%d" % server.server_port
    monkeypatch.setattr("pyensembl.ensembl_release.ENSEMBL_PLATFORM_FTP_SERVER", base)
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    genome = EnsemblRelease("2026_04")
    sources = TinySources(root)
    expected = sources.register(genome, "CHECKED", alt=True)
    for url, content in sources.payloads.items():
        path = root / urlsplit(url).path.lstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    manifest = root / urlsplit(genome.gtf_url).path.lstrip("/")
    manifest = manifest.parent / "md5sum.txt"

    def publish_manifest():
        manifest.write_text("\n".join(
            "%s  %s" % (hashlib.md5(p.read_bytes(), usedforsecurity=False).hexdigest(), p.name)
            for p in sorted(manifest.parent.iterdir()) if p != manifest
        ) + "\n")

    publish_manifest()
    try:
        yield genome, root, manifest, publish_manifest, calls, expected
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_download_verifies_and_reopens_without_network_or_writes(origin, monkeypatch):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    genome.index()
    assert genome.transcript_by_id(expected["transcript"]).sequence == expected["cdna"]
    assert {r["status"] for r in genome.verify_sources(True)} == {"verified"}
    saved = genome.to_json()
    cache = Path(genome.download_cache.cache_directory_path)
    before = snapshot(cache)
    count = len(calls)

    def forbidden(*args, **kwargs):
        pytest.fail("Installed data must reopen offline")

    monkeypatch.setattr(annotation_integrity, "urlopen", forbidden)
    for restored in (EnsemblRelease("2026_04"), EnsemblRelease.from_json(saved),
                     pickle.loads(pickle.dumps(genome))):
        restored.download()
        assert restored.installed()
        assert restored.gene_by_id(expected["gene"]).name == "CHECKED"
        assert {r["status"] for r in restored.verify_sources()} == {"verified"}
        assert len(restored.to_dict()["_resolved"]["source_checksums"]) == 3
        restored.close()
    assert snapshot(cache) == before
    assert len(calls) == count


@pytest.mark.parametrize("legacy", [False, True])
def test_existing_unrecorded_cache_stays_unmodified(origin, legacy, monkeypatch):
    genome, root, manifest, publish, calls, expected = origin
    cache = genome.download_cache
    for url in genome._annotation_sources().values():
        destination = Path(cache._cache_path(url, legacy=legacy))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((root / urlsplit(url).path.lstrip("/")).read_bytes())
    genome.index()
    before = snapshot(Path(cache.cache_directory_path))
    monkeypatch.setattr(annotation_integrity, "urlopen", lambda *a, **kw: pytest.fail("network"))
    genome.download()
    assert genome.installed()
    assert {r["status"] for r in genome.verify_sources()} == {"unrecorded"}
    assert "source_checksums" not in genome.to_dict()["_resolved"]
    assert snapshot(Path(cache.cache_directory_path)) == before
    assert calls == []


def test_bad_transfer_does_not_publish_and_can_retry(origin):
    genome, root, manifest, publish, calls, expected = origin
    source = root / urlsplit(genome.gtf_url).path.lstrip("/")
    original = source.read_bytes()
    source.write_bytes(gzip.compress(b"corrupt GTF transfer"))
    with pytest.raises(datacache.FileValidationError, match="MD5 mismatch"):
        genome.download()
    destination = Path(genome.download_cache.local_path(genome.gtf_url))
    assert not destination.exists()
    assert not annotation_integrity.receipt_path(destination).exists()
    source.write_bytes(original)
    genome.download()
    assert {r["status"] for r in genome.verify_sources()} == {"verified"}


@pytest.mark.parametrize("manifest_text", ["", "bad checksum  genes-including_alt.gtf.gz\n", "duplicate"])
def test_missing_or_ambiguous_checksum_fails_before_download(origin, manifest_text):
    genome, root, manifest, publish, calls, expected = origin
    if manifest_text == "duplicate":
        manifest_text = manifest.read_text() * 2
    manifest.write_text(manifest_text)
    with pytest.raises(annotation_integrity.AnnotationIntegrityError, match="Expected one checksum"):
        genome.download()
    assert all(call.endswith("md5sum.txt") for call in calls)
    assert not Path(genome.download_cache.local_path(genome.gtf_url)).exists()


def test_local_and_upstream_changes_are_separate_and_pins_survive_json(origin, tmp_path):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    saved = json.loads(genome.to_json())
    source = root / urlsplit(genome.gtf_url).path.lstrip("/")
    source.write_bytes(gzip.compress(b"changed upstream"))
    publish()
    report = genome.verify_sources(True)[0]
    assert report["status"] == "verified"
    assert report["upstream_status"] == "changed"
    local = Path(genome.gtf_path)
    local.write_bytes(gzip.compress(b"changed locally"))
    before = snapshot(local.parent)
    assert genome.verify_sources()[0]["status"] == "changed"
    assert snapshot(local.parent) == before
    restored = EnsemblRelease.from_json(json.dumps(saved))
    # Move the serialized selection to an empty cache without changing its pins.
    state = restored.to_dict()
    state["cache_directory_path"] = str(tmp_path / "empty-cache")
    fresh = EnsemblRelease.from_dict(state)
    with pytest.raises(annotation_integrity.AnnotationIntegrityError, match="Upstream checksum changed"):
        fresh.download()
    assert not Path(fresh.download_cache.local_path(fresh.gtf_url)).exists()
    fresh.download(overwrite=True)
    assert fresh.verify_sources(True)[0]["status"] == "verified"
    # The original pin remains authoritative even after a different object refreshes.
    assert restored.verify_sources()[0]["status"] == "changed"


def test_failed_overwrite_preserves_source_receipt_and_index(origin):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    genome.index()
    before = snapshot(Path(genome.download_cache.cache_directory_path))
    source = root / urlsplit(genome.gtf_url).path.lstrip("/")
    source.write_bytes(gzip.compress(b"bad replacement"))
    with pytest.raises(datacache.FileValidationError):
        genome.download(overwrite=True)
    # Lock timestamps may change during an explicit overwrite; saved data and
    # indexes must retain their exact inode, mode, timestamps, and bytes.
    after = snapshot(Path(genome.download_cache.cache_directory_path))
    for name, state in before.items():
        if not name.endswith(".lock"):
            assert after[name] == state, name
    assert genome.gene_by_id(expected["gene"]).name == "CHECKED"


def test_concurrent_downloads_publish_one_consistent_receipt(origin):
    genome, root, manifest, publish, calls, expected = origin
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda _: EnsemblRelease("2026_04").download(), range(2)))
    assert {r["status"] for r in genome.verify_sources(True)} == {"verified"}
    for url in genome._annotation_sources().values():
        assert calls.count(urlsplit(url).path) == 1


def test_cli_verifies_read_only_and_returns_failure_for_changed_sources(origin, monkeypatch, capsys):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "verify", "--release", "2026_04", "--json"])
    before = snapshot(Path(genome.download_cache.cache_directory_path))
    shell.run()
    assert {s["status"] for s in json.loads(capsys.readouterr().out)[0]["sources"]} == {"verified"}
    assert snapshot(Path(genome.download_cache.cache_directory_path)) == before
    Path(genome.gtf_path).write_bytes(b"damaged")
    with pytest.raises(SystemExit) as error:
        shell.run()
    assert error.value.code == 1
    assert json.loads(capsys.readouterr().out)[0]["sources"][0]["status"] == "changed"


def test_cli_verification_cannot_report_success_for_inaccessible_files(origin, monkeypatch, capsys):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    inspect = DownloadCache.inspect

    def inaccessible(cache, source, *args, **kwargs):
        return replace(inspect(cache, source, *args, **kwargs), status="inaccessible",
                       error=PermissionError("shared volume is not readable"))

    monkeypatch.setattr(DownloadCache, "inspect", inaccessible)
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "verify", "--release", "2026_04", "--json"])
    with pytest.raises(SystemExit) as error:
        shell.run()
    assert error.value.code == 1
    assert {s["status"] for s in json.loads(capsys.readouterr().out)[0]["sources"]} == {"inaccessible"}


def test_broken_receipt_is_distinct_from_a_legacy_source_without_receipt(origin):
    genome, root, manifest, publish, calls, expected = origin
    genome.download()
    annotation_integrity.receipt_path(genome.gtf_path).write_text("broken receipt")
    reopened = EnsemblRelease("2026_04")
    before = snapshot(Path(genome.download_cache.cache_directory_path))
    assert reopened.verify_sources()[0]["status"] == "invalid-receipt"
    # A previously serialized acquisition still has its independent pin.
    assert genome.verify_sources()[0]["status"] == "verified"
    reopened.download()
    assert snapshot(Path(genome.download_cache.cache_directory_path)) == before

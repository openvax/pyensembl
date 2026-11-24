"""Tests for OPENVAX_DATA_CACHE fallback, pyensembl dirs, and --offline / OPENVAX_OFFLINE."""

import json
import os
import sys

import pytest

from pyensembl.cache_location import cache_root, resolve_cache_location
from pyensembl import shell
from pyensembl.cli import install, logging_config


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in (
        "PYENSEMBL_CACHE_DIR",
        "OPENVAX_DATA_CACHE",
        "DATACACHE_OFFLINE",
        "OPENVAX_OFFLINE",
        "DATACACHE_NO_PROGRESS",
        "OPENVAX_NO_PROGRESS",
    ):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(logging_config, "configure_logging", lambda **kw: None)


def test_openvax_data_cache_fallback_behind_pyensembl_cache_dir(tmp_path, monkeypatch):
    shared = tmp_path / "shared"
    specific = tmp_path / "specific"

    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    loc = resolve_cache_location()
    assert loc.path == str(shared / "pyensembl")
    assert loc.origin == "shared_env"
    assert loc.detail == "OPENVAX_DATA_CACHE"
    assert cache_root() == str(shared / "pyensembl")
    assert not shared.exists()

    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(specific))
    loc = resolve_cache_location()
    assert loc.path == str(specific / "pyensembl")
    assert loc.origin == "env"
    assert loc.detail == "PYENSEMBL_CACHE_DIR"
    assert cache_root() == str(specific / "pyensembl")


def test_dirs_cli_readonly_and_json(tmp_path, monkeypatch, capsys):
    shared = tmp_path / "shared"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")

    # Missing directory: read-only, creates nothing
    shell.run(["dirs", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["cache_root"] == str(shared / "pyensembl")
    assert payload["origin"] == "shared_env"
    assert payload["detail"] == "OPENVAX_DATA_CACHE"
    assert payload["exists"] is False
    assert payload["offline"] is True
    assert payload["offline_source"] == "OPENVAX_OFFLINE=1"
    assert not shared.exists()

    # Populated directory: reports subdirectories and human sizes
    genome_dir = shared / "pyensembl" / "GRCh38" / "ensembl110"
    genome_dir.mkdir(parents=True)
    (genome_dir / "data.gtf").write_bytes(b"a" * 2048)

    shell.run(["dirs"])
    out = capsys.readouterr().out
    assert "OPENVAX_DATA_CACHE" in out
    assert "GRCh38" in out
    assert "2.0 KiB" in out
    assert "Offline:" in out


def test_offline_flag_blocks_remote_download_and_restores_env(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    with pytest.raises(SystemExit) as exc:
        shell.run([
            "--offline",
            "install",
            "--reference-name", "GRCh38",
            "--annotation-name", "test",
            "--gtf", "https://example.org/test.gtf.gz",
        ])
    assert exc.value.code == 1
    assert "Offline mode is enabled" in capsys.readouterr().err
    assert "DATACACHE_OFFLINE" not in os.environ


def test_progress_respects_openvax_no_progress(monkeypatch):
    monkeypatch.setattr(sys.stderr, "isatty", lambda: True)
    assert install.progress_available() is True
    monkeypatch.setenv("OPENVAX_NO_PROGRESS", "1")
    assert install.progress_available() is False

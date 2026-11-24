"""Tests for the per-action modules behind ``pyensembl.cli.main`` and ``pyensembl.shell.run``."""

import json

import pytest

from pyensembl import EnsemblRelease, shell
from pyensembl.annotation_integrity import AnnotationIntegrityError
from pyensembl.cli import available, delete, install, logging_config, main, parser, selection, verify
from pyensembl.dated_releases import AVAILABLE_EXPIRE_AFTER
from pyensembl.download_cache import MissingLocalFile


def test_validate_arguments_accepts_matching_flags():
    parser.validate_arguments("inspect", True, False, False, [75])
    parser.validate_arguments("available", False, True, False, [])
    parser.validate_arguments("verify", False, False, True, [])
    parser.validate_arguments("list", True, False, False, [])


@pytest.mark.parametrize("action, json_output, refresh, check_remote, releases, message", [
    ("install", False, True, False, [], "--refresh requires available"),
    ("list", False, False, True, [], "--check-remote requires verify"),
    ("delete-all-files", False, False, False, ["latest"], "concrete --release"),
])
def test_validate_arguments_rejects_misplaced_flags(
    action, json_output, refresh, check_remote, releases, message
):
    with pytest.raises(ValueError, match=message):
        parser.validate_arguments(action, json_output, refresh, check_remote, releases)


def test_run_prune_reports_candidates(monkeypatch):
    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", lambda dry_run: [("/a", 1234)])
    assert delete.run_prune(True) == "Would delete /a (1,234 bytes)"
    assert delete.run_prune(False) == "Deleted /a (1,234 bytes)"


def test_run_prune_reports_nothing_and_propagates_errors(monkeypatch):
    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", lambda dry_run: [])
    assert delete.run_prune(False) == "No orphan genome FASTAs"

    def fail(dry_run):
        raise ValueError("unsafe")

    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", fail)
    with pytest.raises(ValueError, match="unsafe"):
        delete.run_prune(False)


def test_selects_catalogue():
    assert not available.selects_catalogue()
    assert available.selects_catalogue(json_output=True)
    assert available.selects_catalogue(include_alt=False)
    assert available.selects_catalogue(species=["human"])


def test_run_available_chooses_listing(monkeypatch):
    monkeypatch.setattr(available, "available_catalogue", lambda **kwargs: "catalogue")
    monkeypatch.setattr("pyensembl.dated_releases.fetch_all_dated_releases", lambda expire_after: expire_after)
    monkeypatch.setattr(
        available,
        "format_available_species",
        lambda dated_releases: "species %s" % dated_releases,
    )
    assert available.run_available(provider="ensembl") == "catalogue"
    assert available.run_available() == "species %s" % AVAILABLE_EXPIRE_AFTER
    assert available.run_available(refresh=True) == "species 0"


@pytest.mark.parametrize("selector", [
    ["--species", "human"], ["--assembly", "GRCh37"], ["--species", "dog"],
])
def test_available_catalogue_text_lists_the_json_releases(
    tmp_path, monkeypatch, capsys, selector
):
    """The text table shows every release --json does, including numbered
    releases, which are ints (#507)."""
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(logging_config, "configure_logging", lambda **kw: None)
    main.main(["available", "--json"] + selector)
    expected = [
        (str(row["release"]), row["latest"])
        for row in json.loads(capsys.readouterr().out)
    ]
    assert any(release.isdigit() for release, _ in expected)
    shell.run(["available"] + selector)
    header, *rows = capsys.readouterr().out.splitlines()
    assert header.split()[:2] == ["Release", "Assembly"]
    assert [(row.split()[0], row.endswith(" yes")) for row in rows] == expected


def test_run_available_conflicting_selectors():
    with pytest.raises(ValueError, match="one assembly selector"):
        available.run_available(reference_name="GRCh37", assembly="GRCh38")


class FakeRelease(EnsemblRelease):
    def __init__(self, statuses):
        self.release = 75  # numbered releases are ints, as in EnsemblRelease
        self.statuses = statuses

    def verify_sources(self, check_remote):
        return [dict(source="gtf", status=s) for s in self.statuses]


def test_run_verify_reports_status():
    text, ok = verify.run_verify([FakeRelease(["verified"])])
    assert ok
    assert text.splitlines()[1].split() == ["75", "gtf", "verified", "not", "checked"]
    text, ok = verify.run_verify([FakeRelease(["mismatch"])], json_output=True)
    assert not ok and '"mismatch"' in text


def test_run_verify_rejects_non_releases():
    with pytest.raises(ValueError, match="requires an Ensembl release"):
        verify.run_verify([object()])


def test_run_install_passes_options(monkeypatch):
    calls = []
    monkeypatch.setattr(install, "install_genome", lambda *a, **kw: calls.append((a, kw)))
    monkeypatch.setattr(install, "progress_available", lambda: False)
    install.run_install("genome", only_genome_fasta=True, overwrite=True)
    assert calls == [(("genome", True, True), {"show_progress": False})]


def test_run_install_wraps_missing_files(monkeypatch):
    def missing(*args, **kwargs):
        raise MissingLocalFile("/nowhere/x.gtf")

    monkeypatch.setattr(install, "install_genome", missing)
    with pytest.raises(install.InstallError, match="file not found: /nowhere/x.gtf"):
        install.run_install("genome")


def test_require_ensembl_releases():
    verify.require_ensembl_releases([FakeRelease(["verified"])])
    with pytest.raises(ValueError, match="requires an Ensembl release"):
        verify.require_ensembl_releases([FakeRelease([]), object()])


def test_run_reports_only_argument_mistakes_as_cli_errors(monkeypatch, capsys):
    """``run`` turns selector conflicts and non-release ``verify`` targets into
    argument errors (exit 2); listing and checksum failures keep tracebacks."""
    monkeypatch.setattr(logging_config, "configure_logging", lambda **kw: None)

    def malformed_listing(**kwargs):
        raise ValueError("malformed listing")

    monkeypatch.setattr(available, "run_available", malformed_listing)
    with pytest.raises(ValueError, match="malformed listing"):
        main.main(["available"])
    with pytest.raises(SystemExit) as exit_info:
        main.main(["available", "--json"])
    assert exit_info.value.code == 2
    assert "malformed listing" in capsys.readouterr().err

    class AmbiguousManifestRelease(FakeRelease):
        def verify_sources(self, check_remote):
            raise AnnotationIntegrityError("Expected one checksum")

    monkeypatch.setattr(
        selection,
        "collect_selected_genomes",
        lambda args: [AmbiguousManifestRelease([])],
    )
    with pytest.raises(AnnotationIntegrityError, match="Expected one checksum"):
        main.main(["verify", "--release", "75", "--check-remote"])
    monkeypatch.setattr(selection, "collect_selected_genomes", lambda args: [object()])
    with pytest.raises(SystemExit) as exit_info:
        main.main(["verify", "--release", "75"])
    assert exit_info.value.code == 2
    assert "requires an Ensembl release" in capsys.readouterr().err

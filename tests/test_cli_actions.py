"""Tests for the per-action functions behind ``pyensembl.shell.run``."""

import json
from types import SimpleNamespace

import pytest

from pyensembl import EnsemblRelease, shell
from pyensembl.annotation_integrity import AnnotationIntegrityError
from pyensembl.dated_releases import AVAILABLE_EXPIRE_AFTER
from pyensembl.download_cache import MissingLocalFile


def test_validate_arguments_accepts_matching_flags():
    shell.validate_arguments("inspect", True, False, False, [75])
    shell.validate_arguments("available", False, True, False, [])
    shell.validate_arguments("verify", False, False, True, [])
    shell.validate_arguments("list", True, False, False, [])


@pytest.mark.parametrize("action, json_output, refresh, check_remote, releases, message", [
    ("install", False, True, False, [], "--refresh requires available"),
    ("list", False, False, True, [], "--check-remote requires verify"),
    ("delete-all-files", False, False, False, ["latest"], "concrete --release"),
])
def test_validate_arguments_rejects_misplaced_flags(
        action, json_output, refresh, check_remote, releases, message):
    with pytest.raises(ValueError, match=message):
        shell.validate_arguments(action, json_output, refresh, check_remote, releases)


def test_run_prune_reports_candidates(monkeypatch):
    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", lambda dry_run: [("/a", 1234)])
    assert shell.run_prune(True) == "Would delete /a (1,234 bytes)"
    assert shell.run_prune(False) == "Deleted /a (1,234 bytes)"


def test_run_prune_reports_nothing_and_propagates_errors(monkeypatch):
    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", lambda dry_run: [])
    assert shell.run_prune(False) == "No orphan genome FASTAs"

    def fail(dry_run):
        raise ValueError("unsafe")

    monkeypatch.setattr("pyensembl.genome_fasta_cache.prune_genome_fastas", fail)
    with pytest.raises(ValueError, match="unsafe"):
        shell.run_prune(False)


def make_available_args(**overrides):
    values = dict(species=[], assembly=None, assembly_accession=None, provider=None,
                  reference_name=None, include_alt=None, json=False, refresh=False)
    values.update(overrides)
    return SimpleNamespace(**values)


def test_selects_catalogue():
    assert not shell.selects_catalogue(make_available_args())
    assert shell.selects_catalogue(make_available_args(json=True))
    assert shell.selects_catalogue(make_available_args(include_alt=False))
    assert shell.selects_catalogue(make_available_args(species=["human"]))


def test_run_available_chooses_listing(monkeypatch):
    monkeypatch.setattr(shell, "available_catalogue", lambda args: "catalogue")
    monkeypatch.setattr("pyensembl.dated_releases.fetch_all_dated_releases", lambda expire_after: expire_after)
    monkeypatch.setattr(shell, "format_available_species",
                        lambda dated_releases: "species %s" % dated_releases)
    assert shell.run_available(make_available_args(provider="ensembl")) == "catalogue"
    assert shell.run_available(make_available_args()) == "species %s" % AVAILABLE_EXPIRE_AFTER
    assert shell.run_available(make_available_args(refresh=True)) == "species 0"


@pytest.mark.parametrize("selector", [
    ["--species", "human"], ["--assembly", "GRCh37"], ["--species", "dog"],
])
def test_available_catalogue_text_lists_the_json_releases(
        tmp_path, monkeypatch, capsys, selector):
    """The text table shows every release --json does, including numbered
    releases, which are ints (#507)."""
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    shell.run(["available", "--json"] + selector)
    expected = [(str(row["release"]), row["latest"])
                for row in json.loads(capsys.readouterr().out)]
    assert any(release.isdigit() for release, _ in expected)
    shell.run(["available"] + selector)
    header, *rows = capsys.readouterr().out.splitlines()
    assert header.split()[:2] == ["Release", "Assembly"]
    assert [(row.split()[0], row.endswith(" yes")) for row in rows] == expected


def test_run_available_conflicting_selectors():
    args = make_available_args(reference_name="GRCh37", assembly="GRCh38")
    with pytest.raises(ValueError, match="one assembly selector"):
        shell.run_available(args)


class FakeRelease(EnsemblRelease):
    def __init__(self, statuses):
        self.release = 75  # numbered releases are ints, as in EnsemblRelease
        self.statuses = statuses

    def verify_sources(self, check_remote):
        return [dict(source="gtf", status=s) for s in self.statuses]


def test_run_verify_reports_status():
    text, ok = shell.run_verify([FakeRelease(["verified"])])
    assert ok
    assert text.splitlines()[1].split() == ["75", "gtf", "verified", "not", "checked"]
    text, ok = shell.run_verify([FakeRelease(["mismatch"])], json_output=True)
    assert not ok and '"mismatch"' in text


def test_run_verify_rejects_non_releases():
    with pytest.raises(ValueError, match="requires an Ensembl release"):
        shell.run_verify([object()])


def test_run_install_passes_options(monkeypatch):
    calls = []
    monkeypatch.setattr(shell, "_install", lambda *a, **kw: calls.append((a, kw)))
    monkeypatch.setattr(shell, "_progress_available", lambda: False)
    shell.run_install("genome", only_genome_fasta=True, overwrite=True)
    assert calls == [(("genome", True, True), {"show_progress": False})]


def test_run_install_wraps_missing_files(monkeypatch):
    def missing(*args, **kwargs):
        raise MissingLocalFile("/nowhere/x.gtf")

    monkeypatch.setattr(shell, "_install", missing)
    with pytest.raises(shell.InstallError, match="file not found: /nowhere/x.gtf"):
        shell.run_install("genome")


def test_require_ensembl_releases():
    shell.require_ensembl_releases([FakeRelease(["verified"])])
    with pytest.raises(ValueError, match="requires an Ensembl release"):
        shell.require_ensembl_releases([FakeRelease([]), object()])


def test_run_reports_only_argument_mistakes_as_cli_errors(monkeypatch, capsys):
    """``run`` turns selector conflicts and non-release ``verify`` targets into
    argument errors (exit 2); listing and checksum failures keep tracebacks."""
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)

    def malformed_listing(args):
        raise ValueError("malformed listing")

    monkeypatch.setattr(shell, "run_available", malformed_listing)
    with pytest.raises(ValueError, match="malformed listing"):
        shell.run(["available"])
    with pytest.raises(SystemExit) as exit_info:
        shell.run(["available", "--json"])
    assert exit_info.value.code == 2
    assert "malformed listing" in capsys.readouterr().err

    class AmbiguousManifestRelease(FakeRelease):
        def verify_sources(self, check_remote):
            raise AnnotationIntegrityError("Expected one checksum")

    monkeypatch.setattr(shell, "collect_selected_genomes",
                        lambda args: [AmbiguousManifestRelease([])])
    with pytest.raises(AnnotationIntegrityError, match="Expected one checksum"):
        shell.run(["verify", "--release", "75", "--check-remote"])
    monkeypatch.setattr(shell, "collect_selected_genomes", lambda args: [object()])
    with pytest.raises(SystemExit) as exit_info:
        shell.run(["verify", "--release", "75"])
    assert exit_info.value.code == 2
    assert "requires an Ensembl release" in capsys.readouterr().err

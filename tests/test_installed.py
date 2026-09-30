"""Genome.installed() only reads the cache, and release selection uses it (#398)."""

from weakref import WeakValueDictionary

import datacache
import pytest

from pyensembl import EnsemblRelease, genome_for_reference_name
from pyensembl.reference_name import find_species_by_reference
from .test_cli_output import ensembl_files, write


@pytest.fixture
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    # Releases memoized elsewhere point at the real cache.
    monkeypatch.setattr(EnsemblRelease, "_genome_cache", WeakValueDictionary())
    return tmp_path


def fail_downloads(monkeypatch):
    def download(*args, **kwargs):
        pytest.fail("installed() must not download")

    monkeypatch.setattr(datacache, "fetch_file", download)
    monkeypatch.setattr("pyensembl.genome_fasta.fetch_file", download)


def test_installed_follows_each_install_step(cache):
    genome = EnsemblRelease(93)
    assert not genome.installed()
    ensembl_files(93, indexes=False)
    assert not genome.installed()
    directory = ensembl_files(93, database=write)  # Interrupted database build.
    assert not genome.installed()
    for database in directory.glob("*.db"):
        database.unlink()
    ensembl_files(93)
    assert genome.installed()


def test_installed_checks_reference_dna_when_configured(cache):
    ensembl_files(93)
    assert EnsemblRelease(93).installed()
    assert not EnsemblRelease(93, genome_fasta=True).installed()


def test_installed_never_downloads_or_creates_files(cache, monkeypatch):
    fail_downloads(monkeypatch)
    for genome in (EnsemblRelease(93), EnsemblRelease(93, genome_fasta=True)):
        assert not genome.installed()
    assert list(cache.iterdir()) == []


def newest_grch38_release():
    species = find_species_by_reference("GRCh38")
    return species.reference_assemblies["GRCh38"][1]


def test_reference_name_prefers_installed_releases(cache, monkeypatch):
    fail_downloads(monkeypatch)
    ensembl_files(93)
    # Neither newer release is ready: one has only downloads, the other an
    # unfinished database, as in #398.
    ensembl_files(112, indexes=False)
    ensembl_files(110, database=write)
    assert genome_for_reference_name("GRCh38").release == 93
    newest = genome_for_reference_name("GRCh38", allow_older_downloaded_release=False)
    assert newest.release == newest_grch38_release()


def test_reference_name_falls_back_to_downloads_then_the_newest_release(cache):
    ensembl_files(93, indexes=False)
    newer = ensembl_files(112, indexes=False)
    assert genome_for_reference_name("GRCh38").release == 112
    for path in newer.iterdir():
        path.unlink()
    assert genome_for_reference_name("GRCh38").release == 93
    for path in cache.glob("pyensembl/GRCh38/*/*"):
        path.unlink()
    assert genome_for_reference_name("GRCh38").release == newest_grch38_release()

"""One cache root on every platform, chosen by datacache.

Windows is simulated through whichever of appdirs and platformdirs is
installed, so these tests pass whichever one datacache uses to find the
platform cache directory.
"""

import importlib
from pathlib import Path
import sys

import datacache
import pytest

from pyensembl import EnsemblRelease
from pyensembl.download_cache import DownloadCache, cache_root, cache_subdirectory
from pyensembl.genome_fasta_cache import SharedGenomeFasta, dna_cache_root
from .test_cli_output import complete_database
from .test_genome_fasta import list_rows


@pytest.fixture
def windows(tmp_path, monkeypatch):
    monkeypatch.delenv("PYENSEMBL_CACHE_DIR", raising=False)
    local_app_data = tmp_path / "Local"
    folder = lambda name: str(local_app_data)  # noqa: E731
    try:
        platformdirs = importlib.import_module("platformdirs")
        windows = importlib.import_module("platformdirs.windows")
        monkeypatch.setattr(platformdirs, "PlatformDirs", windows.Windows)
        monkeypatch.setattr(windows, "get_win_folder", folder)
    except ImportError:
        pass
    try:
        appdirs = importlib.import_module("appdirs")
        monkeypatch.setattr(appdirs, "system", "win32")
        monkeypatch.setattr(appdirs, "_get_win_folder", folder, raising=False)
    except ImportError:
        pass
    return local_app_data


def test_windows_nests_every_genome_under_one_root(windows):
    root = Path(cache_root())
    assert root == windows / "pyensembl" / "pyensembl" / "Cache"
    release = EnsemblRelease(81)
    assert Path(release.download_cache.cache_directory_path) == root / "GRCh38" / "ensembl81"
    assert dna_cache_root() == root / "dna_cache"
    # Releases beside dna_cache share Ensembl DNA.
    assert type(EnsemblRelease(81, genome_fasta=True)._genome_fasta) is SharedGenomeFasta


def test_windows_releases_installed_by_older_versions_stay_put(windows):
    legacy = Path(datacache.get_data_dir(subdir=cache_subdirectory("GRCh38", "ensembl", 82)))
    assert legacy.parent.parent != Path(cache_root())  # The old, unrelated layout.
    legacy.mkdir(parents=True)
    assert Path(EnsemblRelease(82).download_cache.cache_directory_path) == legacy
    new = Path(EnsemblRelease(83).download_cache.cache_directory_path)
    assert new.parent.parent == Path(cache_root())


def test_windows_list_shows_custom_genomes(windows, monkeypatch, capsys):
    complete_database(Path(cache_root()) / "GRCm38" / "custom81" / "annotation.db")
    rows = list_rows(monkeypatch, capsys)
    assert rows["custom81"]["Species"] == "custom"
    assert rows["custom81"]["Annotation"] == "indexed"


@pytest.mark.skipif(sys.platform == "win32", reason="checks Linux and macOS paths")
@pytest.mark.parametrize(
    "genome", [("GRCh38", "ensembl", 81), ("GRCm38", "custom", None), (None, None, None)]
)
def test_linux_and_macos_caches_do_not_move(monkeypatch, genome):
    monkeypatch.delenv("PYENSEMBL_CACHE_DIR", raising=False)
    old_location = datacache.get_data_dir(subdir=cache_subdirectory(*genome))
    assert DownloadCache(*genome).cache_directory_path == old_location

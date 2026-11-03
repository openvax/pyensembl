"""Tab completion through argcomplete for the ``pyensembl`` console script (#514)."""

import os
from pathlib import Path
import subprocess
import sys

from pyensembl import shell

REPOSITORY = Path(shell.__file__).resolve().parents[1]


def complete(line, cache, tmp_path):
    """Ask the console script for completions of ``line``, as a shell would.

    Parameters
    ----------
    line : str
        The command line typed so far, e.g. ``"pyensembl --"``.
    cache : pathlib.Path
        Value for ``PYENSEMBL_CACHE_DIR``.
    tmp_path : pathlib.Path
        Directory for argcomplete's output file.

    Returns
    -------
    list of str
        The completions argcomplete printed.
    """
    command = Path(sys.executable).parent / ("pyensembl.exe" if os.name == "nt" else "pyensembl")
    output = tmp_path / "completions"
    environment = dict(
        os.environ,
        PYENSEMBL_CACHE_DIR=str(cache),
        # Use this checkout, also when tests run from an extracted source archive.
        PYTHONPATH=str(REPOSITORY),
        _ARGCOMPLETE="1",
        COMP_LINE=line,
        COMP_POINT=str(len(line)),
        _ARGCOMPLETE_IFS="\n",
        _ARGCOMPLETE_STDOUT_FILENAME=str(output),
    )
    result = subprocess.run([str(command)], env=environment, capture_output=True,
                            text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return [word.strip() for word in output.read_text().split("\n") if word.strip()]


def make_cache(root):
    """Create empty genome directories the way pyensembl lays them out.

    Parameters
    ----------
    root : pathlib.Path
        Value for ``PYENSEMBL_CACHE_DIR``; the cache is ``root / "pyensembl"``.

    Returns
    -------
    pathlib.Path
        ``root``.
    """
    for directory in (
        "GRCh38/ensembl93", "GRCh38/ensembl110", "GRCh37/ensembl75",
        "ROS_Cfam_1.0/ensembl2020_10", "GRCh38/ensembl-ensembl-GCA_000001405.29-primary-2026_04",
        "GRCm38/custom81", "custom/unannotated", "dna_cache/homo_sapiens",
    ):
        (root / "pyensembl" / directory).mkdir(parents=True)
    (root / "pyensembl" / "GRCh38" / "ensembl999.json").write_text("{}")
    return root


def test_options_and_actions_complete(tmp_path):
    options = complete("pyensembl --", tmp_path / "cache", tmp_path)
    for option in ("--release", "--species", "--reference-name", "--with-genome-fasta", "--json"):
        assert option in options
    actions = complete("pyensembl ", tmp_path / "cache", tmp_path)
    for action in ("install", "list", "inspect", "verify", "available",
                   "delete-index-files", "delete-all-files", "prune"):
        assert action in actions
    assert complete("pyensembl ins", tmp_path / "cache", tmp_path) == ["install", "inspect"]


def test_release_completes_installed_releases(tmp_path):
    cache = make_cache(tmp_path / "cache")
    assert complete("pyensembl list --release ", cache, tmp_path) == [
        "75", "93", "110", "2020_10", "2026_04",
    ]
    assert complete("pyensembl delete-all-files --release 9", cache, tmp_path) == ["93"]
    assert complete("pyensembl list --release ", tmp_path / "empty", tmp_path) == []


def test_species_completes_registered_names(tmp_path):
    assert complete("pyensembl install --species hu", tmp_path / "cache", tmp_path) == ["human"]
    names = complete("pyensembl install --species ", tmp_path / "cache", tmp_path)
    assert {"human", "homo_sapiens", "mouse", "mus_musculus"} <= set(names)


def test_installed_releases_reads_directory_names_only(tmp_path):
    root = make_cache(tmp_path)
    assert shell.installed_releases(str(root / "pyensembl")) == [
        "75", "93", "110", "2020_10", "2026_04",
    ]
    assert shell.installed_releases(str(tmp_path / "missing")) == []


def test_complete_species_lists_common_and_latin_names():
    names = shell.complete_species("")
    assert names == sorted(set(names))
    assert {"human", "homo_sapiens", "dog", "canis_familiaris"} <= set(names)


def test_cli_module_carries_the_argcomplete_marker():
    # Global activation completes only scripts whose first 1024 bytes carry it.
    assert b"PYTHON_ARGCOMPLETE_OK" in Path(shell.__file__).read_bytes()[:1024]

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


def test_options_and_actions_complete(tmp_path):
    options = complete("pyensembl --", tmp_path / "cache", tmp_path)
    for option in ("--release", "--species", "--reference-name", "--with-genome-fasta", "--json"):
        assert option in options
    actions = complete("pyensembl ", tmp_path / "cache", tmp_path)
    for action in ("install", "list", "inspect", "verify", "available",
                   "delete-index-files", "delete-all-files", "prune"):
        assert action in actions
    assert complete("pyensembl ins", tmp_path / "cache", tmp_path) == ["install", "inspect"]


def test_species_completes_registered_names(tmp_path):
    assert complete("pyensembl install --species hu", tmp_path / "cache", tmp_path) == ["human"]
    names = complete("pyensembl install --species ", tmp_path / "cache", tmp_path)
    assert {"human", "homo_sapiens", "mouse", "mus_musculus"} <= set(names)


def test_complete_species_lists_common_and_latin_names():
    names = shell.complete_species("")
    assert names == sorted(set(names))
    assert {"human", "homo_sapiens", "dog", "canis_familiaris"} <= set(names)


def test_cli_module_carries_the_argcomplete_marker():
    # Global activation completes only scripts whose first 1024 bytes carry it.
    assert b"PYTHON_ARGCOMPLETE_OK" in Path(shell.__file__).read_bytes()[:1024]

import logging
import os
import subprocess
import sys

import pytest

from pyensembl import EnsemblRelease, shell
from pyensembl.shell import (
    all_combinations_of_ensembl_genomes,
    configure_logging,
    format_available_species,
    format_installed_genomes,
    parser,
)
from .common import eq_


def test_genome_selection_grch38():
    args = parser.parse_args(["install", "--release", "100", "--species", "human"])
    genomes = all_combinations_of_ensembl_genomes(args)
    assert len(genomes) == 1
    genome = genomes[0]
    eq_(genome.species.latin_name, "homo_sapiens")
    eq_(genome.release, 100)


def test_available_action_parses():
    args = parser.parse_args(["available"])
    eq_(args.action, "available")


def test_format_available_species_includes_human_and_assemblies():
    output = format_available_species(use_color=False)
    # human is registered with common name "human" and three reference assemblies
    assert "homo_sapiens" in output
    assert "human" in output
    assert "GRCh38" in output
    assert "GRCh37" in output
    # mouse should also appear
    assert "mus_musculus" in output
    assert "GRCm38" in output


def test_format_available_species_grouped_by_division():
    output = format_available_species(use_color=False)
    # Every populated division should have a section header.
    assert "── Vertebrates " in output
    assert "── Invertebrates " in output
    assert "── Plants " in output
    assert "── Fungi " in output
    # Section ordering: Vertebrates before Invertebrates before Plants before Fungi.
    v = output.index("── Vertebrates ")
    i = output.index("── Invertebrates ")
    p = output.index("── Plants ")
    f = output.index("── Fungi ")
    assert v < i < p < f
    # Yeast is now classified as fungi; drosophila/c. elegans as metazoa.
    assert output.index("yeast") > f
    drosophila_pos = output.index("drosophila")
    assert i < drosophila_pos < p


def test_format_available_species_no_color_has_no_escape_codes():
    output = format_available_species(use_color=False)
    assert "\x1b[" not in output


def test_format_available_species_collapses_single_release():
    # NCBI36 only exists in Ensembl release 54; verify it renders as "54"
    # rather than "54–54".
    output = format_available_species(use_color=False)
    assert "NCBI36" in output
    ncbi36_line = next(
        line for line in output.splitlines() if "NCBI36" in line
    )
    assert "54–54" not in ncbi36_line
    assert "54" in ncbi36_line


# Regression test for https://github.com/openvax/pyensembl/issues/362:
# importing pyensembl / pyensembl.shell must not reconfigure the root logger
# or disable loggers the host application created before the import. Run in a
# fresh interpreter because logging state is process-global and modules are
# only imported once.
_IMPORT_SIDE_EFFECT_PROBE = """
import logging

created_before = logging.getLogger("created_before")

import pyensembl
import pyensembl.shell

# pyensembl's logging.conf attaches a CRITICAL-level StreamHandler to the root
# logger; it must not be applied merely by importing the package.
root = logging.getLogger()
pyensembl_root_handlers = [
    h for h in root.handlers if getattr(h, "level", None) == logging.CRITICAL
]
assert not pyensembl_root_handlers, (
    "import added pyensembl's console handler to the root logger: %r"
    % (pyensembl_root_handlers,)
)

# A logger created before the import must not be disabled
# (fileConfig(disable_existing_loggers=True) would have disabled it).
assert created_before.disabled is False, "import disabled a pre-existing logger"

# The package logger should carry a NullHandler so library use neither emits
# output nor triggers "No handlers could be found" warnings.
assert any(
    isinstance(h, logging.NullHandler)
    for h in logging.getLogger("pyensembl").handlers
), "pyensembl package logger is missing a NullHandler"

print("ok")
"""


def test_import_does_not_reconfigure_root_logger():
    result = subprocess.run(
        [sys.executable, "-c", _IMPORT_SIDE_EFFECT_PROBE],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("ok")


# Regression tests for https://github.com/openvax/pyensembl/issues/388:
# `pyensembl list` printed raw reprs, counted unindexed downloads as
# installed, and printed nothing for an empty cache.


def _touch(path, contents="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(contents)


@pytest.fixture
def list_cli_cache(monkeypatch, tmp_path):
    """Fake PYENSEMBL_CACHE_DIR with one fully indexed release and one
    downloaded-but-unindexed release (only the four .gz source files)."""
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    indexed = EnsemblRelease(100, species="human")
    unindexed = EnsemblRelease(101, species="mouse")
    for genome in (indexed, unindexed):
        for path in genome.required_local_files():
            _touch(path, "fake-source")
    for path in indexed._index_file_paths():
        _touch(path, "fake-index")
    return indexed, unindexed


def test_index_files_exist(list_cli_cache):
    indexed, unindexed = list_cli_cache
    assert indexed.index_files_exist()
    assert not unindexed.index_files_exist()


def test_format_installed_genomes_table(list_cli_cache):
    indexed, unindexed = list_cli_cache
    output = format_installed_genomes([indexed, unindexed], use_color=False)
    lines = output.splitlines()
    # header row carries the column names in the style of `pyensembl available`
    header = lines[0]
    for column in ("Species", "Assembly", "Release", "Status", "Path"):
        assert column in header
    # no raw Python reprs anywhere in the output
    assert "EnsemblRelease(" not in output
    # species common names, assemblies, and releases are shown
    assert "human" in output
    assert indexed.reference_name in output
    assert "mouse" in output
    assert unindexed.reference_name in output
    indexed_row = next(line for line in lines if indexed.reference_name in line)
    unindexed_row = next(line for line in lines if unindexed.reference_name in line)
    assert str(indexed.release) in indexed_row
    assert str(unindexed.release) in unindexed_row
    # cache directory is shown
    assert str(list_cli_cache[0].download_cache.cache_directory_path) in output
    # data rows are aligned with each other
    assert len({len(line) for line in lines[2:]}) == 1


def test_format_installed_genomes_marks_not_indexed(list_cli_cache):
    indexed, unindexed = list_cli_cache
    output = format_installed_genomes([indexed, unindexed], use_color=False)
    lines = output.splitlines()
    indexed_row = next(line for line in lines if indexed.reference_name in line)
    unindexed_row = next(line for line in lines if unindexed.reference_name in line)
    assert "indexed" in indexed_row
    assert "not indexed" not in indexed_row
    assert "not indexed" in unindexed_row


def test_format_installed_genomes_empty_cache():
    output = format_installed_genomes([], use_color=False)
    assert "No Ensembl genomes are installed yet." in output
    assert "pyensembl install" in output
    assert "pyensembl available" in output


def _run_list_action(monkeypatch, capsys):
    monkeypatch.setattr(shell, "configure_logging", lambda: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "list"])
    shell.run()
    return capsys.readouterr().out


def test_list_action_prints_table(list_cli_cache, monkeypatch, capsys):
    output = _run_list_action(monkeypatch, capsys)
    assert "Species" in output
    assert "Assembly" in output
    assert "Release" in output
    assert "not indexed" in output
    assert "EnsemblRelease(" not in output


def test_list_action_empty_cache_prints_friendly_message(
    monkeypatch, capsys, tmp_path
):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    output = _run_list_action(monkeypatch, capsys)
    assert "No Ensembl genomes are installed yet." in output


def test_configure_logging_preserves_existing_loggers():
    # configure_logging() applies logging.conf, which mutates process-global
    # logging state (root + pyensembl loggers). Snapshot and restore it so this
    # test doesn't leak a live console handler into sibling tests.
    root = logging.getLogger()
    pyensembl_logger = logging.getLogger("pyensembl")
    saved_root_handlers = root.handlers[:]
    saved_root_level = root.level
    saved_pyensembl_handlers = pyensembl_logger.handlers[:]
    try:
        created_before = logging.getLogger("test_configure_logging_preexisting")
        created_before.disabled = False
        configure_logging()
        # The CLI entrypoint applies logging.conf, but with
        # disable_existing_loggers=False so it leaves other loggers alone.
        assert created_before.disabled is False
        # pyensembl's own logger should be wired up to a handler for CLI output.
        assert pyensembl_logger.handlers
    finally:
        root.handlers[:] = saved_root_handlers
        root.level = saved_root_level
        pyensembl_logger.handlers[:] = saved_pyensembl_handlers

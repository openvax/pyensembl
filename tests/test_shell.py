import logging
import subprocess
import sys

from pyensembl.shell import (
    all_combinations_of_ensembl_genomes,
    configure_logging,
    format_available_species,
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

"""Subtle cyan help and table styling with TTY, NO_COLOR, and TERM=dumb rules (#514)."""

import argparse
import io
import os
from pathlib import Path
import re
import subprocess
import sys

import pytest

from pyensembl import shell
from pyensembl.cli import available, formatting, inspect as inspect_action, listing, parser, verify
from pyensembl.cli.help import (
    ANSI_BOLD,
    ANSI_BOLD_CYAN,
    ANSI_CYAN,
    ANSI_DIM,
    ANSI_RESET,
    ANSI_YELLOW,
    HelpArgumentParser,
    HelpFormatter,
    classify_status,
    color_enabled,
    style_help,
    style_status,
    style_table_header,
)
from .test_cli_actions import FakeRelease
from .test_cli_output import ensembl_files, write

ANSI_ESCAPE_RE = re.compile(r"\033\[[0-9;]*m")
ANSI_CODE_RE = re.compile(r"\033\[([0-9;]*)m")
REPOSITORY = Path(shell.__file__).resolve().parents[1]


class TTYBuffer(io.StringIO):
    """In-memory stream that reports ``isatty() == True``."""

    def isatty(self):
        return True


def strip_ansi(text):
    return ANSI_ESCAPE_RE.sub("", text)


def sgr_codes(text):
    return set(ANSI_CODE_RE.findall(text))


def test_color_enabled_respects_tty_no_color_and_term_dumb(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    tty = TTYBuffer()
    plain = io.StringIO()
    assert color_enabled(tty) is True
    assert color_enabled(plain) is False
    assert color_enabled(None) is False

    monkeypatch.setenv("NO_COLOR", "1")
    assert color_enabled(tty) is False

    monkeypatch.setenv("NO_COLOR", "")
    assert color_enabled(tty) is True

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert color_enabled(tty) is False


@pytest.mark.parametrize("columns", [80, 240])
def test_pyensembl_parser_format_help_stays_plain_and_lists_every_flag(monkeypatch, columns):
    monkeypatch.setenv("COLUMNS", str(columns))
    cli_parser = parser.build_parser()
    assert isinstance(cli_parser, HelpArgumentParser)
    assert cli_parser.formatter_class is HelpFormatter

    text = cli_parser.format_help()
    assert "\033[" not in text
    assert "``" not in text
    assert text.startswith("usage: pyensembl ACTION [options]")
    assert "actions:" in text
    assert "examples:" in text
    for action in cli_parser._actions:
        for flag in action.option_strings:
            assert flag in text


def test_help_argument_parser_styles_tty_and_stays_plain_when_piped(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    cli_parser = parser.build_parser()
    plain_stream = io.StringIO()
    tty_stream = TTYBuffer()

    cli_parser.print_help(file=plain_stream)
    cli_parser.print_help(file=tty_stream)

    plain = plain_stream.getvalue()
    styled = tty_stream.getvalue()

    assert plain == cli_parser.format_help()
    assert "\033[" not in plain
    assert ANSI_BOLD + "usage:" + ANSI_RESET in styled
    assert ANSI_BOLD_CYAN + "Ensembl release options:" + ANSI_RESET in styled
    assert ANSI_BOLD_CYAN + "actions:" + ANSI_RESET in styled
    assert ANSI_BOLD_CYAN + "examples:" + ANSI_RESET in styled
    assert ANSI_CYAN + "--release" + ANSI_RESET in styled
    assert ANSI_CYAN + "-v" + ANSI_RESET in styled
    # One cyan hue + bold only; never bright colors (90-97), green, red, or backgrounds.
    assert sgr_codes(styled) <= {"0", "1", "36", "1;36"}
    assert strip_ansi(styled) == plain


@pytest.mark.parametrize("env", [{"NO_COLOR": "1"}, {"TERM": "dumb"}])
def test_style_help_disables_color_for_no_color_and_term_dumb(monkeypatch, env):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    cli_parser = parser.build_parser()
    tty_stream = TTYBuffer()
    cli_parser.print_help(file=tty_stream)
    assert "\033[" not in tty_stream.getvalue()
    assert tty_stream.getvalue() == cli_parser.format_help()
    assert style_help(cli_parser.format_help(), TTYBuffer()) == cli_parser.format_help()


def test_help_argument_parser_passes_color_false_on_python_3_14_plus(monkeypatch):
    captured_kwargs = {}
    original_init = argparse.ArgumentParser.__init__

    def recording_init(self, *args, **kwargs):
        captured_kwargs.update(kwargs)
        kwargs.pop("color", None)
        return original_init(self, *args, **kwargs)

    monkeypatch.setattr(argparse.ArgumentParser, "__init__", recording_init)
    monkeypatch.setattr(sys, "version_info", (3, 14, 0, "final", 0))
    HelpArgumentParser(prog="pyensembl-test")
    assert captured_kwargs["color"] is False
    assert captured_kwargs["formatter_class"] is HelpFormatter


def test_help_formatter_strips_rst_backticks():
    test_parser = HelpArgumentParser(prog="demo")
    test_parser.add_argument("--flag", help="Pass ``--flag`` to enable.")
    text = test_parser.format_help()
    assert "``" not in text
    assert "--flag" in text


def test_subprocess_help_when_piped_has_no_ansi_escapes(tmp_path):
    command = Path(sys.executable).parent / ("pyensembl.exe" if os.name == "nt" else "pyensembl")
    existing_pythonpath = os.environ.get("PYTHONPATH")
    pythonpath = (
        os.pathsep.join([str(REPOSITORY), existing_pythonpath])
        if existing_pythonpath
        else str(REPOSITORY)
    )
    environment = dict(
        os.environ,
        PYENSEMBL_CACHE_DIR=str(tmp_path),
        PYTHONPATH=pythonpath,
        TERM="xterm-256color",
    )
    environment.pop("NO_COLOR", None)
    result = subprocess.run(
        [str(command), "--help"],
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("usage: pyensembl ACTION [options]")
    assert "\033[" not in result.stdout


def test_style_table_header_rules(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    tty = TTYBuffer()
    plain = io.StringIO()

    assert style_table_header("", tty) == ""
    assert style_table_header("Col1  Col2\nval1  val2\n", plain) == "Col1  Col2\nval1  val2\n"
    assert (
        style_table_header("Col1  Col2\nval1  val2\n", tty)
        == ANSI_BOLD_CYAN + "Col1  Col2" + ANSI_RESET + "\nval1  val2\n"
    )
    # Single line without trailing newline does not gain a stray newline.
    assert (
        style_table_header("Col1  Col2", tty)
        == ANSI_BOLD_CYAN + "Col1  Col2" + ANSI_RESET
    )

    monkeypatch.setenv("NO_COLOR", "1")
    assert style_table_header("Col1  Col2\nval1  val2\n", tty) == "Col1  Col2\nval1  val2\n"


def test_classify_and_style_status_rules(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    tty = TTYBuffer()
    plain = io.StringIO()

    for absent in ("-", "?", "not checked", "unrecorded", "not configured", "unsupported"):
        assert classify_status(absent) == "dim"
        assert style_status(absent, tty) == ANSI_DIM + absent + ANSI_RESET
        assert style_status(absent, plain) == absent

    for warning in (
        "not indexed",
        "incomplete",
        "missing",
        "changed",
        "mismatch",
        "corrupt",
        "invalid reference",
        "invalid-receipt",
        "local dna.fa, missing",
        "toplevel, needs index",
        "toplevel, invalid index",
    ):
        assert classify_status(warning) == "warning"
        assert style_status(warning, tty) == ANSI_YELLOW + warning + ANSI_RESET
        assert style_status(warning, plain) == warning

    # Normal values including "OK" stay plain (no green, no red).
    for normal in ("OK", "indexed", "verified", "available", "unchanged", "toplevel, indexed"):
        assert classify_status(normal) == "plain"
        assert style_status(normal, tty) == normal
        assert "\033[" not in style_status(normal, tty)

    # Explicit level overrides and validation.
    assert style_status("custom", tty, level="dim") == ANSI_DIM + "custom" + ANSI_RESET
    assert style_status("custom", tty, level="warning") == ANSI_YELLOW + "custom" + ANSI_RESET
    assert style_status("missing", tty, level="plain") == "missing"
    with pytest.raises(ValueError, match="Unknown status level"):
        style_status("custom", tty, level="green")

    monkeypatch.setenv("NO_COLOR", "1")
    assert style_status("missing", tty) == "missing"
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert style_status("-", tty) == "-"


def test_format_table_bold_cyan_header_and_status_alignment(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    header = ("Release", "Annotation", "Reference DNA")
    rows = [
        ("83", "indexed", "-"),
        ("81", "not indexed", "local dna.fa, missing"),
        ("82", "incomplete", "toplevel, indexed"),
    ]
    plain = formatting.format_table(
        header, rows, use_color=False, status_columns=("Annotation", "Reference DNA")
    )
    styled = formatting.format_table(
        header, rows, use_color=True, status_columns=("Annotation", "Reference DNA")
    )

    assert "\033[" not in plain
    assert strip_ansi(styled) == plain
    lines = styled.splitlines()
    assert lines[0].startswith(ANSI_BOLD_CYAN) and lines[0].endswith(ANSI_RESET)
    assert ANSI_DIM + "-" + ANSI_RESET in lines[1]
    assert "indexed" in lines[1] and ANSI_YELLOW not in lines[1]
    assert ANSI_YELLOW + "not indexed" + ANSI_RESET in lines[2]
    assert ANSI_YELLOW + "local dna.fa, missing" + ANSI_RESET in lines[2]
    assert ANSI_YELLOW + "incomplete" + ANSI_RESET in lines[3]
    assert "toplevel, indexed" in lines[3]
    # Only bold cyan header, dim absent, yellow warning, and reset.
    assert sgr_codes(styled) <= {"0", "1;36", "2", "33"}

    # Auto-detection via stream respects NO_COLOR and TERM=dumb even on a TTY.
    tty = TTYBuffer()
    auto_tty = formatting.format_table(
        header, rows, status_columns=("Annotation", "Reference DNA"), stream=tty
    )
    assert auto_tty == styled

    monkeypatch.setenv("NO_COLOR", "1")
    assert (
        formatting.format_table(
            header, rows, status_columns=("Annotation", "Reference DNA"), stream=tty
        )
        == plain
    )

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert (
        formatting.format_table(
            header, rows, status_columns=("Annotation", "Reference DNA"), stream=tty
        )
        == plain
    )


def test_list_available_verify_and_inspect_tables_use_shared_palette(tmp_path, monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))

    ensembl_files(83)
    ensembl_files(81, indexes=False)
    write(
        tmp_path / "pyensembl" / "GRCm38" / "dna_only" / "genome_fasta.json",
        '{"source": "/data/dna.fa"}',
    )

    list_plain = listing.format_installed_genomes(use_color=False)
    list_styled = listing.format_installed_genomes(use_color=True)
    assert strip_ansi(list_styled) == list_plain
    assert list_styled.splitlines()[0].startswith(ANSI_BOLD_CYAN)
    assert ANSI_YELLOW + "not indexed" + ANSI_RESET in list_styled
    assert ANSI_YELLOW + "local dna.fa, missing" + ANSI_RESET in list_styled
    assert ANSI_DIM + "-" + ANSI_RESET in list_styled
    assert sgr_codes(list_styled) <= {"0", "1;36", "2", "33"}

    avail_plain = available.format_available_species(use_color=False, dated_releases={})
    avail_styled = available.format_available_species(use_color=True, dated_releases={})
    assert strip_ansi(avail_styled) == avail_plain
    assert avail_styled.splitlines()[0].startswith(ANSI_BOLD_CYAN)
    assert ANSI_BOLD_CYAN + "── Vertebrates " in avail_styled
    assert ANSI_DIM + "homo_sapiens" + ANSI_RESET in avail_styled
    assert ANSI_DIM + "?" + ANSI_RESET in avail_styled
    assert sgr_codes(avail_styled) <= {"0", "1;36", "2"}

    # Verify auto-detection respects NO_COLOR even when sys.stdout is a TTY.
    monkeypatch.setattr(sys, "stdout", TTYBuffer())
    monkeypatch.setenv("NO_COLOR", "1")
    assert "\033[" not in available.format_available_species()
    assert "\033[" not in listing.format_installed_genomes()

    monkeypatch.delenv("NO_COLOR", raising=False)
    verify_plain, ok_plain = verify.run_verify(
        [FakeRelease(["verified", "unrecorded", "mismatch"])], use_color=False
    )
    verify_styled, ok_styled = verify.run_verify(
        [FakeRelease(["verified", "unrecorded", "mismatch"])], use_color=True
    )
    assert ok_plain is False and ok_styled is False
    assert strip_ansi(verify_styled) == verify_plain
    assert verify_styled.splitlines()[0].startswith(ANSI_BOLD_CYAN)
    assert ANSI_DIM + "unrecorded" + ANSI_RESET in verify_styled
    assert ANSI_DIM + "not checked" + ANSI_RESET in verify_styled
    assert ANSI_YELLOW + "mismatch" + ANSI_RESET in verify_styled
    assert sgr_codes(verify_styled) <= {"0", "1;36", "2", "33"}

    from pyensembl import Genome

    missing_gtf = tmp_path / "absent.gtf"
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="test_inspect",
        annotation_version="1",
        gtf_path_or_url=str(missing_gtf),
        cache_directory_path=str(tmp_path / "inspect_cache"),
    )
    inspect_plain = inspect_action.inspect_genomes([genome], use_color=False)
    inspect_styled = inspect_action.inspect_genomes([genome], use_color=True)
    assert strip_ansi(inspect_styled) == inspect_plain
    assert ANSI_BOLD_CYAN + "File" in inspect_styled
    assert ANSI_YELLOW + "missing" + ANSI_RESET in inspect_styled
    assert ANSI_DIM + "-" + ANSI_RESET in inspect_styled
    assert sgr_codes(inspect_styled) <= {"0", "1;36", "2", "33"}


def test_cli_help_reexports_datacache_cli_help_symbols():
    import datacache.cli_help as shared_cli_help
    from pyensembl.cli import help as pyensembl_cli_help

    assert set(shared_cli_help.__all__) <= set(pyensembl_cli_help.__all__)
    for name in shared_cli_help.__all__:
        assert getattr(pyensembl_cli_help, name) is getattr(shared_cli_help, name)


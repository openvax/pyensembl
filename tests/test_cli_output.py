"""What `pyensembl list` and `install` show a person at a terminal."""

import logging
from pathlib import Path

import pytest

from pyensembl import EnsemblRelease
from pyensembl import shell
from .test_genome_fasta import list_rows, run_cli

DATA = Path(__file__).parent / "data"


def touch(paths):
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text("")


def test_list_says_when_nothing_is_installed(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    run_cli(monkeypatch, "list")
    assert capsys.readouterr().out.strip() == (
        "No genomes installed in %s" % (tmp_path / "pyensembl")
    )


def test_list_is_a_table_with_index_status_and_custom_genomes(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    touch(EnsemblRelease(81).required_local_files())  # Downloaded only.
    touch(EnsemblRelease(82)._annotation_index_paths())  # e.g. a custom mirror.
    touch([tmp_path / "pyensembl" / "GRCm38" / "cli_test81" / "annotation.db"])
    rows = list_rows(monkeypatch, capsys)
    assert rows["81"] == {
        "Species": "human",
        "Assembly": "GRCh38",
        "Release": "81",
        "Annotation": "not indexed",
        "Reference DNA": "-",
        "Location": str(tmp_path / "pyensembl" / "GRCh38" / "ensembl81"),
    }
    assert rows["82"]["Annotation"] == "indexed"
    assert rows["cli_test81"]["Species"] == "custom"
    assert rows["cli_test81"]["Assembly"] == "GRCm38"
    assert rows["cli_test81"]["Annotation"] == "indexed"
    run_cli(monkeypatch, "list")
    lines = capsys.readouterr().out.splitlines()
    column = lines[0].index("Annotation")
    assert all(line[column - 2:column] == "  " for line in lines)
    assert "\x1b[" not in "".join(lines)  # No styling when piped.


@pytest.fixture
def cli_logging():
    loggers = [logging.getLogger(name) for name in ("pyensembl", "datacache")]
    saved = [(each, each.level, each.handlers[:]) for each in loggers]
    yield
    for each, level, handlers in saved:
        each.setLevel(level)
        each.handlers[:] = handlers
    shell._cli_handler = None


def test_default_logging_is_concise_and_on_stderr(cli_logging, capsys):
    shell.configure_logging()
    progress = logging.getLogger("pyensembl.progress_test")
    progress.info("Installing human GRCh38 release 81")
    progress.debug("hidden detail")
    logging.getLogger("datacache.database").info('Running sqlite query: "CREATE TABLE"')
    logging.getLogger("datacache.download").warning("slow server")
    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.splitlines() == [
        "Installing human GRCh38 release 81",
        "warning: slow server",
    ]
    shell.configure_logging(verbose=True)  # Replaces, not duplicates, output.
    logging.getLogger("datacache.database").info("Running sqlite query")
    lines = capsys.readouterr().err.splitlines()
    assert len(lines) == 1
    assert lines[0].endswith("datacache.database INFO: Running sqlite query")


def test_install_reports_each_genome_once(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    prefix = str(DATA / "mouse.ensembl.81.partial.")
    arguments = [
        "install", "--reference-name", "GRCm38", "--annotation-name", "cli_test",
        "--annotation-version", "81",
        "--gtf", prefix + "ENSMUSG00000017167.gtf",
        "--transcript-fasta", prefix + "ENSMUSG00000017167.fa",
        "--protein-fasta", prefix + "ENSMUSG00000017167.pep",
    ]
    caplog.set_level(logging.INFO, logger="pyensembl")
    run_cli(monkeypatch, *arguments)
    shell_messages = [r.getMessage() for r in caplog.records if r.name == "pyensembl.shell"]
    assert shell_messages == ["Installing GRCm38 cli_test 81"]
    caplog.clear()
    run_cli(monkeypatch, *arguments)
    assert [r.getMessage() for r in caplog.records if r.levelno >= logging.INFO] == [
        "GRCm38 cli_test 81 is already installed"
    ]

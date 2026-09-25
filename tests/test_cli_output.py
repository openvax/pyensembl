"""What `pyensembl list` and `install` show a person at a terminal."""

from contextlib import closing
import logging
import os
from pathlib import Path
import sqlite3

import pytest

from pyensembl import EnsemblRelease
from pyensembl import shell
from pyensembl.database import DATABASE_SCHEMA_VERSION
from .test_genome_fasta import list_rows, run_cli

DATA = Path(__file__).parent / "data"


def complete_database(path):
    """A database that datacache finished building: its version is set last."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(path)) as connection:
        connection.execute('CREATE TABLE "_datacache_metadata" ("version" INT)')
        connection.execute(
            'INSERT INTO "_datacache_metadata" VALUES (%d)' % DATABASE_SCHEMA_VERSION
        )
        connection.commit()


def write(path, text="data"):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)


def ensembl_files(release, sources=True, indexes=True, database=complete_database):
    """Files as a real install names them, derived from Ensembl's URLs."""
    genome = EnsemblRelease(release)
    directory = Path(genome.download_cache.cache_directory_path)
    names = [os.path.basename(url) for url in genome.transcript_fasta_urls]
    names += [os.path.basename(url) for url in genome.protein_fasta_urls]
    gtf = os.path.basename(genome.gtf_url)
    if sources:
        for name in [gtf] + names:
            write(directory / name)
    if indexes:
        database(directory / (gtf[: -len(".gz")] + ".db"))
        for name in names:
            write(directory / (name + ".pickle"))
    return directory


def test_list_says_when_nothing_is_installed(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    run_cli(monkeypatch, "list")
    assert capsys.readouterr().out.strip() == (
        "No genomes installed in %s" % shell._display_path(tmp_path / "pyensembl")
    )


def test_list_is_a_table_with_install_status(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    root = tmp_path / "pyensembl"
    location = ensembl_files(83)
    ensembl_files(81, indexes=False)  # Downloaded, never indexed.
    ensembl_files(82, sources=False)  # Indexes only, e.g. from a custom mirror.
    ensembl_files(84, database=write)  # Interrupted database build.
    write(root / "GRCh38" / "ensembl90" / ".DS_Store")  # Nothing installed.
    complete_database(root / "GRCh38" / "ensembl200" / "future.gtf.db")
    complete_database(root / "GRCm38" / "cli_test81" / "annotation.db")
    write(root / "GRCm38" / "dna_only" / "genome_fasta.json", '{"source": "/data/dna.fa"}')
    rows = list_rows(monkeypatch, capsys)
    assert rows["83"] == {
        "Species": "human",
        "Assembly": "GRCh38",
        "Release": "83",
        "Annotation": "indexed",
        "Reference DNA": "-",
        "Location": shell._display_path(location),
    }
    assert rows["81"]["Annotation"] == "not indexed"
    assert rows["82"]["Annotation"] == "incomplete"
    assert rows["84"]["Annotation"] == "not indexed"
    assert "90" not in rows
    assert rows["200"]["Species"] == "human"  # e.g. from a newer pyensembl
    assert rows["cli_test81"]["Species"] == "custom"
    assert rows["cli_test81"]["Annotation"] == "indexed"
    assert rows["dna_only"]["Annotation"] == "-"
    assert rows["dna_only"]["Reference DNA"] == "local dna.fa, missing"
    run_cli(monkeypatch, "list")
    lines = capsys.readouterr().out.splitlines()
    column = lines[0].index("Annotation")
    assert all(line[column - 2:column] == "  " for line in lines)
    assert "\x1b[" not in "".join(lines)  # No styling when piped.
    releases = [g.release for g in shell.collect_all_installed_ensembl_releases()]
    assert releases == [81, 82, 83, 84]


def test_list_skips_unreadable_directories(tmp_path, monkeypatch, capsys, caplog):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    ensembl_files(83)
    locked = tmp_path / "pyensembl" / "Locked"
    locked.mkdir()
    locked.chmod(0)
    try:
        rows = list_rows(monkeypatch, capsys)
    finally:
        locked.chmod(0o755)
    assert rows["83"]["Annotation"] == "indexed"
    assert "Cannot read" in caplog.text


@pytest.fixture
def cli_logging():
    loggers = [logging.getLogger(name) for name in ("pyensembl", "datacache")]
    saved = [(each, each.level, each.handlers[:], each.propagate) for each in loggers]
    yield
    for each, level, handlers, propagate in saved:
        each.setLevel(level)
        each.handlers[:] = handlers
        each.propagate = propagate
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


def test_messages_print_once_when_the_host_configured_logging(cli_logging, capsys):
    root = logging.getLogger()
    host = logging.StreamHandler()  # e.g. logging.basicConfig() before shell.run()
    root.addHandler(host)
    try:
        shell.configure_logging()
        logging.getLogger("pyensembl.progress_test").warning("once")
    finally:
        root.removeHandler(host)
    assert capsys.readouterr().err.splitlines() == ["warning: once"]


def install_arguments():
    prefix = str(DATA / "mouse.ensembl.81.partial.")
    return [
        "install", "--reference-name", "GRCm38", "--annotation-name", "cli_test",
        "--annotation-version", "81",
        "--gtf", prefix + "ENSMUSG00000017167.gtf",
        "--transcript-fasta", prefix + "ENSMUSG00000017167.fa",
        "--protein-fasta", prefix + "ENSMUSG00000017167.pep",
    ]


def test_install_reports_each_genome_once(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    arguments = install_arguments()
    caplog.set_level(logging.INFO, logger="pyensembl")
    run_cli(monkeypatch, *arguments)
    shell_messages = [r.getMessage() for r in caplog.records if r.name == "pyensembl.shell"]
    assert shell_messages == ["Installing GRCm38 cli_test 81"]
    caplog.clear()
    run_cli(monkeypatch, *arguments)
    assert [r.getMessage() for r in caplog.records if r.levelno >= logging.INFO] == [
        "GRCm38 cli_test 81 is already installed"
    ]
    # An interrupted database build is not "already installed".
    database = next((tmp_path / "pyensembl" / "GRCm38" / "cli_test81").glob("*.db"))
    database.write_bytes(b"")
    caplog.clear()
    run_cli(monkeypatch, *arguments)
    shell_messages = [r.getMessage() for r in caplog.records if r.name == "pyensembl.shell"]
    assert shell_messages == ["Installing GRCm38 cli_test 81"]

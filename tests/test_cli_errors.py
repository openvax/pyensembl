"""Expected CLI mistakes are concise; data and programming failures stay visible."""

import os
from pathlib import Path
import subprocess
import sys

import pytest

from pyensembl import Genome, shell
from pyensembl.annotation_integrity import AnnotationIntegrityError
from pyensembl.download_cache import DownloadCache, MissingGenomeDataFile
from .test_legacy_cache_compatibility import snapshot


def run_cli(arguments, cache):
    command = Path(sys.executable).parent / ("pyensembl.exe" if os.name == "nt" else "pyensembl")
    # Subprocesses must import the same package as their parent, including when
    # these tests run from an extracted source archive beside an editable env.
    environment = dict(os.environ, PYENSEMBL_CACHE_DIR=str(cache),
                       PYTHONPATH=str(Path(shell.__file__).resolve().parents[1]))
    return subprocess.run([str(command), *arguments], env=environment,
                          capture_output=True, text=True, timeout=30)


@pytest.mark.parametrize("option", ["--gtf", "--transcript-fasta", "--protein-fasta",
                                   "--genome-fasta-path", "--custom-mirror", "--server"])
def test_empty_source_is_rejected_before_genome_selection(tmp_path, monkeypatch, option):
    def forbidden(*args, **kwargs):
        pytest.fail("Empty source reached genome selection")

    args = shell.parser.parse_args(["install", option, ""])
    monkeypatch.setattr(shell, "all_combinations_of_ensembl_genomes", forbidden)
    monkeypatch.setattr(shell, "Genome", forbidden)
    with pytest.raises(ValueError, match=option + " must not be empty"):
        shell.collect_selected_genomes(args)
    result = run_cli(["install", option, ""], tmp_path / "cache")
    assert result.returncode == 2, result.stderr
    assert option + " must not be empty" in result.stderr
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("arguments, message", [
    (["--release", "999"], "No genome"),
    (["--species", "martian"], "Species not found"),
    (["--gtf", "x.gtf"], "reference name"),
    (["--reference-name", "GRCh38", "--gtf", "x.gtf"], "annotation source"),
    (["--release", "75", "--reference-name", "GRCh37", "--annotation-name", "custom", "--gtf", "x.gtf"], "specific paths"),
])
def test_argument_errors_keep_argparse_exit_code(tmp_path, arguments, message):
    result = run_cli(["install", *arguments], tmp_path / "cache")
    assert result.returncode == 2, result.stderr
    assert "pyensembl: error:" in result.stderr and message in result.stderr
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("option", ["--gtf", "--transcript-fasta", "--protein-fasta"])
def test_missing_local_annotation_is_one_runtime_error(tmp_path, option):
    missing = tmp_path / "missing source.fa"
    result = run_cli(["install", "--reference-name", "test", "--annotation-name", "custom",
                      option, str(missing)], tmp_path / "cache")
    assert result.returncode == 1
    assert "pyensembl: error: file not found: " + str(missing) in result.stderr
    assert "Traceback" not in result.stderr and "usage:" not in result.stderr
    assert "MissingFile(" not in result.stderr


def test_missing_local_dna_is_actionable_without_annotation_download(tmp_path):
    missing = tmp_path / "missing dna.fa"
    result = run_cli(["install", "--release", "75", "--only-genome-fasta",
                      "--genome-fasta-path", str(missing)], tmp_path / "cache")
    assert result.returncode == 1
    assert "Missing local genome FASTA" in result.stderr and str(missing) in result.stderr
    assert "Traceback" not in result.stderr


def test_missing_remote_cache_exception_preserves_library_compatibility(tmp_path, monkeypatch, capsys):
    cache = DownloadCache("reference", "annotation", cache_directory_path=str(tmp_path / "absent"))
    with pytest.raises(MissingGenomeDataFile) as error:
        cache.local_path_or_install_error("gtf", "https://example.invalid/genes.gtf")
    assert isinstance(error.value, ValueError)
    assert not (tmp_path / "absent").exists()
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "install", "--release", "75"])

    def missing(*args, **kwargs):
        raise error.value

    monkeypatch.setattr(shell, "_install", missing)
    with pytest.raises(SystemExit) as exit_status:
        shell.run()
    assert exit_status.value.code == 1
    assert "pyensembl: error: Missing genome data file" in capsys.readouterr().err


def test_transcript_only_install_preserves_spaces_and_reopens_read_only(tmp_path):
    source = tmp_path / " transcripts with spaces.fa "
    source.write_text(">t\nACGT\n")
    cache = tmp_path / "cache"
    args = ["install", "--reference-name", "test", "--annotation-name", "custom",
            "--transcript-fasta", str(source)]
    first = run_cli(args, cache)
    assert first.returncode == 0, first.stderr
    before = snapshot(cache)
    second = run_cli(args, cache)
    assert second.returncode == 0, second.stderr
    assert snapshot(cache) == before


def test_malformed_gtf_keeps_traceback_and_source_context(tmp_path):
    source = tmp_path / "malformed.gtf"
    source.write_text("not\ta\tvalid\tGTF\n")
    result = run_cli(["install", "--reference-name", "test", "--annotation-name", "custom",
                      "--gtf", str(source)], tmp_path / "cache")
    assert result.returncode == 1
    assert "Traceback (most recent call last)" in result.stderr
    assert "ParsingError" in result.stderr
    assert "malformed.gtf" in result.stderr


@pytest.mark.parametrize("failure", [ValueError("internal indexing failure"),
                                     AnnotationIntegrityError("MD5 mismatch")])
def test_unexpected_value_errors_are_not_converted_to_cli_errors(monkeypatch, failure):
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "install", "--release", "75"])

    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(shell, "_install", fail)
    with pytest.raises(type(failure), match=str(failure)):
        shell.run()


def test_deleting_an_uninstalled_concrete_release_remains_a_noop(tmp_path):
    cache = tmp_path / "cache"
    result = run_cli(["delete-index-files", "--release", "75"], cache)
    assert result.returncode == 0, result.stderr
    assert "Nothing to delete" in result.stdout
    assert not cache.exists()


def test_corrupt_database_without_overwrite_keeps_traceback(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    source = tmp_path / "genes.gtf"
    source.write_text('1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g"; gene_name "G";\n')
    with Genome(reference_name="test", annotation_name="custom", gtf_path_or_url=str(source)) as genome:
        genome.index()
        database = Path(genome.db.local_db_path)
    database.write_bytes(b"not a sqlite database")
    before = snapshot(tmp_path / "cache")
    result = run_cli(["install", "--reference-name", "test", "--annotation-name", "custom",
                      "--gtf", str(source)], tmp_path / "cache")
    assert result.returncode == 1
    assert "Traceback" in result.stderr and "DatabaseError" in result.stderr
    assert snapshot(tmp_path / "cache") == before

"""Progress display and import cost, via datacache's optional features."""

import gzip
import subprocess
import sys

import datacache

from pyensembl import EnsemblRelease, Genome
from pyensembl import genome_fasta, shell
from .test_cli_output import install_arguments
from .test_genome_fasta import DNA, run_cli, serve_downloads

GTF = (
    '1\ttest\tgene\t1\t12\t.\t+\t.\tgene_id "g"; gene_name "g";\n'
    '1\ttest\ttranscript\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    '1\ttest\texon\t1\t12\t.\t+\t.\tgene_id "g"; transcript_id "t"; exon_id "e";\n'
)


def record_progress(monkeypatch):
    """Record show_progress passed to datacache; run without tqdm installed."""
    seen = []
    for name in ("fetch_file", "db_from_dataframes_with_absolute_path"):
        real = getattr(datacache, name)

        def spy(*args, _real=real, _name=name, show_progress=False, **kwargs):
            seen.append((_name, show_progress))
            return _real(*args, **kwargs)

        monkeypatch.setattr(datacache, name, spy)
    return seen


def test_importing_pyensembl_does_not_load_pandas_or_polars():
    probe = (
        "import sys, pyensembl, pyensembl.shell; "
        "print(sorted(m for m in ('pandas', 'polars', 'pyarrow', 'gtfparse', 'requests') "
        "if m in sys.modules))"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "[]"


def test_downloads_and_database_builds_forward_show_progress(tmp_path, monkeypatch):
    source = tmp_path / "annotation.gtf.gz"
    source.write_bytes(gzip.compress(GTF.encode()))
    seen = record_progress(monkeypatch)
    genome = Genome(
        "synthetic", "progress",
        gtf_path_or_url=source.as_uri(),
        cache_directory_path=str(tmp_path / "cache"),
    )
    genome.download(show_progress=True)
    genome.index(show_progress=True)
    assert seen == [
        ("fetch_file", True),
        ("db_from_dataframes_with_absolute_path", True),
    ]
    assert genome.gene_names() == ["g"]
    genome.close()


def test_reference_dna_downloads_forward_show_progress(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        "pyensembl.genome_fasta_cache._remote_identity", lambda source: {"source": source}
    )
    requested = []
    serve_downloads(monkeypatch, tmp_path, lambda url: gzip.compress(DNA))
    served = genome_fasta.fetch_file  # The test server installed above.

    def spy(url, show_progress=False, **kwargs):
        requested.append(show_progress)
        return served(url, **kwargs)

    monkeypatch.setattr("pyensembl.genome_fasta.fetch_file", spy)
    genome = EnsemblRelease(81, genome_fasta=True)
    genome.download_genome_fasta(show_progress=True)
    assert requested == [True]
    assert genome.sequence("MT", 1, 4) == "GCTA"
    genome.close()


def test_cli_shows_progress_only_when_someone_can_see_it(tmp_path, monkeypatch, capsys):
    # capsys replaces stderr with a pipe: no terminal, no progress bars.
    assert shell._progress_available() is False
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(shell, "_progress_available", lambda: True)
    seen = record_progress(monkeypatch)
    run_cli(monkeypatch, *install_arguments())
    assert seen == [("db_from_dataframes_with_absolute_path", True)]

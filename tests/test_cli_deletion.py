"""Deletion commands must select their target explicitly and work offline."""

from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from pyensembl import EnsemblRelease, Genome
from pyensembl import shell
from pyensembl.ensembl_versions import MAX_ENSEMBL_RELEASE
from pyensembl.shell import collect_selected_genomes, parser


@pytest.fixture
def run_cli(monkeypatch, capsys):
    # Exercise the real parser and entrypoint without repeatedly importing the
    # full scientific stack in child interpreters. Logging setup has its own
    # tests; avoid changing process-global handlers during these calls.
    monkeypatch.setattr(shell, "configure_logging", lambda: None)

    def invoke(cache_root, *arguments):
        monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache_root))
        monkeypatch.setattr(sys, "argv", ["pyensembl", *arguments])
        returncode = 0
        try:
            shell.run()
        except SystemExit as error:
            returncode = error.code
        output = capsys.readouterr()
        return SimpleNamespace(
            returncode=returncode, stdout=output.out, stderr=output.err,
        )

    return invoke


@pytest.mark.parametrize("action", ["delete-all-files", "delete-index-files"])
@pytest.mark.parametrize("mirror", [False, True])
def test_deletion_requires_explicit_release(tmp_path, monkeypatch, action, mirror, run_cli):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    genome = EnsemblRelease(MAX_ENSEMBL_RELEASE)
    directory = Path(genome.download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    source = directory / "keep.gtf"
    source.write_bytes(b"preserve this data")
    arguments = [action]
    if mirror:
        arguments += ["--custom-mirror", "https://example.invalid/ensembl"]

    result = run_cli(tmp_path, *arguments)

    assert result.returncode == 2
    assert "requires an explicit --release" in result.stderr
    assert "Traceback" not in result.stderr
    assert source.read_bytes() == b"preserve this data"


def test_install_retains_default_release(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    genomes = collect_selected_genomes(parser.parse_args(["install"]))
    assert len(genomes) == 1
    assert genomes[0].release == MAX_ENSEMBL_RELEASE


def test_delete_all_reports_size_and_preserves_other_releases(tmp_path, monkeypatch, run_cli):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    selected = Path(EnsemblRelease(81).download_cache.cache_directory_path)
    untouched = Path(EnsemblRelease(82).download_cache.cache_directory_path)
    (selected / "nested").mkdir(parents=True)
    (selected / "source.gtf").write_bytes(b"source")
    (selected / "nested" / "index.db").write_bytes(b"index")
    untouched.mkdir(parents=True)
    (untouched / "keep.gtf").write_bytes(b"keep")

    result = run_cli(tmp_path, "delete-all-files", "--release", "81")

    assert result.returncode == 0, result.stderr
    assert not selected.exists()
    assert (untouched / "keep.gtf").read_bytes() == b"keep"
    assert "Deleted %s (11 bytes)" % selected in result.stdout


@pytest.mark.parametrize("action", ["delete-all-files", "delete-index-files"])
def test_empty_cache_reports_every_selected_genome(tmp_path, action, run_cli):
    result = run_cli(
        tmp_path, action, "--release", "81", "82", "--species", "human", "mouse"
    )

    assert result.returncode == 0, result.stderr
    for species, assembly in [("human", "GRCh38"), ("mouse", "GRCm38")]:
        for release in (81, 82):
            assert "Nothing to delete for %s %s release %d" % (
                species, assembly, release,
            ) in result.stdout
    assert "Traceback" not in result.stderr
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("sources_present", [False, True])
def test_delete_indexes_without_loading_sources(tmp_path, monkeypatch, sources_present, run_cli):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    genome = EnsemblRelease(81)
    directory = Path(genome.download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    sources = [
        Path(genome.download_cache.cached_path(url))
        for url in [genome.gtf_url] + genome.transcript_fasta_urls + genome.protein_fasta_urls
    ]
    indexes = [sources[0].with_suffix(".db")]
    indexes += [Path(str(path) + ".pickle") for path in sources[1:]]
    for path in indexes:
        path.write_bytes(b"index")
    if sources_present:
        for path in sources:
            # Deliberately not valid GTF/FASTA: deletion must not parse them.
            path.write_bytes(b"source")
    unrelated = directory / "unrelated.db"
    unrelated.write_bytes(b"keep")

    result = run_cli(tmp_path, "delete-index-files", "--release", "81")

    assert result.returncode == 0, result.stderr
    for path in indexes:
        assert not path.exists()
        assert "Deleted %s (5 bytes)" % path in result.stdout
    assert unrelated.read_bytes() == b"keep"
    for path in sources:
        if sources_present:
            assert path.read_bytes() == b"source"
        else:
            assert not path.exists()
    repeated = run_cli(tmp_path, "delete-index-files", "--release", "81")
    assert repeated.returncode == 0, repeated.stderr
    assert "Nothing to delete for human GRCh38 release 81" in repeated.stdout


@pytest.mark.parametrize("action", ["delete-all-files", "delete-index-files"])
def test_explicit_custom_genome_deletion(tmp_path, monkeypatch, action, run_cli):
    cache_root = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache_root))
    source = tmp_path / "custom.gtf"
    source.write_bytes(b"source")
    genome = Genome("TestAssembly", "custom", gtf_path_or_url=str(source))
    directory = Path(genome.download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    index = directory / "custom.db"
    index.write_bytes(b"index")
    arguments = [
        action, "--reference-name", "TestAssembly", "--annotation-name", "custom",
        "--gtf", str(source),
    ]

    result = run_cli(cache_root, *arguments)

    assert result.returncode == 0, result.stderr
    assert not index.exists()
    assert source.read_bytes() == b"source"
    removed = directory if action == "delete-all-files" else index
    assert "Deleted %s (5 bytes)" % removed in result.stdout
    repeated = run_cli(cache_root, *arguments)
    assert repeated.returncode == 0, repeated.stderr
    assert "Nothing to delete for TestAssembly custom" in repeated.stdout


@pytest.mark.parametrize("remote", [False, True])
@pytest.mark.parametrize("copy_local", [False, True])
@pytest.mark.parametrize("decompress", [False, True])
def test_orphan_index_paths_follow_cache_policy(
    tmp_path, monkeypatch, remote, copy_local, decompress,
):
    prefix = "https://example.invalid/" if remote else str(tmp_path / "sources") + "/"
    cache = tmp_path / "cache"
    cache.mkdir()
    genome = Genome(
        "TestAssembly", "custom",
        gtf_path_or_url=prefix + "test.gtf.gz",
        transcript_fasta_paths_or_urls=[prefix + "test.fa.gz"],
        cache_directory_path=str(cache),
        copy_local_files_to_cache=copy_local,
        decompress_on_download=decompress,
    )

    def no_source_access(*args, **kwargs):
        pytest.fail("Index deletion attempted to download or copy a source")

    monkeypatch.setattr(genome.download_cache, "download_or_copy_if_necessary", no_source_access)
    if decompress and (remote or copy_local):
        indexes = [cache / "test.db", cache / "test.fa.pickle"]
    else:
        indexes = [cache / "test.gtf.db", cache / "test.fa.gz.pickle"]
    for index in indexes:
        index.write_bytes(b"index")

    deleted = genome.delete_index_files()

    assert set(deleted) == {(str(path), 5) for path in indexes}
    assert not any(path.exists() for path in indexes)
    assert genome.delete_index_files() == []
    assert not (tmp_path / "sources").exists()

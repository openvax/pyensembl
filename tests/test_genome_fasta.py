"""Reference DNA tests use real indexed FASTA files and synthetic known bases."""

import gzip
import json
import os
from pathlib import Path
import pickle
import shlex
import sys
from urllib.parse import urlsplit
from uuid import uuid4

import datacache
import pytest

from pyensembl import EnsemblRelease, Genome, MissingGenomeFastaError
from pyensembl.ensembl_url_templates import make_genome_fasta_url
from pyensembl.genome_fasta import GenomeFasta
from pyensembl import shell


DNA = b">1 chromosome\nACgtNN\nTTaacc\nGGTA\n>MT mitochondrion\nGCTA\n>CHR_PATCH\nNNacGT\n"


def serve_downloads(monkeypatch, directory, payload):
    """Route DNA downloads through real datacache staging and validation.

    payload(url) returns the bytes to serve. Returns the requested URLs.
    """
    calls = []

    def fetch(url, **kwargs):
        calls.append(url)
        served = Path(directory) / "served" / uuid4().hex
        served.mkdir(parents=True)
        path = served / os.path.basename(urlsplit(url).path)
        path.write_bytes(payload(url))
        return datacache.fetch_file(path.as_uri(), **kwargs)

    monkeypatch.setattr("pyensembl.genome_fasta.fetch_file", fetch)
    return calls


def forbid_downloads(monkeypatch, message):
    def unexpected(*args, **kwargs):
        pytest.fail(message)

    monkeypatch.setattr("pyensembl.genome_fasta.fetch_file", unexpected)


@pytest.fixture
def dna_path(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "global_cache"))
    # Most tests isolate transport from shared metadata; dedicated cache tests
    # exercise the actual metadata parser and cross-release reuse.
    monkeypatch.setattr(
        "pyensembl.genome_fasta_cache._remote_identity",
        lambda source: {"source": source},
    )
    path = tmp_path / "input" / "dna.fa"
    path.parent.mkdir()
    path.write_bytes(DNA)
    return path


def custom_genome(tmp_path, source, **kwargs):
    return Genome(
        "synthetic",
        "test",
        genome_fasta_path_or_url=source,
        cache_directory_path=str(tmp_path / "cache"),
        **kwargs,
    )


def test_sequence_coordinates_masking_and_line_boundaries(tmp_path, dna_path):
    with custom_genome(tmp_path, dna_path) as genome:
        assert genome.sequence(1, 1, 1) == "A"
        assert genome.sequence("1", 16, 16) == "A"
        assert genome.sequence("1", 3, 10) == "GTNNTTAA"
        assert genome.sequence("1", 3, 10, mask="raw") == "gtNNTTaa"
        assert genome.sequence("1", 1, 16) == "ACGTNNTTAACCGGTA"
        assert genome.sequence("CHR_PATCH", 3, 6) == "ACGT"
        assert genome.fasta["MT"][1:3].seq == "CT"
        assert genome.genome_fasta_path == str(dna_path)


@pytest.mark.parametrize(
    "start,end",
    [
        (0, 1),
        (-1, 2),
        (2, 1),
        (1, 17),
        (17, 17),
        (True, 2),
        (1, False),
        (1.0, 2),
        (1, "2"),
    ],
)
def test_invalid_intervals_fail_without_truncation(tmp_path, dna_path, start, end):
    with custom_genome(tmp_path, dna_path) as genome:
        with pytest.raises(ValueError):
            genome.sequence("1", start, end)


def test_absent_contig_and_invalid_mask(tmp_path, dna_path):
    with custom_genome(tmp_path, dna_path) as genome:
        with pytest.raises(ValueError, match="Contig 'chr1' is absent"):
            genome.sequence("chr1", 1, 2)
        with pytest.raises(ValueError, match="mask"):
            genome.sequence("1", 1, 2, mask="reverse")


def test_missing_and_default_configuration_never_download(
    tmp_path, dna_path, monkeypatch
):
    forbid_downloads(monkeypatch, "Default or lazy sequence lookup accessed the network")
    for genome in (custom_genome(tmp_path, None), EnsemblRelease(81)):
        assert not genome.requires_genome_fasta
        assert genome.genome_fasta_path is None
        assert genome.fasta is None
        with pytest.raises(MissingGenomeFastaError, match="No genome FASTA"):
            genome.sequence("1", 1, 2)
    remote = EnsemblRelease(81, download_genome_fasta=True)
    assert remote.requires_genome_fasta
    assert remote.genome_fasta_path is None
    with pytest.raises(MissingGenomeFastaError, match="download_genome_fasta"):
        remote.sequence("1", 1, 2)
    missing = custom_genome(tmp_path, dna_path.parent / "missing.fa")
    with pytest.raises(MissingGenomeFastaError, match="Missing local"):
        missing.sequence("1", 1, 2)
    # Duck-typed consumers (Varcode) fall back when .fasta is None.
    for genome in (remote, missing):
        assert genome.fasta is None
        assert getattr(genome, "fasta", "absent") is None


def test_missing_dna_error_gives_a_runnable_dna_only_command(dna_path):
    release = EnsemblRelease(
        81,
        download_genome_fasta=True,
        genome_fasta_type="primary_assembly",
        genome_fasta_mask="soft",
    )
    with pytest.raises(MissingGenomeFastaError) as error:
        release.sequence("1", 1, 2)
    command = release.genome_fasta_install_string()
    assert command in str(error.value)
    assert command == (
        "pyensembl install --release 81 --species homo_sapiens --only-genome-fasta "
        "--genome-fasta-type primary_assembly --masked soft"
    )
    args = shell.parser.parse_args(shlex.split(command)[1:])
    assert args.only_genome_fasta
    (selected,) = shell.collect_selected_genomes(args)
    assert selected.genome_fasta_urls == release.genome_fasta_urls
    local = EnsemblRelease(81, genome_fasta_path=dna_path)
    assert local.genome_fasta_install_string().endswith(
        "--only-genome-fasta --genome-fasta-path %s" % shlex.quote(str(dna_path))
    )
    assert "--with-genome-fasta" in release.install_string()


@pytest.mark.parametrize("compressed", [False, True])
def test_local_sources_and_existing_user_indexes_are_untouched(
    tmp_path, dna_path, compressed
):
    if compressed:
        dna_path.write_bytes(gzip.compress(DNA))  # Detect magic, not suffix.
    source_bytes = dna_path.read_bytes()
    user_index = Path(str(dna_path) + ".fai")
    user_index.write_text("user-owned index\n")
    dna_path.chmod(0o444)
    dna_path.parent.chmod(0o555)
    try:
        genome = custom_genome(tmp_path, dna_path)
        genome.download()
        genome.index()
        assert genome.sequence("1", 3, 10) == "GTNNTTAA"
        installed = Path(genome.genome_fasta_path)
        assert genome._genome_fasta.index_path.exists()
        genome.delete_index_files()
        assert not genome._genome_fasta.index_path.exists()
        assert installed.exists()
        assert genome.sequence("MT", 1, 4) == "GCTA"
        genome.close()
        genome.download_cache.delete_cache_directory()
        assert dna_path.read_bytes() == source_bytes
        assert user_index.read_text() == "user-owned index\n"
        assert sorted(x.name for x in dna_path.parent.iterdir()) == [
            "dna.fa",
            "dna.fa.fai",
        ]
    finally:
        dna_path.parent.chmod(0o755)
        dna_path.chmod(0o644)


@pytest.mark.parametrize("compressed", [False, True])
def test_changed_local_source_rebuilds_reader_and_index(tmp_path, dna_path, compressed):
    encode = gzip.compress if compressed else lambda b: b
    dna_path.write_bytes(encode(DNA))
    with custom_genome(tmp_path, dna_path) as genome:
        assert genome.sequence("MT", 1, 4) == "GCTA"
        old_reader = genome.fasta
        old_stat = dna_path.stat()
        dna_path.write_bytes(encode(b">MT\nAACCGG\n"))
        # A preserved/older mtime must not cause a stale FAI to be reused.
        os.utime(dna_path, ns=(old_stat.st_atime_ns, old_stat.st_mtime_ns))
        assert genome.sequence("MT", 1, 6) == "AACCGG"
        assert genome.fasta is not old_reader


def test_close_clear_and_context_manager_reopen(tmp_path, dna_path):
    genome = custom_genome(tmp_path, dna_path)
    reader = genome.fasta
    genome.close()
    assert reader.faidx.file.closed
    assert genome.sequence("MT", 1, 4) == "GCTA"
    # Clearing in-memory caches must not break readers callers still hold,
    # e.g. varcode.Genome stores genome.fasta.
    held = genome.fasta
    genome.clear_cache()
    assert held["MT"][0:4].seq == "GCTA"
    assert genome.fasta is not held
    assert genome.sequence("MT", 1, 4) == "GCTA"
    with genome:
        reader = genome.fasta
    assert reader.faidx.file.closed


def test_same_basename_sources_do_not_reuse_each_others_index(tmp_path, dna_path):
    other = tmp_path / "other" / dna_path.name
    other.parent.mkdir()
    other.write_text(">1\nTT\n")
    with (
        custom_genome(tmp_path, dna_path) as first,
        custom_genome(tmp_path, other) as second,
    ):
        assert first.sequence("1", 1, 2) == "AC"
        assert second.sequence("1", 1, 2) == "TT"
        assert first._genome_fasta.index_path != second._genome_fasta.index_path


def test_remote_gzip_download_is_explicit_reusable_and_independent(
    tmp_path, dna_path, monkeypatch
):
    calls = serve_downloads(monkeypatch, tmp_path, lambda url: gzip.compress(DNA))
    genome = EnsemblRelease(81, download_genome_fasta=True)
    genome.download_genome_fasta()
    genome.index_genome_fasta()
    assert genome.sequence("1", 3, 10) == "GTNNTTAA"
    assert calls == genome.genome_fasta_urls
    assert not Path(genome.download_cache.cached_path(genome.gtf_url)).exists()
    assert not genome.required_local_files_exist()
    genome.close()
    reloaded = EnsemblRelease(81, download_genome_fasta=True)
    reloaded.download_genome_fasta()
    assert reloaded.sequence("MT", 1, 4) == "GCTA"
    assert len(calls) == 1
    reloaded.close()
    different_release = EnsemblRelease(82, download_genome_fasta=True)
    with pytest.raises(MissingGenomeFastaError):
        different_release.sequence("1", 1, 2)


@pytest.mark.parametrize("payload", [b"not FASTA", gzip.compress(DNA)[:-7]])
def test_failed_download_does_not_publish_partial_data(
    tmp_path, dna_path, monkeypatch, payload
):
    serve_downloads(monkeypatch, tmp_path, lambda url: payload)
    genome = custom_genome(tmp_path, "https://example.test/genome.fa.gz")
    with pytest.raises((ValueError, EOFError)):
        genome.download_genome_fasta()
    assert genome.genome_fasta_path is None
    assert not genome._genome_fasta.materialized_path.exists()
    assert not list(genome._genome_fasta.directory.iterdir())


def test_failed_overwrite_preserves_previous_download(tmp_path, dna_path, monkeypatch):
    payload = [gzip.compress(DNA)]
    serve_downloads(monkeypatch, tmp_path, lambda url: payload[0])
    genome = custom_genome(tmp_path, "https://example.test/genome.fa.gz")
    genome.download_genome_fasta()
    assert genome.sequence("MT", 1, 4) == "GCTA"
    payload[0] = gzip.compress(DNA)[:-4]
    with pytest.raises(EOFError):
        genome.download_genome_fasta(overwrite=True)
    assert genome.sequence("MT", 1, 4) == "GCTA"
    genome.close()


def test_published_files_reach_disk_before_rename(tmp_path, dna_path, monkeypatch):
    dna_path.write_bytes(gzip.compress(DNA))
    synced, replaced = set(), []
    real_fsync, real_replace = os.fsync, os.replace

    def fsync(descriptor):
        synced.add(os.fstat(descriptor).st_ino)
        real_fsync(descriptor)

    def replace(source, destination):
        replaced.append((os.stat(source).st_ino, Path(destination).name))
        real_replace(source, destination)

    monkeypatch.setattr(os, "fsync", fsync)
    monkeypatch.setattr(os, "replace", replace)
    with custom_genome(tmp_path, dna_path) as genome:
        genome.index_genome_fasta()
    names = {name for _, name in replaced}
    assert {"sequence.fa", "sequence.fa.fai", "source.json", "index.json"} <= names
    assert all(inode in synced for inode, _ in replaced)


def test_duplicate_contigs_are_rejected(tmp_path, dna_path):
    dna_path.write_text(">1\nAC\n>1\nTT\n")
    genome = custom_genome(tmp_path, dna_path)
    with pytest.raises(ValueError, match="Duplicate"):
        genome.index_genome_fasta()
    assert not genome._genome_fasta.index_path.exists()


def test_serialization_and_cached_releases_preserve_dna_configuration(
    tmp_path, dna_path
):
    genome = custom_genome(tmp_path, dna_path)
    assert genome.sequence("MT", 1, 4) == "GCTA"
    for restored in (
        pickle.loads(pickle.dumps(genome)),
        Genome.from_json(genome.to_json()),
    ):
        assert restored == genome
        assert restored.sequence("MT", 1, 4) == "GCTA"
        restored.close()
    genome.close()
    plain = EnsemblRelease.cached(81)
    local = EnsemblRelease.cached(81, genome_fasta_path=dna_path)
    remote = EnsemblRelease.cached(81, download_genome_fasta=True)
    assert plain is not local and local is not remote
    assert EnsemblRelease.cached(81, genome_fasta_path=str(dna_path)) is local
    assert local.sequence("MT", 1, 4) == "GCTA"
    for genome in (plain, local, remote):
        assert pickle.loads(pickle.dumps(genome)) is genome
        assert EnsemblRelease.from_json(genome.to_json()) is genome
        genome.close()


def test_attached_dna_does_not_change_annotation_equality(tmp_path, dna_path):
    plain = EnsemblRelease(81)
    with_dna = EnsemblRelease(81, download_genome_fasta=True)
    assert plain == with_dna and hash(plain) == hash(with_dna)
    assert plain.to_dict() != with_dna.to_dict()  # Serialization keeps DNA.
    # Gene and Transcript equality compare genomes.
    assert {plain: "annotation"}[with_dna] == "annotation"
    assert custom_genome(tmp_path, dna_path) == custom_genome(tmp_path, None)
    assert EnsemblRelease(82) != plain


@pytest.mark.parametrize(
    "release,species,flavor,mask,expected",
    [
        (
            75,
            "human",
            "toplevel",
            "none",
            "/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz",
        ),
        (
            81,
            "human",
            "primary_assembly",
            "soft",
            "/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz",
        ),
        (
            103,
            "mouse",
            "toplevel",
            "hard",
            "/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.toplevel.fa.gz",
        ),
        (
            58,
            "arabidopsis_thaliana",
            "toplevel",
            "none",
            "/release-58/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
        ),
    ],
)
def test_official_archive_layouts(release, species, flavor, mask, expected):
    url = make_genome_fasta_url(release, species, flavor, mask)
    assert url.endswith(expected)
    if species == "arabidopsis_thaliana":
        assert url.startswith("https://ftp.ensemblgenomes.ebi.ac.uk/")
    else:
        assert url.startswith("https://ftp.ensembl.org/")


def test_local_file_takes_precedence_and_mirrors_remain_flat(tmp_path, dna_path):
    release = EnsemblRelease(81, download_genome_fasta=True, genome_fasta_path=dna_path)
    assert release.genome_fasta_urls == []
    assert release.sequence("MT", 1, 4) == "GCTA"
    release.close()
    args = shell.parser.parse_args(
        [
            "install",
            "--release",
            "81",
            "--with-genome-fasta",
            "--custom-mirror",
            "https://example.test/files",
            "--masked",
            "soft",
        ]
    )
    (genome,) = shell.collect_selected_genomes(args)
    assert genome._genome_fasta_path_or_url == (
        "https://example.test/files/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
    )


def run_cli(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["pyensembl", *args])
    monkeypatch.setattr(shell, "configure_logging", lambda: None)
    shell.run()


def test_cli_only_local_dna_and_inspection_are_offline(
    tmp_path, dna_path, monkeypatch, capsys
):
    def unexpected(*args, **kwargs):
        pytest.fail("DNA-only local install tried to download")

    forbid_downloads(monkeypatch, "DNA-only local install tried to download")
    monkeypatch.setattr(
        "pyensembl.download_cache.DownloadCache._download_if_necessary", unexpected
    )
    run_cli(
        monkeypatch,
        "install",
        "--release",
        "81",
        "--only-genome-fasta",
        "--genome-fasta-path",
        str(dna_path),
    )
    # Listing must discover DNA-only installations and not index/download.
    monkeypatch.setattr(GenomeFasta, "open", unexpected)
    run_cli(monkeypatch, "list", "--check-genome-fasta")
    output = capsys.readouterr().out
    assert "release=81" in output
    assert "Genome FASTA: local, indexed" in output
    assert str(dna_path) in output
    genome = EnsemblRelease(81)  # No DNA constructor flag needed for deletion.
    installed = GenomeFasta.installed_source(genome.download_cache.cache_directory_path)
    assert installed.index_path.exists()
    run_cli(monkeypatch, "delete-index-files", "--release", "81")
    assert not installed.index_path.exists()
    run_cli(monkeypatch, "list", "--check-genome-fasta")
    assert "needs index" in capsys.readouterr().out
    assert not installed.index_path.exists()
    run_cli(monkeypatch, "delete-all-files", "--release", "81")
    assert dna_path.read_bytes() == DNA
    assert not installed.manifest_path.exists()


def test_custom_fasta_only_install_without_gtf(tmp_path, dna_path, monkeypatch):
    run_cli(
        monkeypatch,
        "install",
        "--reference-name",
        "synthetic",
        "--annotation-name",
        "custom",
        "--genome-fasta-path",
        str(dna_path),
    )
    # Also cover the preexisting transcript-only CLI failure (#400).
    args = shell.parser.parse_args(
        [
            "install",
            "--reference-name",
            "synthetic",
            "--annotation-name",
            "custom",
            "--transcript-fasta",
            str(dna_path),
        ]
    )
    (genome,) = shell.collect_selected_genomes(args)
    assert not genome.requires_gtf
    genome.download()
    genome.index()
    assert genome.transcript_sequence("MT") == "GCTA"


def test_local_annotation_coverage_warning_and_intronic_sequence(tmp_path, dna_path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1\t12\t.\t-\t.\tgene_id "g"; gene_name "g";\n'
        '1\ttest\ttranscript\t1\t12\t.\t-\t.\tgene_id "g"; transcript_id "t";\n'
        '1\ttest\texon\t1\t2\t.\t-\t.\tgene_id "g"; transcript_id "t"; exon_id "e1";\n'
        '1\ttest\texon\t11\t12\t.\t-\t.\tgene_id "g"; transcript_id "t"; exon_id "e2";\n'
        'absent\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "missing"; gene_name "missing";\n'
    )
    with custom_genome(tmp_path, dna_path, gtf_path_or_url=str(gtf)) as genome:
        with pytest.warns(UserWarning, match="lacks 1 annotation contigs"):
            genome.index()
        assert genome.transcript_by_id("t").strand == "-"
        assert genome.sequence("1", 3, 10) == "GTNNTTAA"  # Intron, always plus strand.
        assert genome.sequence("1", 13, 16) == "GGTA"  # Outside the annotated gene.


def test_status_detects_missing_stale_and_invalid_indexes(tmp_path, dna_path):
    genome = custom_genome(tmp_path, dna_path)
    genome.index_genome_fasta()
    genome.close()
    source = GenomeFasta.installed_source(genome.download_cache.cache_directory_path)
    assert "indexed" in source.status(check=True)
    source.index_path.write_text("broken index\n")
    assert "invalid index" in source.status(check=True)
    genome.index_genome_fasta(overwrite=True)
    genome.close()
    dna_path.write_bytes(b">1\nGG\n")
    assert "needs index" in source.status(check=True)
    dna_path.unlink()
    assert "missing" in source.status(check=True)
    source.manifest_path.write_text("{invalid json")
    assert (
        GenomeFasta.installed_source(genome.download_cache.cache_directory_path) is None
    )
    source.manifest_path.write_text(json.dumps({"source": None}))
    assert (
        GenomeFasta.installed_source(genome.download_cache.cache_directory_path) is None
    )

"""Reference requests must never act on another assembly (issue #223)."""

import sys

import pytest

from pyensembl import EnsemblRelease, Genome, shell
from pyensembl.ensembl_versions import MAX_ENSEMBL_GENOMES_RELEASE, MAX_ENSEMBL_RELEASE


def select(*arguments):
    return shell.collect_selected_genomes(
        shell.parser.parse_args(["install", *arguments])
    )


@pytest.mark.parametrize("reference", ["GRCh37", "grch37", " GrCh37 "])
@pytest.mark.parametrize("species", [None, "human", "Homo sapiens"])
def test_grch37_selects_latest_matching_release(reference, species):
    arguments = ["--reference-name", reference]
    if species:
        arguments += ["--species", species]
    genome, = select(*arguments)

    assert genome.reference_name == "GRCh37"
    assert genome.release == 75
    assert genome.species.latin_name == "homo_sapiens"
    assert "/release-75/" in genome.gtf_url
    assert "Homo_sapiens.GRCh37.75.gtf.gz" in genome.gtf_url
    assert all(
        ".GRCh37." in url
        for url in genome.transcript_fasta_urls + genome.protein_fasta_urls
    )


@pytest.mark.parametrize(
    "reference, species, release",
    [("GRCm38", "mus_musculus", 102), ("NCBI36", "homo_sapiens", 54),
     ("GRCh38", "homo_sapiens", MAX_ENSEMBL_RELEASE),
     ("TAIR10", "arabidopsis_thaliana", MAX_ENSEMBL_GENOMES_RELEASE)],
)
def test_reference_infers_species_and_release(reference, species, release):
    genome, = select("--reference-name", reference)
    assert genome.reference_name == reference
    assert genome.species.latin_name == species
    assert genome.release == release


def test_matching_explicit_releases_are_preserved():
    genomes = select("--reference-name", "GRCh37", "--release", "55", "74", "75")
    assert [genome.release for genome in genomes] == [55, 74, 75]
    assert all(genome.reference_name == "GRCh37" for genome in genomes)


def test_reference_install_does_not_prefer_older_cached_release(monkeypatch):
    monkeypatch.setattr(EnsemblRelease, "required_local_files_exist", lambda self: self.release == 74)
    genome, = select("--reference-name", "GRCh37")
    assert genome.release == 75


def test_reference_selection_applies_to_custom_mirror():
    genome, = select(
        "--reference-name", "grch37", "--custom-mirror", "https://example.invalid"
    )
    assert genome.reference_name == "GRCh37"
    assert genome.annotation_version == 75
    sources = genome.to_dict()
    assert sources["gtf_path_or_url"] == "https://example.invalid/Homo_sapiens.GRCh37.75.gtf.gz"
    assert all(
        url.startswith("https://example.invalid/") and ".GRCh37." in url
        for url in sources["transcript_fasta_paths_or_urls"] + sources["protein_fasta_paths_or_urls"]
    )


def test_custom_source_reference_is_not_resolved_as_ensembl():
    genome, = select(
        "--reference-name", "MyAssembly", "--annotation-name", "custom",
        "--gtf", "https://example.invalid/custom.gtf",
    )
    assert type(genome) is Genome
    assert genome.reference_name == "MyAssembly"
    assert genome.annotation_name == "custom"


def test_custom_source_install_without_gtf():
    # A FASTA-only (or protein-only) custom install must get past argument
    # handling even when no GTF is supplied (issue #400): the Genome class
    # supports omitted GTF sources, so the CLI must not crash on a None
    # --gtf value before the install begins.
    genome, = select(
        "--reference-name", "MyAssembly", "--annotation-name", "custom",
        "--transcript-fasta", "/tmp/transcripts.fa",
    )
    assert type(genome) is Genome
    assert not genome.requires_gtf
    sources = genome.to_dict()
    assert sources["gtf_path_or_url"] is None
    assert sources["transcript_fasta_paths_or_urls"] == ["/tmp/transcripts.fa"]


def test_no_reference_preserves_species_release_combinations():
    genomes = select("--species", "human", "mouse", "--release", "75", "76")
    assert [(g.species.latin_name, g.release) for g in genomes] == [
        ("homo_sapiens", 75), ("homo_sapiens", 76),
        ("mus_musculus", 75), ("mus_musculus", 76),
    ]


@pytest.fixture
def cli_calls(monkeypatch, tmp_path):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(shell, "configure_logging", lambda: None)
    calls = []
    monkeypatch.setattr(Genome, "download", lambda self, **kw: calls.append(("download", self)))
    monkeypatch.setattr(Genome, "index", lambda self, **kw: calls.append(("index", self)))
    monkeypatch.setattr(shell, "_delete_genome_files", lambda *args: calls.append(("delete", args)))
    return calls


def test_install_acts_on_requested_assembly(monkeypatch, cli_calls):
    monkeypatch.setattr(sys, "argv", ["pyensembl", "install", "--reference-name", "grch37"])
    shell.run()
    assert [action for action, genome in cli_calls] == ["download", "index"]
    assert all(
        isinstance(genome, EnsemblRelease) and genome.reference_name == "GRCh37"
        and genome.release == 75
        for action, genome in cli_calls
    )


@pytest.mark.parametrize("action", ["install", "delete-all-files", "delete-index-files"])
@pytest.mark.parametrize(
    "arguments, message",
    [
        (["--reference-name", "unknown", "--release", "75"], "not found"),
        (["--reference-name", "GRCh37", "--release", "54"], "55-75"),
        (["--reference-name", "GRCh37", "--release", "76"], "55-75"),
        (["--reference-name", "GRCh37", "--release", "75", "76"], "55-75"),
        (["--reference-name", "GRCh37", "--release", "75", "--species", "mouse"], "mus_musculus"),
        (["--reference-name", "GRCh37", "--release", "75", "--species", "human", "mouse"], "mus_musculus"),
        (["--reference-name", "GRCh37", "--release", "76", "--custom-mirror", "https://example.invalid"], "55-75"),
    ],
)
def test_conflicting_selection_fails_before_any_action(
    monkeypatch, capsys, cli_calls, action, arguments, message,
):
    monkeypatch.setattr(sys, "argv", ["pyensembl", action, *arguments])
    with pytest.raises(SystemExit) as error:
        shell.run()
    assert error.value.code == 2
    assert message in capsys.readouterr().err
    assert cli_calls == []


@pytest.mark.parametrize("action", ["delete-all-files", "delete-index-files"])
def test_reference_does_not_replace_explicit_deletion_release(
    monkeypatch, capsys, cli_calls, action,
):
    monkeypatch.setattr(sys, "argv", ["pyensembl", action, "--reference-name", "GRCh37"])
    with pytest.raises(SystemExit) as error:
        shell.run()
    assert error.value.code == 2
    assert "requires an explicit --release" in capsys.readouterr().err
    assert cli_calls == []

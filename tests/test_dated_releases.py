"""Dated releases: annotation dates on the new Ensembl platform (#447)."""

import gzip
import os
from pathlib import Path
import re
import urllib.request

import pytest

from pyensembl import EnsemblRelease, shell
from pyensembl.ensembl_url_templates import make_dated_release_urls
from pyensembl.ensembl_versions import check_release_number
from pyensembl.species import Species

from .test_ensembl_annotation import bgzf

PLATFORM = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"


def test_install_index_query_and_reopen_dated_release(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    mirror = tmp_path / "mirror"
    root = mirror / "GCA/000/001/405/29/ensembl/2026_04"
    (root / "geneset").mkdir(parents=True)
    (root / "genome").mkdir()
    attrs = 'gene_id "ENSG1"; gene_name "TEST"; transcript_id "ENST1"; transcript_version "2";'
    # 2026-era dumps are coordinate-sorted: exons can precede their gene.
    rows = (
        f'1\tensembl\texon\t3\t6\t.\t-\t.\t{attrs} exon_id "e2"; exon_number "2";\n'
        f'1\tensembl\texon\t9\t12\t.\t-\t.\t{attrs} exon_id "e1"; exon_number "1";\n'
        f'1\tensembl\ttranscript\t3\t12\t.\t-\t.\t{attrs}\n'
        '1\tensembl\tgene\t3\t12\t.\t-\t.\tgene_id "ENSG1"; gene_name "TEST";\n'
    )
    # GRCh38 selects the counterpart of the complete patch/haplotype GTF.
    (root / "geneset/genes-including_alt.gtf.gz").write_bytes(gzip.compress(rows.encode()))
    (root / "geneset/cdna.fa.bgz").write_bytes(bgzf(b">ENST1.2\nCCCCGGTT\n"))
    (root / "geneset/pep.fa.bgz").write_bytes(bgzf(b">ENSP1.1\nMP\n"))
    (root / "genome/unmasked.fa.bgz").write_bytes(bgzf(b">1\nAAAACCCCGGGG\n"))

    genome = EnsemblRelease("2026_04", server=mirror.as_uri(), genome_fasta=True)
    assert genome.reference_name == "GRCh38"
    directory = Path(genome.download_cache.cache_directory_path)
    assert directory.parts[-2:] == ("GRCh38", "ensembl2026_04")
    assert not directory.exists()  # Construction is offline.
    genome.download()
    genome.index()
    assert genome.gene_by_id("ENSG1").strand == "-"
    transcript = genome.transcript_by_id("ENST1")
    assert transcript.sequence == "CCCCGGTT"
    assert [exon.id for exon in transcript.exons] == ["e1", "e2"]
    assert genome.sequence("1", 9, 12, strand="-") == "CCCC"
    reopened = EnsemblRelease.from_json(genome.to_json())
    assert reopened.release == "2026_04"
    assert reopened.transcript_by_id("ENST1").sequence == "CCCCGGTT"


@pytest.mark.parametrize(
    "species, release, reference, directory, gtf",
    [
        ("human", "2026_04", "GRCh38", "GCA/000/001/405/29/ensembl/2026_04",
         "genes-including_alt.gtf.gz"),
        ("mouse", "2025_12", "GRCm39", "GCA/000/001/635/9/ensembl/2025_12", "genes.gtf.gz"),
        # The new platform labels this accession BDGP6.46; keep pyensembl's name.
        ("fly", "2022_07", "BDGP6.54", "GCA/000/001/215/4/flybase/2022_07", "genes.gtf.gz"),
        ("arabidopsis_thaliana", "2010_09", "TAIR10",
         "GCA/000/001/735/1/community_araport11/2010_09", "genes.gtf.gz"),
    ],
)
def test_dated_release_files(species, release, reference, directory, gtf):
    genome = EnsemblRelease(release, species=species, genome_fasta=True, genome_fasta_mask="soft")
    base = "%s/%s/" % (PLATFORM, directory)
    assert genome.reference_name == reference
    assert genome.gtf_url == base + "geneset/" + gtf
    # cDNA covers every biotype, so there is no separate ncRNA FASTA.
    assert genome.transcript_fasta_urls == [base + "geneset/cdna.fa.bgz"]
    assert genome.protein_fasta_urls == [base + "geneset/pep.fa.bgz"]
    assert genome.genome_fasta_urls == [base + "genome/softmasked.fa.bgz"]


def test_dated_and_numbered_releases_share_the_reference_directory():
    dated, numbered = EnsemblRelease("2026_04"), EnsemblRelease(116)
    assert dated != numbered
    assert dated.reference_name == numbered.reference_name == "GRCh38"
    dated_directory = Path(dated.download_cache.cache_directory_path)
    numbered_directory = Path(numbered.download_cache.cache_directory_path)
    assert dated_directory.parent == numbered_directory.parent
    assert (dated_directory.name, numbered_directory.name) == ("ensembl2026_04", "ensembl116")


def test_dated_release_descriptions_quote_the_date():
    genome = EnsemblRelease("2026_04", species="mouse")
    assert str(genome) == "EnsemblRelease(release='2026_04', species='mus_musculus')"
    assert genome.install_string() == "pyensembl install --release 2026_04 --species mus_musculus"
    # A bare 2026_04 would be the Python integer 202604.
    assert "EnsemblRelease('2026_04', species='mus_musculus', genome_fasta=True)" in (
        genome._genome_fasta_setup_hint()
    )


@pytest.mark.parametrize("release", ["2026-07", "2026-09-22"])
def test_website_release_labels_are_not_annotation_dates(release):
    with pytest.raises(ValueError, match="website release label"):
        EnsemblRelease(release)


@pytest.mark.parametrize("release", ["2026_13", "2026_4", "26_04"])
def test_invalid_annotation_dates_rejected(release):
    with pytest.raises(ValueError):
        EnsemblRelease(release)


def test_release_numbers_do_not_read_dates_as_integers():
    # int("2026_04") == 202604 because Python accepts digit separators.
    with pytest.raises(ValueError, match="Invalid Ensembl release"):
        check_release_number("2026_04")
    assert EnsemblRelease("116").release == 116


def test_species_without_dated_releases_rejected():
    with pytest.raises(ValueError, match="No dated Ensembl releases"):
        EnsemblRelease("2026_04", species="toxoplasma_gondii")


@pytest.mark.parametrize(
    "options",
    [{"genome_fasta_type": "primary_assembly"}, {"genome_fasta_mask": "invalid"}],
)
def test_unpublished_reference_dna_rejected(options):
    with pytest.raises(ValueError):
        EnsemblRelease("2026_04", genome_fasta=True, **options)


def select(*arguments):
    return shell.collect_selected_genomes(shell.parser.parse_args(["install", *arguments]))


def test_cli_installs_numbered_and_dated_releases():
    genomes = select("--release", "116", "2026_04", "--species", "human")
    assert [genome.release for genome in genomes] == [116, "2026_04"]
    genome, = select("--reference-name", "GRCh38", "--release", "2026_04")
    assert genome.reference_name == "GRCh38"
    with pytest.raises(ValueError, match="provides GRCh38, not GRCh37"):
        select("--reference-name", "GRCh37", "--release", "2026_04")


def test_cli_rejects_website_release_labels(capsys):
    with pytest.raises(SystemExit):
        shell.parser.parse_args(["install", "--release", "2026-07"])
    assert "website release label" in capsys.readouterr().err


def test_list_labels_dated_release_directories():
    assert shell._other_genome_labels("GRCh38", "ensembl2026_04") == ("human", "2026_04")
    assert shell._other_genome_labels("GRCh38", "ensembl116") == ("human", "116")


def _dated_release_problems(species):
    accession, provider = species.dated_releases
    directory = make_dated_release_urls(accession, provider, "").directory
    try:
        with urllib.request.urlopen(directory, timeout=60) as response:
            dates = re.findall(r'href="(\d{4}_\d{2})/"', response.read().decode())
    except OSError as error:
        return ["%s: %s" % (directory, error)]
    if not dates:
        return ["%s: no dated releases" % directory]
    genome = EnsemblRelease(max(dates), species=species, genome_fasta=True)
    problems = []
    for url in [genome.gtf_url] + genome.transcript_fasta_urls + genome.protein_fasta_urls + (
        genome.genome_fasta_urls
    ):
        try:
            urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60)
        except OSError as error:
            problems.append("%s: %s" % (url, error))
    return problems


@pytest.mark.skipif(
    not os.environ.get("PYENSEMBL_NETWORK_TESTS"),
    reason="set PYENSEMBL_NETWORK_TESTS=1 to check the new Ensembl platform",
)
def test_every_dated_release_table_entry_publishes_its_files():
    from concurrent.futures import ThreadPoolExecutor

    species = [s for s in Species._latin_names_to_species.values() if s.dated_releases]
    with ThreadPoolExecutor(6) as pool:
        problems = [p for found in pool.map(_dated_release_problems, species) for p in found]
    assert problems == []

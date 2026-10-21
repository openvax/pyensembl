"""Assembly release ranges agree with Ensembl's archive (#396).

Expected filenames were read from the archive's GTF and FASTA directories.
Every configured assembly boundary is checked against the archive URL layout.
"""

import os
import re
from urllib.parse import urlsplit

import pytest

from pyensembl import EnsemblRelease
from pyensembl.species import Species, find_species_by_name


@pytest.mark.parametrize(
    "species,release,gtf,cdna",
    [
        # The reported case: release 97 is GRCg6a, not Gallus_gallus-5.0.
        ("chicken", 94, "Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz",
         "Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz"),
        ("chicken", 97, "Gallus_gallus.GRCg6a.97.gtf.gz",
         "Gallus_gallus.GRCg6a.cdna.all.fa.gz"),
        ("chicken", 107, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.107.gtf.gz",
         "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"),
        ("cat", 114, "Felis_catus.F.catus_Fca126_mat1.0.114.gtf.gz",
         "Felis_catus.F.catus_Fca126_mat1.0.cdna.all.fa.gz"),
        ("macaca_mulatta", 97, "Macaca_mulatta.Mmul_8.0.1.97.gtf.gz",
         "Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz"),
        ("pig", 89, "Sus_scrofa.Sscrofa10.2.89.gtf.gz",
         "Sus_scrofa.Sscrofa10.2.cdna.all.fa.gz"),
        ("drosophila_melanogaster", 111, "Drosophila_melanogaster.BDGP6.46.111.gtf.gz",
         "Drosophila_melanogaster.BDGP6.46.cdna.all.fa.gz"),
        ("naked_mole_rat", 110,
         "Heterocephalus_glaber_female.Naked_mole-rat_maternal.110.gtf.gz",
         "Heterocephalus_glaber_female.Naked_mole-rat_maternal.cdna.all.fa.gz"),
        # #407: Ensembl renamed dog to canis_lupus_familiaris in release 100.
        ("dog", 99, "Canis_familiaris.CanFam3.1.99.gtf.gz",
         "Canis_familiaris.CanFam3.1.cdna.all.fa.gz"),
        ("dog", 100, "Canis_lupus_familiaris.CanFam3.1.100.gtf.gz",
         "Canis_lupus_familiaris.CanFam3.1.cdna.all.fa.gz"),
        ("dog", 105, "Canis_lupus_familiaris.ROS_Cfam_1.0.105.gtf.gz",
         "Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz"),
        ("plasmodium_falciparum", 60, "Plasmodium_falciparum.ASM276v2.60.gtf.gz",
         "Plasmodium_falciparum.ASM276v2.cdna.all.fa.gz"),
        ("plasmodium_falciparum", 61, "Plasmodium_falciparum.GCA000002765v3.61.gtf.gz",
         "Plasmodium_falciparum.GCA000002765v3.cdna.all.fa.gz"),
        ("toxoplasma_gondii", 61, "Toxoplasma_gondii.GCA000006565v2.61.gtf.gz",
         "Toxoplasma_gondii.GCA000006565v2.cdna.all.fa.gz"),
        # Ensembl 116 republished Ensembl Genomes 63 GTFs for non-vertebrates.
        ("drosophila_melanogaster", 116, "Drosophila_melanogaster.BDGP6.54.63.gtf.gz",
         "Drosophila_melanogaster.BDGP6.54.cdna.all.fa.gz"),
        ("yeast", 116, "Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz",
         "Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"),
        ("nematode", 116, "Caenorhabditis_elegans.WBcel235.63.gtf.gz",
         "Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"),
    ],
)
def test_assembly_transitions_name_archived_files(species, release, gtf, cdna):
    genome = EnsemblRelease(release, species=species)
    assert os.path.basename(genome.gtf_url) == gtf
    assert os.path.basename(genome.transcript_fasta_urls[0]) == cdna


@pytest.mark.parametrize(
    "species,release",
    [
        ("syrian_hamster", 89),  # Added in release 90.
        ("meriones_unguiculatus", 95),  # Added in release 96.
        ("mus_musculus_balbcj", 87),  # 87-91 kept release 86's filenames.
    ],
)
def test_releases_missing_from_the_archive_are_rejected(species, release):
    with pytest.raises(ValueError, match="No genome for"):
        EnsemblRelease(release, species=species)


@pytest.mark.parametrize(
    "species,release,gtf",
    [
        ("dingo", 95, "Canis_lupus_dingo.ASM325472v1.95.gtf.gz"),
        ("dog_basenji", 99, "Canis_lupus_familiarisbasenji.Basenji_breed-1.1.99.gtf.gz"),
        ("dog_great_dane", 99, "Canis_lupus_familiarisgreatdane.UMICH_Zoey_3.1.99.gtf.gz"),
        ("dog_boxer", 105, "Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.105.gtf.gz"),
        ("dog_german_shepherd", 107, "Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.107.gtf.gz"),
    ],
)
def test_every_dog_genome_ensembl_publishes(species, release, gtf):
    assert os.path.basename(EnsemblRelease(release, species=species).gtf_url) == gtf
    with pytest.raises(ValueError, match="No genome for"):
        EnsemblRelease(release - 1, species=species)


def test_renamed_species_are_found_by_every_name():
    dog = find_species_by_name("dog")
    for name in ("canis_familiaris", "canis_lupus_familiaris", "Canis lupus familiaris"):
        assert find_species_by_name(name) is dog
    assert dog.ensembl_name(99) == "canis_familiaris"
    assert dog.ensembl_name(100) == "canis_lupus_familiaris"


def _boundary_urls():
    for species in Species._latin_names_to_species.values():
        releases = set()
        for start, end in species.reference_assemblies.values():
            releases |= {start, end}
        for renamed_in in species.ensembl_names:  # Both sides of a rename.
            releases |= {renamed_in - 1, renamed_in}
        for release in sorted(releases):
            genome = EnsemblRelease(release, species=species)
            yield genome.gtf_url
            yield genome.transcript_fasta_urls[0]
            yield genome.protein_fasta_urls[0]


def test_every_assembly_boundary_uses_archive_layout():
    urls = list(_boundary_urls())
    assert len(urls) % 3 == 0
    assert len(urls) == len(set(urls))
    for gtf, cdna, pep in zip(urls[::3], urls[1::3], urls[2::3]):
        gtf_url, cdna_url, pep_url = map(urlsplit, (gtf, cdna, pep))
        assert {url.scheme for url in (gtf_url, cdna_url, pep_url)} == {"https"}
        assert len({url.netloc for url in (gtf_url, cdna_url, pep_url)}) == 1
        host = gtf_url.netloc
        assert host in {"ftp.ensembl.org", "ftp.ensemblgenomes.ebi.ac.uk"}

        root, gtf_file = gtf_url.path.split("/gtf/")
        if host == "ftp.ensembl.org":
            assert re.fullmatch(r"/pub/release-\d+", root)
        else:
            assert re.fullmatch(r"/pub/release-\d+/(plants|fungi|metazoa|protists)", root)
        species, gtf_name = gtf_file.split("/")
        assert gtf_name.startswith(species.capitalize() + ".")
        assert re.search(r"\.\d+(?:\.chr_patch_hapl_scaff)?\.gtf\.gz$", gtf_name)
        assert cdna_url.path.startswith(f"{root}/fasta/{species}/cdna/")
        assert pep_url.path.startswith(f"{root}/fasta/{species}/pep/")

        gtf_stem = gtf_name.removesuffix(".gtf.gz")
        cdna_stem = cdna_url.path.rsplit("/", 1)[-1].removesuffix(".cdna.all.fa.gz")
        pep_stem = pep_url.path.rsplit("/", 1)[-1].removesuffix(".pep.all.fa.gz")
        assert cdna_url.path.endswith(".cdna.all.fa.gz")
        assert pep_url.path.endswith(".pep.all.fa.gz")
        assert cdna_stem == pep_stem
        assert gtf_stem == cdna_stem or gtf_stem.startswith(cdna_stem + ".")

"""
Tests for #298: more species support across Ensembl Genomes divisions
(plants / fungi / metazoa / protists) and the URL-template branch that
routes them through ``ftp.ensemblgenomes.ebi.ac.uk``.
"""

from pyensembl import EnsemblRelease
from pyensembl.ensembl_url_templates import (
    ENSEMBL_FTP_SERVER,
    ENSEMBL_GENOMES_FTP_SERVER,
)
from pyensembl.ensembl_versions import MAX_ENSEMBL_GENOMES_RELEASE
from pyensembl.species import (
    Species,
    anopheles_gambiae,
    arabidopsis_thaliana,
    aspergillus_nidulans,
    candida_albicans,
    fission_yeast,
    maize,
    plasmodium_falciparum,
    rice,
    tomato,
    toxoplasma_gondii,
    wheat,
)


GENOMES_SPECIES_BY_DIVISION = {
    "plants": [arabidopsis_thaliana, rice, wheat, maize, tomato],
    "fungi": [fission_yeast, aspergillus_nidulans, candida_albicans],
    "metazoa": [anopheles_gambiae],
    "protists": [plasmodium_falciparum, toxoplasma_gondii],
}


def test_ensembl_genomes_species_have_expected_divisions():
    for division, species_list in GENOMES_SPECIES_BY_DIVISION.items():
        for species in species_list:
            assert species.division == division, (
                "%s should be tagged division=%s, got %s"
                % (species.latin_name, division, species.division)
            )
            assert species.ensembl_genomes is True


def test_ensembl_genomes_species_route_through_genomes_server():
    for species_list in GENOMES_SPECIES_BY_DIVISION.values():
        for species in species_list:
            release = EnsemblRelease(
                release=MAX_ENSEMBL_GENOMES_RELEASE, species=species
            )
            assert release.server == ENSEMBL_GENOMES_FTP_SERVER
            assert release.gtf_url.startswith(ENSEMBL_GENOMES_FTP_SERVER)
            assert "/%s/gtf/%s/" % (species.division, species.latin_name) in (
                release.gtf_url
            )
            for url in release.transcript_fasta_urls + release.protein_fasta_urls:
                assert url.startswith(ENSEMBL_GENOMES_FTP_SERVER)
                assert "/%s/fasta/%s/" % (species.division, species.latin_name) in url


def test_genomes_filenames_use_modern_layout_regardless_of_release():
    # Older Ensembl releases (<= 75) used a different filename pattern on
    # main Ensembl. Ensembl Genomes uses the modern pattern at every release.
    release = EnsemblRelease(release=40, species=wheat)
    for url in release.transcript_fasta_urls + release.protein_fasta_urls:
        assert ".40." not in url, (
            "Ensembl Genomes FASTA filenames should not embed the release"
            " number, got %s" % url
        )


def test_existing_non_vertebrate_species_stay_on_main_ensembl():
    """yeast / drosophila / C. elegans are tagged taxonomically as
    fungi/metazoa but their reference_assemblies are registered against main
    Ensembl release numbers, so they must keep routing through ftp.ensembl.org.
    """
    for latin_name in (
        "saccharomyces_cerevisiae",
        "drosophila_melanogaster",
        "caenorhabditis_elegans",
    ):
        species = Species._latin_names_to_species[latin_name]
        assert species.ensembl_genomes is False
        release = EnsemblRelease(release=110, species=species)
        assert release.server == ENSEMBL_FTP_SERVER
        assert release.gtf_url.startswith(ENSEMBL_FTP_SERVER)


def test_is_plant_kwarg_back_compat_implies_ensembl_genomes():
    species = Species.register(
        latin_name="_test_legacy_is_plant_species",
        synonyms=["_test_legacy"],
        reference_assemblies={"FakeAsm": (40, MAX_ENSEMBL_GENOMES_RELEASE)},
        is_plant=True,
    )
    try:
        assert species.is_plant is True
        assert species.division == "plants"
        assert species.ensembl_genomes is True
    finally:
        # Remove the test fixture so it doesn't leak into other tests.
        Species._latin_names_to_species.pop(species.latin_name, None)
        for synonym in species.synonyms:
            Species._common_names_to_species.pop(synonym, None)
        for ref in species.reference_assemblies:
            Species._reference_names_to_species.pop(ref, None)

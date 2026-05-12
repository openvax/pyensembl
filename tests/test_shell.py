from pyensembl.shell import (
    all_combinations_of_ensembl_genomes,
    format_available_species,
    parser,
)
from .common import eq_


def test_genome_selection_grch38():
    args = parser.parse_args(["install", "--release", "100", "--species", "human"])
    genomes = all_combinations_of_ensembl_genomes(args)
    assert len(genomes) == 1
    genome = genomes[0]
    eq_(genome.species.latin_name, "homo_sapiens")
    eq_(genome.release, 100)


def test_available_action_parses():
    args = parser.parse_args(["available"])
    eq_(args.action, "available")


def test_format_available_species_includes_human_and_assemblies():
    output = format_available_species(use_color=False)
    # human is registered with common name "human" and three reference assemblies
    assert "homo_sapiens" in output
    assert "human" in output
    assert "GRCh38" in output
    assert "GRCh37" in output
    # mouse should also appear
    assert "mus_musculus" in output
    assert "GRCm38" in output


def test_format_available_species_grouped_by_division():
    output = format_available_species(use_color=False)
    # Every populated division should have a section header.
    assert "── Vertebrates " in output
    assert "── Invertebrates " in output
    assert "── Plants " in output
    assert "── Fungi " in output
    # Section ordering: Vertebrates before Invertebrates before Plants before Fungi.
    v = output.index("── Vertebrates ")
    i = output.index("── Invertebrates ")
    p = output.index("── Plants ")
    f = output.index("── Fungi ")
    assert v < i < p < f
    # Yeast is now classified as fungi; drosophila/c. elegans as metazoa.
    assert output.index("yeast") > f
    drosophila_pos = output.index("drosophila")
    assert i < drosophila_pos < p


def test_format_available_species_no_color_has_no_escape_codes():
    output = format_available_species(use_color=False)
    assert "\x1b[" not in output


def test_format_available_species_collapses_single_release():
    # NCBI36 only exists in Ensembl release 54; verify it renders as "54"
    # rather than "54–54".
    output = format_available_species(use_color=False)
    assert "NCBI36" in output
    ncbi36_line = next(
        line for line in output.splitlines() if "NCBI36" in line
    )
    assert "54–54" not in ncbi36_line
    assert "54" in ncbi36_line

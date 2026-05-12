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
    output = format_available_species()
    # human is registered with common name "human" and three reference assemblies
    assert "homo_sapiens" in output
    assert "human" in output
    assert "GRCh38" in output
    assert "GRCh37" in output
    # mouse should also appear
    assert "mus_musculus" in output
    assert "GRCm38" in output

from nose.tools import eq_, with_setup
from pyensembl.shell import parser, all_combinations_of_ensembl_genomes
from pyensembl import ensembl_grch38


def test_genome_selection_grch38():
    args = parser.parse_args(["install", "--release", "100", "--species", "human"])
    genomes = all_combinations_of_ensembl_genomes(args)
    assert len(genomes) == 1
    genome = genomes[0]
    eq_(genome.species.latin_name, "homo_sapiens")
    eq_(genome.release, 100)

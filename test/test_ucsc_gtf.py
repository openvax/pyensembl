from __future__ import absolute_import

from pyensembl import Genome, GenomeSource
from nose.tools import eq_
from os.path import join, dirname

UCSC_PATH = join(dirname(__file__), "data/ucsc.small.gtf")

def test_ucsc_parsing():
    genome = Genome(GenomeSource(gtf_path=UCSC_PATH))
    genome.index()
    eq_(len(genome.genes()), 7)
    eq_(len(genome.transcripts()), 7)

    gene_uc001aak4 = genome.gene_by_id("uc001aak.4")
    eq_(gene_uc001aak4.id, "uc001aak.4")
    eq_(gene_uc001aak4.name, None)
    eq_(gene_uc001aak4.biotype, None)

    gene_1_17369 = genome.genes_at_locus(1, 17369)
    eq_(gene_1_17369[0].id, "uc031tla.1")

    transcript_1_30564 = genome.transcripts_at_locus(1, 30564)
    eq_(transcript_1_30564[0].id, "uc057aty.1")

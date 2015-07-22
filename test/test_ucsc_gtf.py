from __future__ import absolute_import

from . import data_path
from pyensembl import Genome
from nose.tools import eq_

UCSC_GENCODE_PATH = data_path("gencode.ucsc.small.gtf")
UCSC_REFSEQ_PATH = data_path("refseq.ucsc.small.gtf")

def test_ucsc_gencode():
    """
    Testing with a small GENCODE GTF file downloaded from 
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    genome = Genome("GRCh38", gtf_path_or_url=UCSC_GENCODE_PATH)
    genome.install()
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

def test_ucsc_refseq():
    """
    Testing with a small RefSeq GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    genome = Genome("GRCh38", gtf_path_or_url=UCSC_REFSEQ_PATH)
    genome.install()
    eq_(len(genome.genes()), 2)
    eq_(len(genome.transcripts()), 2)

    genes_at_locus = genome.genes_at_locus(1, 67092176)
    eq_(len(genes_at_locus), 2)
    ids = set([gene.id for gene in genes_at_locus])
    eq_(set(["NM_001276352", "NR_075077"]), ids)

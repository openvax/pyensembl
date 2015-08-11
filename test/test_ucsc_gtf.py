from __future__ import absolute_import

from pyensembl import Genome, GTF
from nose.tools import eq_

from .data import data_path

UCSC_GENCODE_PATH = data_path("gencode.ucsc.small.gtf")
UCSC_REFSEQ_PATH = data_path("refseq.ucsc.small.gtf")

def test_ucsc_gencode_gtf():
    gtf = GTF(UCSC_GENCODE_PATH)
    df = gtf.dataframe(save_to_disk=False)
    exons = df[df["feature"] == "exon"]
    # expect 12 exons from the dataframe
    eq_(len(exons), 12)

def test_ucsc_gencode_genome():
    """
    Testing with a small GENCODE GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="ucsc_test",
        gtf_path_or_url=UCSC_GENCODE_PATH)
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

def test_ucsc_refseq_gtf():
    """
    Test GTF object with a small RefSeq GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    gtf = GTF(UCSC_REFSEQ_PATH)
    df = gtf.dataframe(save_to_disk=False)
    exons = df[df["feature"] == "exon"]
    # expect 16 exons from the GTF
    eq_(len(exons), 16)

def test_ucsc_refseq_genome():
    """
    Test Genome object with a small RefSeq GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="ucsc_test",
        gtf_path_or_url=UCSC_REFSEQ_PATH)
    eq_(len(genome.genes()), 2)
    eq_(len(genome.transcripts()), 2)

    genes_at_locus = genome.genes_at_locus(1, 67092176)
    eq_(len(genes_at_locus), 2)
    ids = set([gene.id for gene in genes_at_locus])
    eq_(set(["NM_001276352", "NR_075077"]), ids)

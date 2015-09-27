from __future__ import absolute_import

from pyensembl import Genome, GTF
from nose.tools import eq_

from .common import TemporaryDirectory
from .data import data_path

UCSC_GENCODE_PATH = data_path("gencode.ucsc.small.gtf")
UCSC_REFSEQ_PATH = data_path("refseq.ucsc.small.gtf")

def test_ucsc_gencode_gtf():
    with TemporaryDirectory() as tmpdir:
        gtf = GTF(
            UCSC_GENCODE_PATH,
            cache_directory_path=tmpdir)
        df = gtf.dataframe(save_to_disk=False)
        exons = df[df["feature"] == "exon"]
        # expect 12 exons from the dataframe
        assert len(exons) == 12, \
            "Expected 12 exons, got %d: %s" % (len(exons), exons)
        df2 = gtf.dataframe(save_to_disk=True)
        assert len(df) == len(df2), "Got different length DataFrame"
        assert list(df.columns) == list(df2.columns)


def test_ucsc_gencode_genome():
    """
    Testing with a small GENCODE GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    with TemporaryDirectory() as tmpdir:
        genome = Genome(
            reference_name="GRCh38",
            annotation_name="ucsc_test",
            gtf_path_or_url=UCSC_GENCODE_PATH,
            cache_directory_path=tmpdir)
        genome.index()
        genes = genome.genes()
        for gene in genes:
            assert gene.id, \
                "Gene with missing ID in %s" % (genome.gtf.dataframe(),)
        assert len(genes) == 7, \
            "Expected 7 genes, got %d: %s" % (
                len(genes), genes)
        transcripts = genome.transcripts()
        for transcript in transcripts:
            assert transcript.id, \
                "Transcript with missing ID in %s" % (genome.gtf.dataframe(),)
        assert len(transcripts) == 7, \
            "Expected 7 transcripts, got %d: %s" % (
                len(transcripts), transcripts)

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
    with TemporaryDirectory() as tmpdir:
        gtf = GTF(
            UCSC_REFSEQ_PATH,
            cache_directory_path=tmpdir)
        df = gtf.dataframe(save_to_disk=False)
        exons = df[df["feature"] == "exon"]
        # expect 16 exons from the GTF
        assert len(exons) == 16, \
            "Expected 16 exons, got %d: %s" % (
                len(exons), exons)
        df2 = gtf.dataframe(save_to_disk=True)
        assert len(df) == len(df2), "Got different length DataFrame"
        assert list(df.columns) == list(df2.columns)

def test_ucsc_refseq_genome():
    """
    Test Genome object with a small RefSeq GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    with TemporaryDirectory() as tmpdir:
        genome = Genome(
            reference_name="GRCh38",
            annotation_name="ucsc_test",
            gtf_path_or_url=UCSC_REFSEQ_PATH,
            cache_directory_path=tmpdir)
        genome.index()
        genes = genome.genes()
        for gene in genes:
            assert gene.id, \
                "Gene with missing ID in %s" % (genome.gtf.dataframe(),)
        assert len(genes) == 2, \
            "Expected 2 genes, got %d: %s" % (
                len(genes), genes)
        transcripts = genome.transcripts()
        for transcript in transcripts:
            assert transcript.id, \
                 "Transcript with missing ID in %s" % (genome.gtf.dataframe(),)
        assert len(transcripts) == 2, \
            "Expected 2 transcripts, got %d: %s" % (
                len(transcripts), transcripts)
        genes_at_locus = genome.genes_at_locus(1, 67092176)
        assert len(genes_at_locus) == 2, \
            "Expected 2 genes at locus chr1:67092176, got %d: %s" % (
                len(genes_at_locus), genes_at_locus)
        ids = set([gene.id for gene in genes_at_locus])
        eq_(set(["NM_001276352", "NR_075077"]), ids)

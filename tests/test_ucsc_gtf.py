from pyensembl import Genome, Database

from .common import TemporaryDirectory, eq_
from .data import data_path

UCSC_GENCODE_PATH = data_path("gencode.ucsc.small.gtf")
UCSC_REFSEQ_PATH = data_path("refseq.ucsc.small.gtf")


def test_ucsc_gencode_gtf():
    with TemporaryDirectory() as tmpdir:
        db = Database(UCSC_GENCODE_PATH, cache_directory_path=tmpdir)
        df = db._load_gtf_as_dataframe()
        exons = df[df["feature"] == "exon"]
        # expect 12 exons from the dataframe
        assert len(exons) == 12, "Expected 12 exons, got %d: %s" % (len(exons), exons)


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
            cache_directory_path=tmpdir,
        )
        genome.index()
        genes = genome.genes()
        for gene in genes:
            assert gene.id, "Gene with missing ID in %s" % (genome,)
        assert len(genes) == 7, "Expected 7 genes, got %d: %s" % (len(genes), genes)
        transcripts = genome.transcripts()
        for transcript in transcripts:
            assert transcript.id, "Transcript with missing ID in %s" % (genome,)
        assert len(transcripts) == 7, "Expected 7 transcripts, got %d: %s" % (
            len(transcripts),
            transcripts,
        )

        gene_uc001aak4 = genome.gene_by_id("uc001aak.4")
        eq_(gene_uc001aak4.id, "uc001aak.4")
        eq_(gene_uc001aak4.name, None)
        eq_(gene_uc001aak4.biotype, None)

        gene_1_17369 = genome.genes_at_locus("chr1", 17369)
        eq_(gene_1_17369[0].id, "uc031tla.1")

        transcript_1_30564 = genome.transcripts_at_locus("chr1", 30564)
        eq_(transcript_1_30564[0].id, "uc057aty.1")


def test_transcript_exons_without_exon_id():
    """
    Regression test: older Ensembl releases (e.g. release 54) and other
    GTFs omit the exon_id attribute while still providing exon_number.
    Transcript.exons previously raised
    sqlite3.OperationalError: no such column: exon_id; it now falls back
    to building Exon objects directly from the exon row with a
    synthesized per-transcript ID.
    """
    import os

    with TemporaryDirectory() as tmpdir:
        gtf_path = os.path.join(tmpdir, "no_exon_id.gtf")
        with open(gtf_path, "w") as f:
            # minimal Ensembl-style GTF: transcript + 2 ordered exons,
            # with exon_number but no exon_id (as in Ensembl release 54)
            f.write(
                '1\ttest\ttranscript\t100\t500\t.\t+\t.\t'
                'gene_id "G1"; transcript_id "T1"; gene_name "FN1";\n'
                '1\ttest\texon\t100\t200\t.\t+\t.\t'
                'gene_id "G1"; transcript_id "T1"; exon_number "1"; gene_name "FN1";\n'
                '1\ttest\texon\t300\t500\t.\t+\t.\t'
                'gene_id "G1"; transcript_id "T1"; exon_number "2"; gene_name "FN1";\n'
            )
        genome = Genome(
            reference_name="GRCh38",
            annotation_name="no_exon_id_test",
            gtf_path_or_url=gtf_path,
            cache_directory_path=tmpdir,
        )
        genome.index()
        assert not genome.db.column_exists("exon", "exon_id"), (
            "Test fixture unexpectedly has exon_id — update this test"
        )
        transcript = genome.transcript_by_id("T1")
        exons = transcript.exons
        eq_(len(exons), 2)
        eq_(exons[0].id, "T1_exon_1")
        eq_(exons[0].start, 100)
        eq_(exons[0].end, 200)
        eq_(exons[1].id, "T1_exon_2")
        eq_(exons[1].start, 300)
        eq_(exons[1].end, 500)


def test_ucsc_refseq_gtf():
    """
    Test GTF object with a small RefSeq GTF file downloaded from
    http://genome.ucsc.edu/cgi-bin/hgTables
    """
    with TemporaryDirectory() as tmpdir:
        db = Database(UCSC_REFSEQ_PATH, cache_directory_path=tmpdir)
        df = db._load_gtf_as_dataframe()
        exons = df[df["feature"] == "exon"]
        # expect 16 exons from the GTF
        assert len(exons) == 16, "Expected 16 exons, got %d: %s" % (len(exons), exons)


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
            cache_directory_path=tmpdir,
        )
        genome.index()
        genes = genome.genes()
        for gene in genes:
            assert gene.id, "Gene with missing ID in %s" % (
                genome.db._load_gtf_as_dataframe(),
            )
        assert len(genes) == 2, "Expected 2 genes, got %d: %s" % (len(genes), genes)
        transcripts = genome.transcripts()
        for transcript in transcripts:
            assert transcript.id, "Transcript with missing ID in %s" % (
                genome.db._load_gtf_as_dataframe(),
            )
        assert len(transcripts) == 2, "Expected 2 transcripts, got %d: %s" % (
            len(transcripts),
            transcripts,
        )
        genes_at_locus = genome.genes_at_locus("chr1", 67092176)
        assert (
            len(genes_at_locus) == 2
        ), "Expected 2 genes at locus chr1:67092176, got %d: %s" % (
            len(genes_at_locus),
            genes_at_locus,
        )
        ids = set([gene.id for gene in genes_at_locus])
        eq_(set(["NM_001276352", "NR_075077"]), ids)

from pyensembl import Locus, Gene, ensembl_grch37, Transcript, Exon
from nose.tools import eq_


def test_Locus_string_representation():
    locus = Locus("X", 1000, 1010, "+")
    string_repr = str(locus)
    expected = "Locus(contig='X', start=1000, end=1010, strand='+')"
    eq_(string_repr, expected)


def test_Gene_string_representation():
    gene = Gene(
        gene_id="ENSG0001", gene_name="CAPITALISM",
        biotype="protein_coding", contig="Y", start=1, end=5, strand="+",
        genome=ensembl_grch37)
    string_repr = str(gene)
    expected = (
        "Gene(gene_id='ENSG0001',"
        " gene_name='CAPITALISM',"
        " biotype='protein_coding',"
        " contig='Y',"
        " start=1, end=5, strand='+', genome='GRCh37')")
    eq_(string_repr, expected)


def test_Transcript_string_representation():
    transcript = Transcript(
        transcript_id="ENST0001",
        transcript_name="CAPITALISM-001",
        gene_id="ENSG0001",
        biotype="protein_coding",
        contig="Y",
        start=1,
        end=5,
        strand="+",
        genome=ensembl_grch37)

    expected = (
        "Transcript(transcript_id='ENST0001',"
        " transcript_name='CAPITALISM-001',"
        " gene_id='ENSG0001',"
        " biotype='protein_coding',"
        " contig='Y',"
        " start=1,"
        " end=5, strand='+', genome='GRCh37')"
    )
    string_repr = str(transcript)
    eq_(string_repr, expected)


def test_Exon_string_representation():
    exon = Exon(
        exon_id="ENSE0001",
        gene_id="ENSG0001",
        gene_name="CAPITALISM",
        contig="Y",
        start=1,
        end=5,
        strand="+")

    expected = (
        "Exon(exon_id='ENSE0001',"
        " gene_id='ENSG0001',"
        " gene_name='CAPITALISM',"
        " contig='Y',"
        " start=1,"
        " end=5, strand='+')"
    )
    string_repr = str(exon)
    eq_(string_repr, expected)

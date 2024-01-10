from pyensembl import Genome
from pytest import raises
from .common import eq_
from .data import data_path

MOUSE_ENSMUSG00000017167_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.gtf"
)
MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.fa"
)
# MOUSE_ENSMUSG00000088969_NCRNA_FASTA_PATH = data_path(
#    "mouse.ensembl.81.partial.ncrna.ENSMUSG00000017167.fa")
MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.pep"
)


def no_gtf_(e):
    print("Testing for 'GTF' in %s : %s" % (type(e), e))
    assert "GTF" in str(e)


def no_transcript_(e):
    print("Testing for 'transcript' in %s : %s" % (type(e), e))
    assert "transcript" in str(e)


def no_protein_(e):
    print("Testing for 'protein' in %s : %s" % (type(e), e))
    assert "protein" in str(e)


def test_transcript_fasta_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        transcript_fasta_paths_or_urls=[MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH],
    )
    genome.index()

    eq_(2, len(genome.transcript_sequences.fasta_dictionary))

    with raises(ValueError) as e:
        genome.genes()
    no_gtf_(e)

    with raises(ValueError) as e:
        genome.gene_ids()
    no_gtf_(e)

    with raises(ValueError) as e:
        genome.gene_ids_of_gene_name("test")
    no_gtf_(e)

    with raises(ValueError) as e:
        genome.transcript_names()
    no_gtf_(e)

    with raises(ValueError) as e:
        genome.protein_sequence("test")
    no_protein_(e)


def test_protein_fasta_only():
    genome_only_proteins = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        protein_fasta_paths_or_urls=[MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH],
    )
    genome_only_proteins.index()

    eq_(4, len(genome_only_proteins.protein_sequences.fasta_dictionary))

    with raises(ValueError) as e:
        genome_only_proteins.genes()
    no_gtf_(e)

    with raises(ValueError) as e:
        genome_only_proteins.transcript_sequence("DOES_NOT_EXIST")
    no_transcript_(e)


def test_gtf_only():
    genome_only_gtf = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
    )
    genome_only_gtf.index()

    eq_(1, len(genome_only_gtf.genes()))

    with raises(ValueError) as e:
        genome_only_gtf.transcript_sequence("DOES_NOT_EXIST")

    no_transcript_(e)

    with raises(ValueError) as cm:
        genome_only_gtf.protein_sequence("genome_only_gtf")

    no_protein_(cm)


def test_gtf_transcript_only():
    genome_gtf_with_cdna = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        transcript_fasta_paths_or_urls=[MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH],
    )
    genome_gtf_with_cdna.index()

    eq_(1, len(genome_gtf_with_cdna.genes()))

    transcript = genome_gtf_with_cdna.transcripts()[0]
    assert transcript.sequence is not None

    with raises(ValueError) as e:
        transcript.protein_sequence
    no_protein_(e)


def test_gtf_protein_only():
    genome_gtf_with_proteins = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        protein_fasta_paths_or_urls=[MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH],
    )
    genome_gtf_with_proteins.index()

    eq_(1, len(genome_gtf_with_proteins.genes()))

    transcript = genome_gtf_with_proteins.transcripts()[0]
    assert transcript.protein_sequence is not None

    with raises(ValueError) as e:
        transcript.sequence
    no_transcript_(e)

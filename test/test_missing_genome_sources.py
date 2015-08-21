from __future__ import absolute_import

from pyensembl import Genome
from nose.tools import eq_, ok_, assert_raises

from .data import data_path

MOUSE_ENSMUSG00000017167_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.gtf")
MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.fa")
MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.pep")

def no_gtf_(cm):
    ok_("GTF" in str(cm.exception))

def no_transcript_(cm):
    ok_("transcript" in str(cm.exception))

def no_protein_(cm):
    ok_("protein" in str(cm.exception))

def test_transcript_fasta_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        transcript_fasta_path_or_url=MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH)
    genome.index()

    eq_(2, len(genome.transcript_sequences.fasta_dictionary))

    with assert_raises(ValueError) as cm:
        genome.genes()
    no_gtf_(cm)

    with assert_raises(ValueError) as cm:
        genome.gene_ids()
    no_gtf_(cm)

    with assert_raises(ValueError) as cm:
        genome.gene_ids_of_gene_name("test")
    no_gtf_(cm)

    with assert_raises(ValueError) as cm:
        genome.transcript_names()
    no_gtf_(cm)

    with assert_raises(ValueError) as cm:
        genome.protein_sequence("test")
    no_protein_(cm)

def test_protein_fasta_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        protein_fasta_path_or_url=MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH)
    genome.index()

    eq_(4, len(genome.protein_sequences.fasta_dictionary))

    with assert_raises(ValueError) as cm:
        genome.genes()
    no_gtf_(cm)
    with assert_raises(ValueError) as cm:
        genome.transcript_sequence("test")
    no_transcript_(cm)

def test_gtf_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH)
    genome.index()

    eq_(1, len(genome.genes()))

    with assert_raises(ValueError) as cm:
        genome.transcript_sequence("test")
    no_transcript_(cm)
    with assert_raises(ValueError) as cm:
        genome.protein_sequence("test")
    no_protein_(cm)

def test_gtf_transcript_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        transcript_fasta_path_or_url=MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH)
    genome.index()

    eq_(1, len(genome.genes()))

    transcript = genome.transcripts()[0]
    ok_(transcript.sequence)

    with assert_raises(ValueError) as cm:
        transcript.protein_sequence
    no_protein_(cm)

def test_gtf_protein_only():
    genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        protein_fasta_path_or_url=MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH)
    genome.index()

    eq_(1, len(genome.genes()))

    transcript = genome.transcripts()[0]
    ok_(transcript.protein_sequence)

    with assert_raises(ValueError) as cm:
        transcript.sequence
    no_transcript_(cm)

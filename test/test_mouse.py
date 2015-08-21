from __future__ import absolute_import

from pyensembl import Genome
from nose.tools import eq_, with_setup

from .data import data_path

MOUSE_ENSMUSG00000017167_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.gtf")
MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.fa")
MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH = data_path(
    "mouse.ensembl.81.partial.ENSMUSG00000017167.pep")


mouse_genome = None

def setup_create_genome():
    global mouse_genome
    mouse_genome = Genome(
        reference_name="GRCm38",
        annotation_name="_test_mouse_ensembl81_subset",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        transcript_fasta_path_or_url=MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH,
        protein_fasta_path_or_url=MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH)
    mouse_genome.clear_cache()
    mouse_genome.index()

@with_setup(setup=setup_create_genome)
def test_mouse_ENSMUSG00000017167():
    """
    GTF cropped from ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/
    Mus_musculus.GRCm38.81.gtf.gz via:
    grep "ENSMUSG00000017167" Mus_musculus.GRCm38.81.gtf

    Transcript FASTA cropped from ftp://ftp.ensembl.org/pub/release-81/
    fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz via:
    grep "ENSMUSG00000017167" Mus_musculus.GRCm38.cdna.all.fa -A 50

    Protein FASTA cropped from ftp://ftp.ensembl.org/pub/release-81/fasta/
    mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz via:
    grep "ENSMUSG00000017167" Mus_musculus.GRCm38.pep.all.fa -A 50

    Tested against:
    http://useast.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000017167
    """
    genes_cntnap1 = mouse_genome.genes_by_name("Cntnap1")
    eq_(len(genes_cntnap1), 1)
    gene_cntnap1 = genes_cntnap1[0]
    transcripts_cntnap1 = gene_cntnap1.transcripts
    eq_(len(transcripts_cntnap1), 2)
    transcripts_coding_cntnap1 = [
        transcript
        for transcript in transcripts_cntnap1
        if transcript.biotype == "protein_coding"
    ]
    eq_(len(transcripts_coding_cntnap1), 1)
    transcript_cntnap1 = transcripts_coding_cntnap1[0]
    eq_(transcript_cntnap1.sequence[:120],
        ("GAGAGAAGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"
         "GAGAGAGAGAGATTGGGGGTAGGAGAGAGGGAAGGGTGGATAAGGACGGAAAAAAGCTTT"))
    eq_(transcript_cntnap1.protein_sequence[:120],
        ("MMSLRLFSILLATVVSGAWGWGYYGCNEELVGPLYARSLGASSYYGLFTTARFARLHGIS"
         "GWSPRIGDPNPWLQIDLMKKHRIRAVATQGAFNSWDWVTRYMLLYGDRVDSWTPFYQKGH"))

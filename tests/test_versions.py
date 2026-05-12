"""
Issue #242: expose Ensembl annotation versions on Gene, Transcript, Exon,
and Protein (via Transcript.protein). Each entity has both a prefixed
canonical name (``gene_version``, ``versioned_gene_id`` …) and a generic
alias (``version``, ``versioned_id``).
"""

from pyensembl import Genome, Protein

from .common import eq_, ok_
from .data import (
    MOUSE_ENSMUSG00000017167_PATH,
    MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH,
    MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH,
    data_path,
)


# Versions confirmed against the GTF in tests/data/ and the source
# Ensembl 81 release: gene v6, the protein-coding transcript v3 (its protein
# ENSMUSP00000099398 also v3) and the noncoding transcript v1. All exons v1.
ENSEMBL_MOUSE_GTF = MOUSE_ENSMUSG00000017167_PATH
ENSEMBL_MOUSE_TRANSCRIPT_FASTA = MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH
ENSEMBL_MOUSE_PROTEIN_FASTA = MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH

TAIR10_GTF = data_path("arabidopsis.tair10.partial.gtf")
TAIR10_CDNA_FASTA = data_path("arabidopsis.tair10.partial.cdna.fa")


mouse_genome = Genome(
    reference_name="GRCm38",
    annotation_name="_test_versions_mouse_ensembl81_subset",
    gtf_path_or_url=ENSEMBL_MOUSE_GTF,
    transcript_fasta_paths_or_urls=[ENSEMBL_MOUSE_TRANSCRIPT_FASTA],
    protein_fasta_paths_or_urls=[ENSEMBL_MOUSE_PROTEIN_FASTA],
)

tair_genome = Genome(
    reference_name="TAIR10",
    annotation_name="_test_versions_arabidopsis_tair10_subset",
    gtf_path_or_url=TAIR10_GTF,
    transcript_fasta_paths_or_urls=[TAIR10_CDNA_FASTA],
)


def setup_module(module):
    mouse_genome.clear_cache()
    mouse_genome.index()
    tair_genome.clear_cache()
    tair_genome.index()


def test_gene_version_and_versioned_id():
    gene = mouse_genome.gene_by_id("ENSMUSG00000017167")
    eq_(gene.gene_version, 6)
    eq_(gene.version, 6)
    eq_(gene.versioned_gene_id, "ENSMUSG00000017167.6")
    eq_(gene.versioned_id, "ENSMUSG00000017167.6")


def test_transcript_version_and_versioned_id():
    transcript = mouse_genome.transcript_by_id("ENSMUST00000103109")
    eq_(transcript.transcript_version, 3)
    eq_(transcript.version, 3)
    eq_(transcript.versioned_transcript_id, "ENSMUST00000103109.3")
    eq_(transcript.versioned_id, "ENSMUST00000103109.3")


def test_transcript_protein_view():
    transcript = mouse_genome.transcript_by_id("ENSMUST00000103109")
    protein = transcript.protein
    ok_(isinstance(protein, Protein))
    eq_(protein.protein_id, "ENSMUSP00000099398")
    eq_(protein.id, "ENSMUSP00000099398")
    eq_(protein.protein_version, 3)
    eq_(protein.version, 3)
    eq_(protein.versioned_protein_id, "ENSMUSP00000099398.3")
    eq_(protein.versioned_id, "ENSMUSP00000099398.3")


def test_noncoding_transcript_has_no_protein():
    transcript = mouse_genome.transcript_by_id("ENSMUST00000138942")
    eq_(transcript.transcript_version, 1)
    eq_(transcript.protein_id, None)
    eq_(transcript.protein, None)


def test_exon_version_and_versioned_id():
    # ENSMUSE00000243064 is exon 1 of ENSMUST00000103109 at version 2
    # (other exons in this transcript are at version 1).
    exon = mouse_genome.exon_by_id("ENSMUSE00000243064")
    eq_(exon.exon_version, 2)
    eq_(exon.version, 2)
    eq_(exon.versioned_exon_id, "ENSMUSE00000243064.2")
    eq_(exon.versioned_id, "ENSMUSE00000243064.2")


def test_versions_are_none_when_gtf_lacks_them():
    # TAIR10 partial GTF has no *_version attributes; all version
    # accessors should return None and versioned_id should equal id.
    transcript = tair_genome.transcript_by_id("AT1G01010.1")
    eq_(transcript.transcript_version, None)
    eq_(transcript.version, None)
    eq_(transcript.versioned_transcript_id, "AT1G01010.1")
    eq_(transcript.versioned_id, "AT1G01010.1")

    gene = transcript.gene
    eq_(gene.gene_version, None)
    eq_(gene.version, None)
    eq_(gene.versioned_gene_id, gene.gene_id)
    eq_(gene.versioned_id, gene.gene_id)

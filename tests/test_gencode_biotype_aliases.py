"""
Issue #335 part 1: pyensembl now passes gtfparse's GENCODE_BIOTYPE_ALIASES
into read_gtf, so a vanilla GENCODE GTF (which carries gene_type /
transcript_type) gets normalized onto the Ensembl gene_biotype /
transcript_biotype column names at parse time. Without this, GENCODE
genomes silently looked like every transcript had no biotype and
Transcript.is_protein_coding returned False for everything (see #335
issue body).
"""

import os

from pyensembl import Genome

from .common import TemporaryDirectory, eq_


GENCODE_STYLE_GTF = """\
1\tHAVANA\tgene\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000000001.4"; gene_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\ttranscript\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\texon\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; exon_number "1"; exon_id "ENSETEST00000000001.2"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\tCDS\t100\t798\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; exon_number "1"; exon_id "ENSETEST00000000001.2"; protein_id "ENSPTEST00000000001.3"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\tstart_codon\t100\t102\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\tstop_codon\t799\t801\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "FAKE1";
1\tHAVANA\tgene\t2000\t3000\t.\t+\t.\tgene_id "ENSGTEST00000000002.1"; gene_type "lncRNA"; gene_name "FAKE2";
1\tHAVANA\ttranscript\t2000\t3000\t.\t+\t.\tgene_id "ENSGTEST00000000002.1"; transcript_id "ENSTTEST00000000002.1"; gene_type "lncRNA"; transcript_type "lncRNA"; gene_name "FAKE2";
1\tHAVANA\texon\t2000\t3000\t.\t+\t.\tgene_id "ENSGTEST00000000002.1"; transcript_id "ENSTTEST00000000002.1"; exon_number "1"; exon_id "ENSETEST00000000002.1"; gene_type "lncRNA"; transcript_type "lncRNA"; gene_name "FAKE2";
"""


def _make_genome(tmpdir):
    gtf_path = os.path.join(tmpdir, "gencode_style.gtf")
    with open(gtf_path, "w") as f:
        f.write(GENCODE_STYLE_GTF)
    return Genome(
        reference_name="GRCh38",
        annotation_name="_test_gencode_biotype_aliases_335p1",
        gtf_path_or_url=gtf_path,
        cache_directory_path=tmpdir,
    )


def test_gencode_gene_type_normalized_to_gene_biotype():
    """A GENCODE GTF (with gene_type / transcript_type attributes) is
    parsed so that the sqlite database stores gene_biotype /
    transcript_biotype, and querying by Ensembl biotype names works."""
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir)
        genome.index()

        # Biotype columns exist in the gene + transcript tables
        assert genome.db.column_exists("gene", "gene_biotype")
        assert genome.db.column_exists("transcript", "transcript_biotype")

        # Filter by biotype works against the renamed columns
        coding_genes = genome.genes(biotype="protein_coding")
        eq_(len(coding_genes), 1)
        eq_(coding_genes[0].id, "ENSGTEST00000000001.4")
        eq_(coding_genes[0].biotype, "protein_coding")

        ncrna_genes = genome.genes(biotype="lncRNA")
        eq_(len(ncrna_genes), 1)
        eq_(ncrna_genes[0].id, "ENSGTEST00000000002.1")


def test_gencode_transcript_type_normalized_to_transcript_biotype():
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir)
        genome.index()

        coding = genome.transcripts(biotype="protein_coding")
        eq_(len(coding), 1)
        eq_(coding[0].id, "ENSTTEST00000000001.5")
        eq_(coding[0].biotype, "protein_coding")


def test_gencode_transcript_is_protein_coding_true_after_alias():
    """The original #335 repro: Variant(...).effects() returned
    NoncodingTranscript for every transcript because is_protein_coding
    was False. With attribute_aliases the biotype column lands on the
    Ensembl name and is_protein_coding returns True for coding rows."""
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir)
        genome.index()
        coding = genome.transcript_by_id("ENSTTEST00000000001.5")
        assert coding.is_protein_coding is True
        noncoding = genome.transcript_by_id("ENSTTEST00000000002.1")
        assert noncoding.is_protein_coding is False


def test_pure_ensembl_gtf_unaffected_by_aliasing():
    """Vanilla Ensembl GTFs already use gene_biotype / transcript_biotype.
    Passing the alias dict must be a no-op for them (gtfparse skips
    missing alias sources). Reuses the existing mouse partial fixture
    which is Ensembl 81 style."""
    from .data import (
        custom_mouse_genome_grcm38_subset,
        setup_init_custom_mouse_genome,
    )
    setup_init_custom_mouse_genome()
    coding = custom_mouse_genome_grcm38_subset.transcripts(biotype="protein_coding")
    eq_(len(coding), 1)
    eq_(coding[0].id, "ENSMUST00000103109")
    assert coding[0].is_protein_coding is True

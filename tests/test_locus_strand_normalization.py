"""Regression tests for strand normalization in locus queries."""

import pytest

from pyensembl import Genome


@pytest.fixture
def stranded_genome(tmp_path):
    gtf = tmp_path / "stranded.gtf"
    gtf.write_text(
        '17\ttest\tgene\t100\t200\t.\t-\t.\tgene_id "G_NEG"; gene_name "NEG_GENE";\n'
        '17\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "G_POS"; gene_name "POS_GENE";\n'
        '17\ttest\ttranscript\t100\t200\t.\t-\t.\tgene_id "G_NEG"; transcript_id "T_NEG";\n'
        '17\ttest\ttranscript\t100\t200\t.\t+\t.\tgene_id "G_POS"; transcript_id "T_POS";\n'
    )
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="_test_locus_strand_normalization",
        gtf_path_or_url=str(gtf),
        cache_directory_path=str(tmp_path / "cache"),
    )
    genome.index()
    try:
        yield genome
    finally:
        genome.close()


@pytest.mark.parametrize("strand", ["+", "+1", "1", 1])
def test_locus_queries_normalize_forward_strand(stranded_genome, strand):
    assert [g.id for g in stranded_genome.genes_at_locus("17", 150, strand=strand)] == ["G_POS"]
    assert stranded_genome.transcript_ids_at_locus("17", 150, strand=strand) == ["T_POS"]
    assert stranded_genome.db.column_values_at_locus(
        "gene_id", "gene", "17", 150, strand=strand
    ) == ["G_POS"]


@pytest.mark.parametrize("strand", ["-", "-1", -1])
def test_locus_queries_normalize_reverse_strand(stranded_genome, strand):
    assert [g.id for g in stranded_genome.genes_at_locus("17", 150, strand=strand)] == ["G_NEG"]
    assert stranded_genome.transcript_ids_at_locus("17", 150, strand=strand) == ["T_NEG"]
    assert stranded_genome.db.column_values_at_locus(
        "gene_id", "gene", "17", 150, strand=strand
    ) == ["G_NEG"]


@pytest.mark.parametrize("invalid_strand", ["invalid", "forward", "", 0])
def test_locus_queries_reject_invalid_strand(stranded_genome, invalid_strand):
    with pytest.raises(ValueError, match="Invalid strand"):
        stranded_genome.genes_at_locus("17", 150, strand=invalid_strand)

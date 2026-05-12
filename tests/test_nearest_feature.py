"""
Issue #283: surface find_nearest_locus as Genome.nearest_gene /
Genome.nearest_transcript so callers can locate the closest annotated
feature to a position even when no feature overlaps it directly.
"""

from pyensembl import Gene, Transcript

from .common import eq_
from .data import (
    custom_mouse_genome_grcm38_subset,
    setup_init_custom_mouse_genome,
)


def setup_module(module):
    setup_init_custom_mouse_genome()


def test_nearest_gene_overlapping_returns_distance_zero():
    # ENSMUSG00000017167 spans chr11:101170523-101190724 on +; query
    # inside the gene must return distance 0 and that gene.
    distance, gene = custom_mouse_genome_grcm38_subset.nearest_gene(
        contig="11", position=101180000
    )
    eq_(distance, 0)
    assert isinstance(gene, Gene)
    eq_(gene.id, "ENSMUSG00000017167")


def test_nearest_gene_upstream_position_reports_gap():
    # 1000 bases upstream of gene start (101170523)
    distance, gene = custom_mouse_genome_grcm38_subset.nearest_gene(
        contig="11", position=101169523
    )
    eq_(distance, 1000)
    eq_(gene.id, "ENSMUSG00000017167")


def test_nearest_gene_downstream_interval_reports_gap():
    # interval entirely downstream of gene end (101190724); distance
    # is gap from end of interval to start of next feature backward,
    # i.e. position 101200724 - gene.end 101190724 = 10000.
    distance, gene = custom_mouse_genome_grcm38_subset.nearest_gene(
        contig="11", position=101200724, end=101201000
    )
    eq_(distance, 10000)
    eq_(gene.id, "ENSMUSG00000017167")


def test_nearest_gene_on_empty_contig_returns_none():
    distance, gene = custom_mouse_genome_grcm38_subset.nearest_gene(
        contig="22", position=100
    )
    eq_(distance, float("inf"))
    eq_(gene, None)


def test_nearest_transcript_overlapping_returns_distance_zero():
    distance, transcript = (
        custom_mouse_genome_grcm38_subset.nearest_transcript(
            contig="11", position=101180000
        )
    )
    eq_(distance, 0)
    assert isinstance(transcript, Transcript)
    # both transcripts overlap; either is acceptable, just confirm we got one
    assert transcript.id in {"ENSMUST00000103109", "ENSMUST00000138942"}


def test_nearest_transcript_upstream_picks_closer_endpoint():
    # ENSMUST00000138942 starts at 101170523, ENSMUST00000103109 at 101176041;
    # query at 101170000 is 523 bp from the closer transcript.
    distance, transcript = (
        custom_mouse_genome_grcm38_subset.nearest_transcript(
            contig="11", position=101170000
        )
    )
    eq_(distance, 523)
    eq_(transcript.id, "ENSMUST00000138942")


def test_nearest_gene_strand_filter():
    # Restricting to the wrong strand (the only gene is on +) should
    # produce no candidates.
    distance, gene = custom_mouse_genome_grcm38_subset.nearest_gene(
        contig="11", position=101180000, strand="-"
    )
    eq_(distance, float("inf"))
    eq_(gene, None)

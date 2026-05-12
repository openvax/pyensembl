"""
Issue #177: biotype filter pushed into the SQL on Genome.genes(),
Genome.transcripts(), and the underlying *_ids queries.
"""

import pytest

from .common import eq_
from .data import (
    custom_mouse_genome_grcm38_subset,
    setup_init_custom_mouse_genome,
)


def setup_module(module):
    setup_init_custom_mouse_genome()


def test_transcripts_filter_by_biotype():
    all_transcripts = custom_mouse_genome_grcm38_subset.transcripts()
    eq_(len(all_transcripts), 2)
    biotypes = {t.biotype for t in all_transcripts}
    # mouse partial fixture has exactly one protein_coding transcript and
    # one processed_transcript transcript for Cntnap1
    assert "protein_coding" in biotypes
    assert biotypes != {"protein_coding"}

    coding = custom_mouse_genome_grcm38_subset.transcripts(
        biotype="protein_coding"
    )
    eq_(len(coding), 1)
    eq_(coding[0].biotype, "protein_coding")


def test_transcript_ids_filter_by_biotype():
    coding_ids = custom_mouse_genome_grcm38_subset.transcript_ids(
        biotype="protein_coding"
    )
    eq_(coding_ids, ["ENSMUST00000103109"])


def test_genes_filter_by_biotype():
    coding_genes = custom_mouse_genome_grcm38_subset.genes(
        biotype="protein_coding"
    )
    eq_(len(coding_genes), 1)
    eq_(coding_genes[0].biotype, "protein_coding")
    # Filtering by a biotype that's absent should return an empty list, not raise.
    eq_(
        custom_mouse_genome_grcm38_subset.genes(biotype="miRNA"),
        [],
    )


def test_biotype_combines_with_contig_and_strand():
    # protein_coding transcript ENSMUST00000103109 is on chr11, strand +
    transcripts = custom_mouse_genome_grcm38_subset.transcripts(
        contig="11", strand="+", biotype="protein_coding"
    )
    eq_(len(transcripts), 1)
    eq_(transcripts[0].id, "ENSMUST00000103109")
    # wrong strand should yield zero results
    eq_(
        custom_mouse_genome_grcm38_subset.transcripts(
            contig="11", strand="-", biotype="protein_coding"
        ),
        [],
    )


def test_biotype_on_feature_without_biotype_column_raises():
    # The 'exon' table has no exon_biotype column; if a user were to
    # call _all_feature_values with biotype= for it we should fail loudly
    # rather than silently return everything.
    with pytest.raises(ValueError, match="biotype"):
        custom_mouse_genome_grcm38_subset._all_feature_values(
            column="exon_id", feature="exon", biotype="protein_coding"
        )

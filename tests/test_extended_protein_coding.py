"""
Issue #169: three-tier protein-coding biotype ontology.

* ``is_protein_coding`` — strict, the canonical ``"protein_coding"`` only.
  Unchanged from previous releases for back-compat with downstream effect
  predictors (varcode).
* ``is_protein_coding_extended`` — widens to IG/TR gene segments and
  translated pseudogenes.  Excludes NMD/NSD because those are degraded.
* ``is_translated`` — widens further to include NMD/NSD; covers any
  biotype that goes through the ribosome at all.
"""

from pyensembl.locus_with_genome import (
    EXTENDED_PROTEIN_CODING_BIOTYPES,
    LocusWithGenome,
    PROTEIN_CODING_BIOTYPES,
    TRANSLATED_BIOTYPES,
)

from .common import eq_
from .data import (
    custom_mouse_genome_grcm38_subset,
    setup_init_custom_mouse_genome,
)


class _BiotypeFixture:
    """Stand-in that exposes only ``.biotype`` so we can invoke the three
    LocusWithGenome property accessors without constructing a Genome."""

    is_protein_coding = LocusWithGenome.is_protein_coding
    is_protein_coding_extended = LocusWithGenome.is_protein_coding_extended
    is_translated = LocusWithGenome.is_translated

    def __init__(self, biotype):
        self.biotype = biotype


# -----------------------------
# Set constants
# -----------------------------


def test_biotype_set_ontology_is_monotonic():
    """strict ⊂ extended ⊂ translated."""
    assert PROTEIN_CODING_BIOTYPES.issubset(EXTENDED_PROTEIN_CODING_BIOTYPES)
    assert EXTENDED_PROTEIN_CODING_BIOTYPES.issubset(TRANSLATED_BIOTYPES)


def test_extended_set_adds_immunoglobulin_and_tcr_and_translated_pseudogenes():
    extras = EXTENDED_PROTEIN_CODING_BIOTYPES - PROTEIN_CODING_BIOTYPES
    eq_(
        extras,
        {
            "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
            "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
            "polymorphic_pseudogene",
            "translated_processed_pseudogene",
            "translated_unprocessed_pseudogene",
        },
    )


def test_translated_set_only_adds_nmd_and_nsd_over_extended():
    extras = TRANSLATED_BIOTYPES - EXTENDED_PROTEIN_CODING_BIOTYPES
    eq_(extras, {"nonsense_mediated_decay", "non_stop_decay"})


def test_strict_excludes_everything_but_canonical():
    eq_(PROTEIN_CODING_BIOTYPES, {"protein_coding"})


# -----------------------------
# Property dispatch
# -----------------------------


def test_canonical_protein_coding_is_true_in_every_tier():
    f = _BiotypeFixture("protein_coding")
    assert f.is_protein_coding
    assert f.is_protein_coding_extended
    assert f.is_translated


def test_ig_gene_segments_are_extended_and_translated_but_not_strict():
    for biotype in ("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene"):
        f = _BiotypeFixture(biotype)
        assert not f.is_protein_coding, biotype
        assert f.is_protein_coding_extended, biotype
        assert f.is_translated, biotype


def test_tcr_gene_segments_are_extended_and_translated_but_not_strict():
    for biotype in ("TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"):
        f = _BiotypeFixture(biotype)
        assert not f.is_protein_coding, biotype
        assert f.is_protein_coding_extended, biotype
        assert f.is_translated, biotype


def test_translated_pseudogenes_are_extended_and_translated_but_not_strict():
    for biotype in (
        "polymorphic_pseudogene",
        "translated_processed_pseudogene",
        "translated_unprocessed_pseudogene",
    ):
        f = _BiotypeFixture(biotype)
        assert not f.is_protein_coding, biotype
        assert f.is_protein_coding_extended, biotype
        assert f.is_translated, biotype


def test_nmd_and_nsd_are_translated_only():
    """NMD/NSD biotypes are translated by the ribosome but don't yield a
    stable protein product, so they're True for `is_translated` and
    False for the other two tiers."""
    for biotype in ("nonsense_mediated_decay", "non_stop_decay"):
        f = _BiotypeFixture(biotype)
        assert not f.is_protein_coding, biotype
        assert not f.is_protein_coding_extended, biotype
        assert f.is_translated, biotype


def test_noncoding_biotypes_are_false_at_every_tier():
    """Sanity: typical non-coding biotypes (lncRNA, miRNA, snRNA, plain
    pseudogenes, regulatory) are False across all three tiers."""
    for biotype in (
        "lncRNA",
        "miRNA",
        "snRNA",
        "snoRNA",
        "rRNA",
        "processed_pseudogene",
        "unprocessed_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "transcribed_unitary_pseudogene",
        "TEC",
        "retained_intron",
        "processed_transcript",
        "antisense",
    ):
        f = _BiotypeFixture(biotype)
        assert not f.is_protein_coding, biotype
        assert not f.is_protein_coding_extended, biotype
        assert not f.is_translated, biotype


def test_none_biotype_is_false_at_every_tier():
    """Some GTFs (older Ensembl releases, TAIR) don't carry biotype
    attributes at all. `biotype=None` must not crash and must be False
    across the board."""
    f = _BiotypeFixture(None)
    assert not f.is_protein_coding
    assert not f.is_protein_coding_extended
    assert not f.is_translated


# -----------------------------
# Backward compat with existing callers
# -----------------------------


def setup_module(module):
    setup_init_custom_mouse_genome()


def test_existing_mouse_protein_coding_transcript_is_strict_protein_coding():
    """The mouse partial fixture has exactly one transcript with
    biotype=='protein_coding' and one with biotype=='processed_transcript'.
    is_protein_coding must continue to return True for the first and
    False for the second; downstream effect prediction depends on this."""
    coding = custom_mouse_genome_grcm38_subset.transcripts(biotype="protein_coding")
    eq_(len(coding), 1)
    assert coding[0].is_protein_coding
    assert coding[0].is_protein_coding_extended
    assert coding[0].is_translated

    all_transcripts = custom_mouse_genome_grcm38_subset.transcripts()
    noncoding = [t for t in all_transcripts if t.biotype != "protein_coding"]
    eq_(len(noncoding), 1)
    assert not noncoding[0].is_protein_coding
    assert not noncoding[0].is_protein_coding_extended
    assert not noncoding[0].is_translated

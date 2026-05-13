# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from .locus import Locus


# Reference: http://www.gencodegenes.org/pages/biotypes.html
#
# Three-tier ontology of "does this transcript make a polypeptide?":
#
#   strict  ⊂  extended  ⊂  translated
#
# Each tier widens the previous one.

# Tier 1 — strict canonical protein-coding biotype only. This is what
# Ensembl uses for the bulk of reference proteome work; varcode and other
# downstream effect predictors expect this set to drive `is_protein_coding`.
PROTEIN_CODING_BIOTYPES = frozenset({"protein_coding"})

# Tier 2 — anything Ensembl/GENCODE marks as producing a stable, functional
# polypeptide:
#  * IG_{C,D,J,V}_gene  - immunoglobulin gene segments. They produce
#    protein after V(D)J recombination; the biotype reflects the productive
#    pre-recombination form.
#  * TR_{C,D,J,V}_gene  - T-cell receptor gene segments, same situation.
#  * polymorphic_pseudogene  - codes for protein in some individuals,
#    pseudogene in others. Genuinely protein-producing where it's coded.
#  * translated_processed_pseudogene / translated_unprocessed_pseudogene -
#    pseudogenes with ribosome occupancy / mass-spec evidence of being
#    translated into a stable product.
#
# Excludes NMD / NSD because those are translated transiently and the
# product is targeted for degradation, not stable accumulation.
EXTENDED_PROTEIN_CODING_BIOTYPES = frozenset(PROTEIN_CODING_BIOTYPES) | frozenset({
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
    "polymorphic_pseudogene",
    "translated_processed_pseudogene",
    "translated_unprocessed_pseudogene",
})

# Tier 3 — every biotype that gets translated on a ribosome, regardless of
# whether the product survives. Adds nonsense_mediated_decay and
# non_stop_decay to the extended set.
#
# Useful when you care about whether a variant lands in a translated frame
# (e.g. picking a top variant effect, RNA-seq peptide analysis) rather
# than whether it perturbs a stably-expressed protein product.
TRANSLATED_BIOTYPES = frozenset(EXTENDED_PROTEIN_CODING_BIOTYPES) | frozenset({
    "nonsense_mediated_decay",
    "non_stop_decay",
})


class LocusWithGenome(Locus):
    """
    Common base class for Gene and Transcript to avoid copying
    their shared logic.
    """

    def __init__(self, contig, start, end, strand, biotype, genome):
        Locus.__init__(self, contig, start, end, strand)
        self.genome = genome
        self.db = self.genome.db
        self.biotype = biotype

    def to_dict(self):
        return dict(
            contig=self.contig,
            start=self.start,
            end=self.end,
            strand=self.strand,
            biotype=self.biotype,
            genome=self.genome,
        )

    @property
    def is_protein_coding(self):
        """
        True iff this entry's biotype is the canonical ``"protein_coding"``.
        Conservative by design - this is what downstream effect predictors
        (e.g. varcode) read to decide whether a variant lands in a
        translatable transcript. See :attr:`is_protein_coding_extended`
        for a wider definition that also covers IG/TR gene segments and
        translated pseudogenes, and :attr:`is_translated` for the widest
        definition that additionally includes NMD/NSD targets.
        """
        return self.biotype in PROTEIN_CODING_BIOTYPES

    @property
    def is_protein_coding_extended(self):
        """
        True for any biotype that makes a stable, functional protein
        product:

        * ``protein_coding`` (canonical)
        * ``IG_{C,D,J,V}_gene`` and ``TR_{C,D,J,V}_gene`` (immunoglobulin
          and T-cell receptor gene segments — produce protein after
          V(D)J recombination)
        * ``polymorphic_pseudogene`` (codes in some individuals)
        * ``translated_{processed,unprocessed}_pseudogene`` (pseudogenes
          with translation evidence)

        Excludes ``nonsense_mediated_decay`` and ``non_stop_decay`` —
        those biotypes are translated but the product is targeted for
        degradation rather than stable accumulation. Use
        :attr:`is_translated` if you want those.
        """
        return self.biotype in EXTENDED_PROTEIN_CODING_BIOTYPES

    @property
    def is_translated(self):
        """
        True for any biotype that gets translated on a ribosome, even
        transiently. Equivalent to :attr:`is_protein_coding_extended`
        plus ``nonsense_mediated_decay`` and ``non_stop_decay``.

        Useful when you care about whether a variant lands in a
        translated frame at all (e.g. picking a top variant effect,
        RNA-seq peptide analysis) rather than whether the encoded
        protein product is stably expressed.
        """
        return self.biotype in TRANSLATED_BIOTYPES

# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

"""
Ensembl/GENCODE biotype classifications, for more information:
http://useast.ensembl.org/Help/Faq?id=468

Definitions for GENCODE biotypes from:
http://www.gencodegenes.org/gencode_biotypes.html

(not all of these are necessarily used in Ensembl)

IG_C_gene
IG_D_gene
IG_J_gene
IG_V_gene
TR_C_gene
TR_J_gene
TR_V_gene
TR_D_gene
----------------
Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes imported
or annotated according to the IMGT. IG_C_pseudogene

IG_J_pseudogene
IG_V_pseudogene
TR_V_pseudogene
TR_J_pseudogene
----------------
Inactivated immunoglobulin gene.

Mt_rRNA
Mt_tRNA
miRNA
misc_RNA
rRNA
snRNA
snoRNA
----------------
Non-coding RNA predicted using sequences from RFAM and miRBase

Mt_tRNA_pseudogene
tRNA_pseudogene
snoRNA_pseudogene
snRNA_pseudogene
scRNA_pseudogene
rRNA_pseudogene
misc_RNA_pseudogene
miRNA_pseudogene
Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline

TEC
----------------
To be Experimentally Confirmed. This is used for non-spliced EST clusters that
have polyA features. This category has been specifically created for the ENCODE
project to highlight regions that could indicate the presence of protein coding
genes that require experimental validation, either by 5' RACE or RT-PCR to
extend the transcripts, or by confirming expression of the putatively-encoded
peptide with specific antibodies.

nonsense_mediated_decay
----------------
If the coding sequence (following the appropriate reference) of a transcript
finishes >50bp from a downstream splice site then it is tagged as NMD.
If the variant does not cover the full reference coding sequence then it is
annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure
of the missing portion is the transcript will be subject to NMD.

non_stop_decay
----------------
Transcripts that have polyA features (including signal) without a prior stop
codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS
without 3' UTR. These transcripts are subject to degradation.


retained_intron
----------------
Alternatively spliced transcript believed to contain intronic sequence relative
to other, coding, variants. protein_coding Contains an open reading frame (ORF).

processed_transcript
----------------
Doesn't contain an ORF.

non_coding
----------------
Transcript which is known from the literature to not be protein coding.

ambiguous_orf
----------------
Transcript believed to be protein coding, but with more than one possible
open reading frame.

sense_intronic
----------------
Long non-coding transcript in introns of a coding gene that does not overlap any
exons.

sense_overlapping
----------------
Long non-coding transcript that contains a coding gene in its intron on the same
strand.

antisense
----------------
Has transcripts that overlap the genomic span (i.e. exon or introns) of a
protein-coding locus on the opposite strand.

known_ncrna
----------------

pseudogene
----------------
Have homology to proteins but generally suffer from a disrupted coding sequence
and an active homologous gene can be found at another locus. Sometimes these
entries have an intact coding sequence or an open but truncated ORF, in which
case there is other evidence used (for example genomic polyA stretches at
the 3' end) to classify them as a pseudogene. Can be further classified as one
of the following.

processed_pseudogene
----------------
Pseudogene that lack introns and is thought to arise from reverse transcription
of mRNA followed by reinsertion of DNA into the genome.

polymorphic_pseudogene
----------------
Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the
gene is translated.

retrotransposed
----------------
Pseudogene owing to a reverse transcribed and re-inserted sequence.

transcribed_processed_pseudogene
transcribed_unprocessed_pseudogene
transcribed_unitary_pseudogene
----------------
Pseudogene where protein homology or genomic structure indicates a pseudogene,
but the presence of locus-specific transcripts indicates expression.

translated_unprocessed_pseudogene
----------------
Pseudogene that has mass spec data suggesting that it is also translated.

unitary_pseudogene
----------------
A species specific unprocessed pseudogene without a parent gene, as it has an
active orthologue in another species.

unprocessed_pseudogene
----------------
Pseudogene that can contain introns since produced by gene duplication.

artifact
----------------
Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)

lincRNA
----------------
Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily
conserved, intergenic regions.

LRG_gene
----------------
Gene in a "Locus Reference Genomic" region known to have disease-related
sequence variations.

3prime_overlapping_ncrna
----------------
Transcripts where ditag and/or published experimental data strongly supports the
existence of short non-coding transcripts transcribed from the 3'UTR.

disrupted_domain
----------------
Otherwise viable coding region omitted from this alternatively spliced
transcript because the splice variation affects a region coding for a
protein domain.
"""

TCR_biotypes = {
    'TR_C_gene',
    'TR_D_gene',
    'TR_gene',
    'TR_J_gene',
    'TR_V_gene'
}

IG_biotypes = {
    'IG_C_gene',
    'IG_D_gene',
    'IG_gene',
    'IG_J_gene',
    'IG_LV_gene',
    'IG_M_gene',
    'IG_V_gene',
    'IG_Z_gene',
}

non_immune_protein_coding = {
    # premature stop codon will cause degradation
    'nonsense_mediated_decay',
    # TODO: find out what this is
    'nontranslating_CDS',
    # no stop codon within transcript, should get degraded at translation
    'non_stop_decay',
    # ordinary protein coding transcript
    'protein_coding',
    # usually a non-coding pseudogene but can be translated in some individuals
    # depending on common genetic variation
    'polymorphic_pseudogene',
    # Otherwise viable coding region omitted from this alternatively spliced
    # transcript because the splice variation affects a region coding for a
    # protein domain.
    'disrupted_domain',
    # Gene in a "Locus Reference Genomic" region known to have disease-related
    # sequence variations.
    'LRG_gene',
}

protein_coding = set.union(
    TCR_biotypes,
    IG_biotypes,
    non_immune_protein_coding
)


coding_pseudogenes = {
    # TODO: why is disrupted_domain a pseudogene rather than
    # a translated protein we expect to be dysfunctional?
    'disrupted_domain',
    'IG_C_pseudogene',
    'IG_J_pseudogene',
    'IG_pseudogene',
    'IG_V_pseudogene',
    # Found in ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus
    'IG_D_pseudogene',
    'processed_pseudogene',
    'pseudogene',
    'transcribed_processed_pseudogene',
    'transcribed_unitary_pseudogene',
    'transcribed_unprocessed_pseudogene',
    'translated_processed_pseudogene',
    'TR_J_pseudogene',
    'TR_pseudogene',
    'TR_V_pseudogene',
    'unitary_pseudogene',
    'unprocessed_pseudogene',
    # to be experimentally confirmed
    'TEC',
    # TODO: should this be here or considered protein_coding?
    'translated_unprocessed_pseudogene',
    # pseudogene owing to a reverse transcribed and re-inserted sequence.
    "retrotransposed",
}

long_noncoding = {
    'lincRNA',
    'ncrna_host',
    '3prime_overlapping_ncrna',
    # why is an ambiguous ORF noncoding? Isn't rather coding but we can't
    # yet determine the sequence?
    'ambiguous_orf',
    'antisense',
    'antisense_RNA',
    'non_coding',
    'processed_transcript',
    'retained_intron',
    'sense_intronic',
    'sense_overlapping',
    'known_ncrna',
    # unspliced lncRNAs that are several kb in size.
    'macro_lncRNA',
}

mitochondrial = {
    'Mt_rRNA',
    'Mt_tRNA',
    'Mt_tRNA_pseudogene',
}

# short RNAs homologous to functional but which themselves are expected to
# be inert
short_noncoding_pseudogene = {
    'miRNA_pseudogene',
    'misc_RNA_pseudogene',
    'ncRNA_pseudogene',
    'rRNA_pseudogene',
    'scRNA_pseudogene',
    'snoRNA_pseudogene',
    'snRNA_pseudogene',
    'tRNA_pseudogene'
}

# short RNAs expected to have some effect or function in the nucleus
short_noncoding_functional = {
    'miRNA',
    'misc_RNA',
    'ncRNA',
    'rRNA',
    'scRNA',
    'snlRNA',
    'snoRNA',
    'snRNA',
    'tRNA',
    'sRNA',
    # Small Cajal body-specific RNA
    'scaRNA',
    # Vault RNA (http://en.wikipedia.org/wiki/Vault_RNA)
    'vaultRNA',
    'ribozyme',
}

short_noncoding = set.union(
    short_noncoding_functional,
    short_noncoding_pseudogene,
    mitochondrial)

valid_biotypes = set.union(
    protein_coding,
    coding_pseudogenes,
    long_noncoding,
    short_noncoding,
    # used to tag mistakes in the annotation database
    {"artifact"},
)

def is_valid_biotype(biotype):
    return biotype in valid_biotypes

def is_coding_biotype(biotype):
    return biotype in protein_coding

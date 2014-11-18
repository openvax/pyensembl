
"""
Ensembl/GENCODE biotype classifications from:
http://useast.ensembl.org/Help/Faq?id=468
"""

TCR_biotype_set = {
    'TR_C_gene',
    'TR_D_gene',
    'TR_gene',
    'TR_J_gene',
    'TR_V_gene'
}

IG_biotype_set = {
    'IG_C_gene',
    'IG_D_gene',
    'IG_gene',
    'IG_J_gene',
    'IG_LV_gene',
    'IG_M_gene',
    'IG_V_gene',
    'IG_Z_gene',
}

non_immune_protein_coding_biotype_set = {
    'nonsense_mediated_decay',
    'nontranslating_CDS',
    'non_stop_decay',
    'polymorphic',
    'polymorphic_pseudogene',
    'protein_coding',
}

protein_coding_biotype_set = set.union(
    TCR_biotype_set,
    IG_biotype_set,
    non_immune_protein_coding_biotype_set
)


pseudogene_biotype_set = {
    'disrupted_domain',
    'IG_C_pseudogene',
    'IG_J_pseudogene',
    'IG_pseudogene',
    'IG_V_pseudogene',
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
    'unprocessed_pseudogene'
}

long_noncoding_biotype_set = {
    '3prime_overlapping_ncrna',
    'ambiguous_orf',
    'antisense',
    'antisense_RNA',
    'lincRNA',
    'ncrna_host',
    'non_coding',
    'processed_transcript',
    'retained_intron',
    'sense_intronic',
    'sense_overlapping'
}

short_noncoding_biotype_set = {
    'miRNA',
    'miRNA_pseudogene',
    'misc_RNA',
    'misc_RNA_pseudogene',
    'Mt_rRNA',
    'Mt_tRNA',
    'Mt_tRNA_pseudogene',
    'ncRNA',
    'ncRNA_pseudogene',
    'rRNA',
    'rRNA_pseudogene',
    'scRNA',
    'scRNA_pseudogene',
    'snlRNA',
    'snoRNA',
    'snoRNA_pseudogene',
    'snRNA',
    'snRNA_pseudogene',
    'tRNA',
    'tRNA_pseudogene'
}

all_valid_biotypes = set.union(
    protein_coding_biotype_set,
    pseudogene_biotype_set,
    long_noncoding_biotype_set,
    short_noncoding_biotype_set,
)


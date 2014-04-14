from lazy_loader import LazyLoader

GENE_HEADER = [
    'gene_id', 'biotype', 'analysis_id',
    'seq_region_id', 'seq_region_start', 'seq_region_end',
    'seq_region_strand', 'display_xref_id',
    'source', 'status', 'description', 'is_current',
    'canonical_transcript_id', 'stable_id',
    'version', 'created_date', 'modified_date'
]

loader = LazyLoader("gene.txt.gz", GENE_HEADER)


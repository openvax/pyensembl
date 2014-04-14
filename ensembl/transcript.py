from lazy_loader import LazyLoader

TRANSCRIPT_HEADER = [
    "transcript_id", "gene_id", "analysis_id",
    "seq_region_id", "seq_region_start", "seq_region_end",
    "seq_region_strand", "display_xref_id", "biotype", "status",
    "description", "is_current", "canonical_translation_id", "stable_id",
    "version", "created_date", "modified_date"
]

loader = LazyLoader("transcript.txt.gz", TRANSCRIPT_HEADER)

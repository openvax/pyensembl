from lazy_loader import LazyLoader

TRANSLATION_HEADER = [
    "translation_id", "transcript_id", "seq_start",
    "start_exon_id", "seq_end", "end_exon_id",
    "stable_id","version", "created_date", "modified_date"
]
loader = LazyLoader("translation.txt.gz", TRANSLATION_HEADER)
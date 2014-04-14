from lazy_loader import LazyLoader

EXON_HEADER = [
    "exon_id", "seq_region_id", "seq_region_start",
    "seq_region_end", "seq_region_strand", "phase", "end_phase",
    "is_current", "is_constitutive", "stable_id", "version",
    "created_date", "modified_date"
]

loader = LazyLoader("exon.txt.gz", EXON_HEADER)

EXON_TRANSCRIPT_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-RELEASE/mysql/homo_sapiens_core_RELEASE_37/exon_transcript.txt.gz"

TRANSLATION_HEADER = [
    "translation_id", "transcript_id", "seq_start",
    "start_exon_id", "seq_end", "end_exon_id",
    "stable_id","version", "created_date", "modified_date"
]

TRANSLATION_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-RELEASE/mysql/homo_sapiens_core_RELEASE_37/translation.txt.gz"

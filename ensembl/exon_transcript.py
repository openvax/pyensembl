from lazy_loader import LazyLoader 

EXON_TRANSCRIPT_HEADER = ['exon_id', 'transcript_id', 'rank']
loader = LazyLoader("exon_transcript.txt.gz", EXON_TRANSCRIPT_HEADER)



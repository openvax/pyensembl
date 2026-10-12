# 2.27.2

Transcript.protein_id, Transcript.protein and Transcript.protein_sequence now
return None when a noncoding-only or custom annotation has no CDS protein_id
column. This fixes the missing-table error reported in #459 without changing
canonical protein identities in ordinary annotations.

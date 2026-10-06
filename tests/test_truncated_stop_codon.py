"""
Transcript.complete must return False, not raise, when a stop codon is
annotated over fewer than three bases (e.g. truncated at a contig edge).
"""

import os
from tempfile import TemporaryDirectory

from pyensembl import Genome

from .common import eq_

ATTRIBUTES = (
    'gene_id "GTEST1"; transcript_id "TTEST1"; gene_name "TRUNC1"; '
    'gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
)

# Exon 100..800 on the plus strand; the stop codon row covers only 799..800.
TRUNCATED_STOP_GTF = "".join(
    "1\ttest\t%s\t%d\t%d\t.\t+\t%s\t%s%s\n" % (feature, start, end, frame, ATTRIBUTES, extra)
    for feature, start, end, frame, extra in [
        ("transcript", 100, 800, ".", ""),
        ("exon", 100, 800, ".", ' exon_number "1"; exon_id "ETEST1";'),
        ("CDS", 100, 798, "0", ' exon_number "1"; protein_id "PTEST1";'),
        ("start_codon", 100, 102, "0", ' protein_id "PTEST1";'),
        ("stop_codon", 799, 800, "0", ' protein_id "PTEST1";'),
    ]
)


def test_truncated_stop_codon_is_not_complete():
    with TemporaryDirectory() as tmpdir:
        gtf_path = os.path.join(tmpdir, "truncated_stop.gtf")
        with open(gtf_path, "w") as f:
            f.write(TRUNCATED_STOP_GTF)
        transcript_fasta = os.path.join(tmpdir, "truncated_stop.cdna.fa")
        with open(transcript_fasta, "w") as f:
            f.write(">TTEST1\n" + "ATG" + "A" * 698 + "\n")
        genome = Genome(
            reference_name="GRCh38",
            annotation_name="_test_truncated_stop_codon",
            gtf_path_or_url=gtf_path,
            transcript_fasta_paths_or_urls=[transcript_fasta],
            cache_directory_path=tmpdir,
        )
        genome.index()
        transcript = genome.transcript_by_id("TTEST1")
        eq_(transcript.contains_stop_codon, True)
        eq_(transcript.stop_codon_complete, False)
        eq_(transcript.complete, False)
        genome.close()

"""
Regression tests for Transcript.complete and Transcript.coding_sequence on
Ensembl Plants / TAIR10-style data.

Covers issue #252 (Transcript.complete returning the wrong result for TAIR
transcripts — fixed by #268 which stopped stripping the `.N` isoform suffix
from non-ENSEMBL FASTA headers) and a related defect where coding_sequence
raised ValueError for transcripts missing a start_codon or stop_codon
feature.

GTF and cDNA FASTA fragments were taken from Ensembl Plants release 57:
    http://ftp.ensemblgenomes.org/pub/plants/release-57/gtf/arabidopsis_thaliana/
        Arabidopsis_thaliana.TAIR10.57.gtf.gz
    http://ftp.ensemblgenomes.org/pub/plants/release-57/fasta/arabidopsis_thaliana/cdna/
        Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

Three transcripts were selected:
  * AT1G01010.1 (NAC001) — fully annotated: 6 exons, start_codon + stop_codon.
  * AT1G03325.1 — protein-coding fragment: no start_codon feature.
  * AT1G24475.1 — protein-coding fragment: no stop_codon feature.
"""
from __future__ import absolute_import

from pyensembl import Genome

from .common import eq_, ok_
from .data import data_path


TAIR10_GTF_PATH = data_path("arabidopsis.tair10.partial.gtf")
TAIR10_CDNA_FASTA_PATH = data_path("arabidopsis.tair10.partial.cdna.fa")


custom_tair10_genome_subset = Genome(
    reference_name="TAIR10",
    annotation_name="_test_arabidopsis_tair10_subset",
    gtf_path_or_url=TAIR10_GTF_PATH,
    transcript_fasta_paths_or_urls=[TAIR10_CDNA_FASTA_PATH],
)


def setup_module(module):
    custom_tair10_genome_subset.clear_cache()
    custom_tair10_genome_subset.index()


def test_complete_transcript_with_start_and_stop():
    # AT1G01010.1 has both start_codon and stop_codon features and a cDNA
    # in the FASTA — .complete must be True and the coding sequence must
    # start with ATG and end with a stop codon.
    transcript = custom_tair10_genome_subset.transcript_by_id("AT1G01010.1")
    ok_(transcript.contains_start_codon)
    ok_(transcript.contains_stop_codon)
    ok_(transcript.sequence is not None)
    ok_(transcript.complete)
    cds = transcript.coding_sequence
    ok_(cds is not None)
    eq_(len(cds) % 3, 0)
    eq_(cds[:3], "ATG")
    ok_(cds[-3:] in ("TAA", "TAG", "TGA"))


def test_transcript_missing_start_codon_is_not_complete():
    # AT1G03325.1 has stop_codon but no start_codon. .complete must be
    # False and .coding_sequence must return None rather than raise.
    transcript = custom_tair10_genome_subset.transcript_by_id("AT1G03325.1")
    eq_(transcript.contains_start_codon, False)
    eq_(transcript.contains_stop_codon, True)
    eq_(transcript.complete, False)
    eq_(transcript.coding_sequence, None)


def test_transcript_missing_stop_codon_is_not_complete():
    # AT1G24475.1 has start_codon but no stop_codon. .complete must be
    # False and .coding_sequence must return None rather than raise.
    transcript = custom_tair10_genome_subset.transcript_by_id("AT1G24475.1")
    eq_(transcript.contains_start_codon, True)
    eq_(transcript.contains_stop_codon, False)
    eq_(transcript.complete, False)
    eq_(transcript.coding_sequence, None)

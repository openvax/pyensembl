from pyensembl import find_nearest_locus
from .common import eq_
from .common import run_multiple_genomes


@run_multiple_genomes()
def test_find_nearest_BRAF_exon(genome):
    braf = genome.genes_by_name("BRAF")[0]
    braf_transcripts = braf.transcripts
    exons = braf_transcripts[0].exons
    for exon in exons:
        # immediately before exon
        result_before = find_nearest_locus(
            start=exon.start - 2, end=exon.start - 1, loci=exons
        )
        eq_(result_before, (1, exon))

        # overlapping with exon
        result_overlap = find_nearest_locus(
            start=exon.start - 2, end=exon.start + 1, loci=exons
        )
        eq_(result_overlap, (0, exon))

        # immediately after exon
        result_after = find_nearest_locus(
            start=exon.end + 1, end=exon.end + 2, loci=exons
        )
        eq_(result_after, (1, exon))


@run_multiple_genomes()
def test_find_nearest_BRAF_transcript(genome):
    braf_transcript = genome.genes_by_name("BRAF")[0].transcripts[0]
    egfr_transcript = genome.genes_by_name("EGFR")[0].transcripts[0]
    transcripts = [braf_transcript, egfr_transcript]
    for transcript in transcripts:
        # immediately before transcript
        result_before = find_nearest_locus(
            start=transcript.start - 2, end=transcript.start - 1, loci=transcripts
        )
        eq_(result_before, (1, transcript))

        # overlapping with transcript
        result_overlap = find_nearest_locus(
            start=transcript.start - 2, end=transcript.start + 1, loci=transcripts
        )
        eq_(result_overlap, (0, transcript))

        # immediately after transcript
        # may overlap with other transcripts
        result_after = find_nearest_locus(
            start=transcript.end + 1, end=transcript.end + 2, loci=transcripts
        )
        eq_(result_after, (1, transcript))

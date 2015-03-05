from test_common import test_ensembl_releases
from pyensembl import find_nearest_locus
from nose.tools import eq_

@test_ensembl_releases
def test_find_nearest_BRAF_exon(ensembl):
    braf = ensembl.genes_by_name("BRAF")[0]
    braf_transcripts = braf.transcripts
    exons = braf_transcripts[0].exons
    for exon in exons:
        # immediately before exon
        result = find_nearest_locus(
            start=exon.start-2,
            end=exon.end-1,
            loci=exons)
        eq_(result, (1, exon))

        # overlapping with exon
        result = find_nearest_locus(
            start=exon.start-2,
            end=exon.start+1,
            loci=exons)
        eq_(result, (0, exon))

        # immediately after exon
        result = find_nearest_locus(
            start=exon.end+1,
            end=exon.end+2,
            loci=exons)
        eq_(result, (1, exon))

@test_ensembl_releases
def test_find_nearest_BRAF_transcript(ensembl):
    braf = ensembl.genes_by_name("BRAF")[0]
    braf_transcripts = braf.transcripts
    for transcript in braf_transcripts:
        # immediately before transcript
        result = find_nearest_locus(
            start=transcript.start-2,
            end=transcript.end-1,
            loci=braf_transcripts)
        eq_(result, (1, transcript))

        # overlapping with transcript
        result = find_nearest_locus(
            start=transcript.start-2,
            end=transcript.start+1,
            loci=braf_transcripts)
        eq_(result, (0, transcript))

        # immediately after transcript
        result = find_nearest_locus(
            start=transcript.end+1,
            end=transcript.end+2,
            loci=braf_transcripts)
        eq_(result, (1, transcript))

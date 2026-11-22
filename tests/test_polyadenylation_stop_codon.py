"""
Mitochondrial reading frames whose stop codon the poly(A) tail completes (#527).

NCBI's rCRS record (NC_012920.1) notes on six human CDS that the "TAA stop
codon is completed by the addition of 3' A residues to the mRNA". Ensembl
gives these transcripts no stop_codon row.
"""

import pytest

from pyensembl import Genome

from .common import run_multiple_genomes

POLY_A_STOPS = {"MT-ND1", "MT-ND2", "MT-CO3", "MT-ND3", "MT-ND4", "MT-CYB"}


@run_multiple_genomes
def test_human_mitochondrial_transcripts(genome):
    transcripts = [t for t in genome.transcripts(contig="MT") if t.biotype == "protein_coding"]
    assert len(transcripts) == 13
    for t in transcripts:
        cds = t.coding_sequence
        if t.gene_name in POLY_A_STOPS:
            assert t.stop_codon_completed_by_polyadenylation
            assert not t.contains_stop_codon and not t.complete
            assert len(cds) % 3 in (1, 2)
            assert (cds + "A" * (-len(cds) % 3))[-3:] == "TAA"
            assert len(cds) // 3 == len(t.protein_sequence)
            assert t.three_prime_utr_sequence == ""
        else:
            assert not t.stop_codon_completed_by_polyadenylation
            assert t.complete
            assert len(cds) // 3 == len(t.protein_sequence) + 1


CODONS = "ATG" + "GCC" * 5


def transcript(tmp_path, tail, contig="MT", strand="+", utr="", start_codon=True,
               transcript_attributes="", fasta_id="T1"):
    """One single-exon transcript: CODONS + tail, then utr, with no stop_codon row."""
    sequence = CODONS + tail + utr
    length, coding = len(sequence), len(CODONS + tail)
    if strand == "+":
        cds, start = (1, coding), (1, 3)
    else:
        cds, start = (length - coding + 1, length), (length - 2, length)
    ids = 'gene_id "G1"; transcript_id "T1"; gene_name "G";'
    rows = [
        ("gene", (1, length), ".", ids + ' gene_biotype "protein_coding";'),
        ("transcript", (1, length), ".",
         ids + ' transcript_biotype "protein_coding"; ' + transcript_attributes),
        ("exon", (1, length), ".", ids + ' exon_number "1"; exon_id "E1";'),
        ("CDS", cds, "0", ids + ' exon_number "1"; protein_id "P1";'),
    ]
    if start_codon:
        rows.append(("start_codon", start, "0", ids + ' protein_id "P1";'))
    (tmp_path / "t.gtf").write_text("".join(
        "%s\ttest\t%s\t%d\t%d\t.\t%s\t%s\t%s\n" % (contig, feature, a, b, strand, frame, attributes)
        for feature, (a, b), frame, attributes in rows))
    (tmp_path / "t.cdna.fa").write_text(">%s\n%s\n" % (fasta_id, sequence))
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="_test_poly_a_stop",
        gtf_path_or_url=str(tmp_path / "t.gtf"),
        transcript_fasta_paths_or_urls=[str(tmp_path / "t.cdna.fa")],
        cache_directory_path=str(tmp_path / "cache"),
    )
    genome.index()
    return genome, genome.transcript_by_id("T1")


@pytest.mark.parametrize("tail", ["TA", "T", "ta"])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_poly_a_completes_t_or_ta(tmp_path, tail, strand):
    genome, t = transcript(tmp_path, tail, strand=strand)
    with genome:
        assert t.stop_codon_completed_by_polyadenylation
        assert t.coding_sequence == CODONS + tail
        assert t.five_prime_utr_sequence == ""
        assert t.three_prime_utr_sequence == ""
        assert not t.complete


@pytest.mark.parametrize("why,kwargs", [
    ("partial codon AG; poly(A) makes AGA", dict(tail="AG")),
    ("ends on a full codon", dict(tail="")),
    ("ends on a full TAA without a stop_codon row", dict(tail="TAA")),
    ("the CDS stops before a 3' UTR", dict(tail="TA", utr="GGG")),
    ("nuclear contig", dict(tail="TA", contig="1")),
    ("no start codon", dict(tail="TA", start_codon=False)),
    ("TAA is tyrosine in NCBI table 14", dict(tail="TA", transcript_attributes='transl_table "14";')),
    ("TAA is tyrosine in NCBI table 33", dict(tail="TA", transcript_attributes='transl_table "33";')),
])
def test_other_missing_stops_stay_unknown(tmp_path, why, kwargs):
    genome, t = transcript(tmp_path, **kwargs)
    with genome:
        assert not t.stop_codon_completed_by_polyadenylation, why
        assert t.coding_sequence is None
        assert t.three_prime_utr_sequence is None
        assert not t.complete


def test_vertebrate_mitochondrial_table_is_explicit(tmp_path):
    genome, t = transcript(tmp_path, "TA", transcript_attributes='transl_table "2";')
    with genome:
        assert t.stop_codon_completed_by_polyadenylation


def test_unavailable_sequence_is_not_flagged(tmp_path):
    genome, t = transcript(tmp_path, "TA", fasta_id="OTHER")
    with genome:
        assert t.sequence is None
        assert not t.stop_codon_completed_by_polyadenylation
        assert t.coding_sequence is None

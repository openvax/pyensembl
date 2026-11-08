"""Lookups in annotations that omit a feature type or an attribute.

The database has a table only for the features a GTF contains and a column
only for the attributes it uses, so a lookup on either finds nothing instead
of raising an SQLite error (#499, #502).
"""

import pytest

from pyensembl import Genome

GENE = 'gene_id "g"; gene_name "G"; gene_biotype "protein_coding"; '
TRANSCRIPT = 'transcript_id "t"; transcript_biotype "protein_coding"; '
EXON = 'exon_id "e"; exon_number "1";'


def indexed_genome(directory, attributes):
    """Index a one-gene GTF with gene, transcript and exon rows only."""
    gtf = directory / "partial.gtf"
    gtf.write_text("".join(
        "1\ttest\t%s\t100\t200\t.\t+\t.\t%s\n" % (feature, attributes)
        for feature in ("gene", "transcript", "exon")
    ))
    genome = Genome(reference_name="test", annotation_name="partial",
                    gtf_path_or_url=str(gtf), cache_directory_path=str(directory))
    genome.index()
    return genome


def test_transcript_without_codon_or_cds_tables_is_incomplete(tmp_path):
    with indexed_genome(tmp_path, GENE + TRANSCRIPT + EXON) as genome:
        assert not any(genome.db.columns(feature)
                       for feature in ("start_codon", "stop_codon", "CDS"))
        transcript = genome.transcript_by_id("t")
        assert not transcript.contains_start_codon
        assert not transcript.contains_stop_codon
        assert not transcript.start_codon_complete
        assert not transcript.stop_codon_complete
        assert not transcript.complete


@pytest.mark.parametrize("accessor, feature", [
    ("start_codon_positions", "start_codon"),
    ("stop_codon_positions", "stop_codon"),
    ("coding_sequence_position_ranges", "CDS"),
])
def test_required_transcript_feature_still_raises(tmp_path, accessor, feature):
    with indexed_genome(tmp_path, GENE + TRANSCRIPT + EXON) as genome:
        transcript = genome.transcript_by_id("t")
        with pytest.raises(ValueError, match="does not contain feature %s" % feature):
            getattr(transcript, accessor)


def test_transcript_name_lookup_without_transcript_names(tmp_path):
    with indexed_genome(tmp_path, GENE + TRANSCRIPT + EXON) as genome:
        assert not genome.db.column_exists("transcript", "transcript_name")
        # The same error as for a name absent from an annotation with names;
        # the gene name "G" is not used as a transcript name.
        for name in ("CFTR", "G"):
            with pytest.raises(ValueError, match="No results found"):
                genome.transcripts_by_name(name)
            with pytest.raises(ValueError, match="No results found"):
                genome.exon_ids_of_transcript_name(name)
        assert genome.transcript_ids_of_gene_name("G") == ["t"]
        assert genome.transcript_by_id("t").transcript_name is None


def test_transcript_name_lookup_with_transcript_names(tmp_path):
    attributes = GENE + TRANSCRIPT + 'transcript_name "G-201"; ' + EXON
    with indexed_genome(tmp_path, attributes) as genome:
        assert [t.id for t in genome.transcripts_by_name("G-201")] == ["t"]
        with pytest.raises(ValueError, match="No results found"):
            genome.transcripts_by_name("CFTR")

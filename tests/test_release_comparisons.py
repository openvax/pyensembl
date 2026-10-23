"""Independent source excerpts and real grouped/coordinate-sorted indexing."""

import gzip
from pathlib import Path

import pytest

from pyensembl import Genome
from scripts import capture_release_comparisons as capture


@pytest.mark.parametrize("stem", ["human-116", "human-2026_04", "dog-116", "dog-2020_10"])
def test_capture_parser_is_order_independent_but_keeps_coordinates_and_attributes(tmp_path, stem):
    source = Path(__file__).parent / "data/release_comparisons" / (stem + ".gtf")
    rows = [line for line in source.read_text().splitlines(True) if not line.startswith("#")]

    def summarize(name, lines):
        path = tmp_path / (name + ".gtf.gz")
        path.write_bytes(gzip.compress("".join(lines).encode()))
        return capture.summarize(dict(path=str(path), species=name, release=1))

    original = summarize("original", rows)
    reordered = summarize("reordered", list(reversed(rows)))
    assert original["normalized_gtf_sha256"] == reordered["normalized_gtf_sha256"]
    assert original["gene_ids_sha256"] == reordered["gene_ids_sha256"]
    assert original["genes"] > 0
    assert sum(original["feature_counts"].values()) == len(rows)
    changed = rows.copy()
    fields = changed[0].split("\t")
    fields[3] = str(int(fields[3]) + 1)
    changed[0] = "\t".join(fields)
    altered = summarize("altered", changed)
    assert original["normalized_gtf_sha256"] != altered["normalized_gtf_sha256"]
    assert original["gene_ids_sha256"] == altered["gene_ids_sha256"]
    duplicated = summarize("duplicated", rows + [rows[0]])
    assert original["normalized_gtf_sha256"] != duplicated["normalized_gtf_sha256"]


@pytest.mark.parametrize("ordering", ["coordinates", "children-first"])
def test_interleaved_gtf_retains_gene_transcript_exon_and_coding_queries(tmp_path, ordering):
    fixture = Path(__file__).parent / "data/mouse.ensembl.81.partial.ENSMUSG00000017167.gtf"
    rows = fixture.read_text().splitlines(True)

    def key(line):
        fields = line.split("\t")
        # Equal-coordinate child rows precede their parent; transcripts from
        # overlapping isoforms interleave instead of forming contiguous blocks.
        rank = {"exon": 0, "CDS": 1, "transcript": 2, "gene": 3}.get(fields[2], 0)
        return (fields[0], int(fields[3]), rank, int(fields[4]))

    ordered = sorted(rows, key=key)
    if ordering == "children-first":
        ordered.sort(key=lambda line: line.split("\t")[2] in ("transcript", "gene"))
    assert ordered != rows and ordered[0].split("\t")[2] == "exon"

    def query(name, lines):
        path = tmp_path / (name + ".gtf")
        path.write_text("".join(lines))
        with Genome(reference_name="GRCm38", annotation_name=name,
                    gtf_path_or_url=str(path), cache_directory_path=str(tmp_path / name)) as genome:
            genome.index()
            transcripts = genome.transcripts()
            assert len(transcripts) >= 2
            result = {}
            for transcript in transcripts:
                result[transcript.id] = (
                    transcript.gene_id, transcript.start, transcript.end,
                    [(e.id, e.start, e.end) for e in transcript.exons],
                    transcript.coding_sequence_position_ranges if transcript.biotype == "protein_coding" else [],
                    transcript.start_codon_positions if transcript.biotype == "protein_coding" else [],
                    transcript.stop_codon_positions if transcript.biotype == "protein_coding" else [],
                )
            assert genome.gene_by_id("ENSMUSG00000017167").name == "Cntnap1"
            assert result["ENSMUST00000103109"][4]
            return result

    assert query("grouped", rows) == query(ordering, ordered)

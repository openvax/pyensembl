"""Exercise pandas parsing through the complete annotation indexing path."""

import gtfparse
import pandas as pd
import pytest

from pyensembl import Genome


ROWS = (
    'chr1\ttest\tgene\t100\t400\t.\t-\t.\tgene_id "g1"; gene_name "GENE1"; gene_type "protein_coding"; gene_version "7";\n',
    'chr1\ttest\ttranscript\t100\t400\t.\t-\t.\tgene_id "g1"; gene_name "GENE1"; transcript_id "t1"; transcript_type "protein_coding"; transcript_version "3";\n',
    'chr1\ttest\texon\t100\t200\t.\t-\t.\tgene_id "g1"; gene_name "GENE1"; transcript_id "t1"; exon_id "e1"; exon_number "2"; gene_type "protein_coding"; transcript_type "protein_coding";\n',
    'chr1\ttest\texon\t300\t400\t.\t-\t.\tgene_id "g1"; gene_name "GENE1"; transcript_id "t1"; exon_id "e2"; exon_number "1"; gene_type "protein_coding"; transcript_type "protein_coding";\n',
)


@pytest.mark.parametrize("missing_features", [False, True])
@pytest.mark.parametrize("show_progress", [False, True])
def test_pandas_indexing_preserves_annotation(tmp_path, monkeypatch, missing_features, show_progress):
    source = tmp_path / "annotation.gtf"
    source.write_text("".join(ROWS[2:] if missing_features else ROWS))
    parsed = []
    read_gtf = gtfparse.read_gtf

    def record_parse(*args, **kwargs):
        frame = read_gtf(*args, **kwargs)
        parsed.append((frame, kwargs["progress_callback"]))
        return frame

    monkeypatch.setattr(gtfparse, "read_gtf", record_parse)
    with Genome("synthetic", "pandas", gtf_path_or_url=str(source),
                cache_directory_path=str(tmp_path / "cache")) as genome:
        genome.index(show_progress=show_progress)
        assert isinstance(parsed[0][0], pd.DataFrame)
        assert (parsed[0][1] is not None) == show_progress
        gene = genome.gene_by_id("g1")
        transcript = genome.transcript_by_id("t1")
        assert (gene.name, gene.contig, gene.start, gene.end, gene.strand) == (
            "GENE1", "chr1", 100, 400, "-"
        )
        assert gene.biotype == transcript.biotype == "protein_coding"
        assert transcript.gene_id == "g1"
        assert [exon.id for exon in transcript.exons] == ["e2", "e1"]
        if not missing_features:
            assert gene.version == 7
            assert transcript.version == 3
            assert parsed[0][0].loc[0, "gene_version"] == "7"
            assert parsed[0][0].loc[1, "transcript_version"] == "3"

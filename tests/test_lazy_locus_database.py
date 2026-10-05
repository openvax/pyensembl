"""Restoring complete locus metadata must not require annotation resources."""

import pytest
from serializable import from_json, to_json

from pyensembl import Gene, Genome, Transcript
from pyensembl.download_cache import MissingLocalFile
from .data import data_path


@pytest.mark.parametrize("cls,identity", [
    (Gene, dict(gene_id="gene", gene_name="GENE")),
    (Transcript, dict(transcript_id="tx", transcript_name="TX", gene_id="gene")),
])
def test_unavailable_reference_allows_identity_roundtrip_but_not_annotation(tmp_path, cls, identity):
    missing_gtf = str(tmp_path / "original-location" / "reference.gtf")
    genome = Genome(reference_name="custom", annotation_name="pinned", annotation_version=1,
                    gtf_path_or_url=missing_gtf, cache_directory_path=str(tmp_path / "cache"))
    locus = cls(contig="1", start=10, end=20, strand="+", biotype="protein_coding",
                genome=genome, **identity)
    restored = from_json(to_json(locus))
    assert restored == locus
    assert (restored.contig, restored.start, restored.end, restored.strand) == ("1", 10, 20, "+")
    assert restored.gene_id == "gene"
    assert restored.genome.to_dict() == genome.to_dict()
    assert restored.is_protein_coding
    with pytest.raises(MissingLocalFile, match="reference.gtf"):
        restored.db
    with pytest.raises(MissingLocalFile, match="reference.gtf"):
        restored.exons if isinstance(restored, Transcript) else restored.transcripts


def test_available_reference_resolves_real_database_after_metadata_roundtrip(tmp_path):
    genome = Genome(reference_name="GRCm38", annotation_name="test", annotation_version=81,
        gtf_path_or_url=data_path("mouse.ensembl.81.partial.ENSMUSG00000017167.gtf"),
        cache_directory_path=str(tmp_path))
    genome.index()
    transcript = genome.transcripts()[0]
    restored = from_json(to_json(transcript))
    assert restored.db is restored.genome.db
    assert [(e.contig, e.start, e.end) for e in restored.exons] == [
        (e.contig, e.start, e.end) for e in transcript.exons]
    assert restored.gene.id == transcript.gene.id
    assert from_json(to_json(transcript.gene)).transcripts == transcript.gene.transcripts

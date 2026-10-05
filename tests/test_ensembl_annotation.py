"""Install assembly/date datasets through the actual Genome pipeline."""

import gzip
from pathlib import Path
import struct
import zlib

import pytest

from pyensembl import EnsemblAnnotation


def bgzf(data):
    """One small BGZF member followed by its empty end-of-file member."""
    def block(payload):
        compressor = zlib.compressobj(wbits=-15)
        body = compressor.compress(payload) + compressor.flush()
        size = 18 + len(body) + 8
        header = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff\x06\x00BC\x02\x00"
        return header + struct.pack("<H", size - 1) + body + struct.pack(
            "<II", zlib.crc32(payload), len(payload)
        )
    return block(data) + block(b"")


@pytest.mark.parametrize("include_alt", [False, True])
def test_install_index_query_and_reopen_dated_annotation(tmp_path, include_alt):
    mirror = tmp_path / "mirror"
    root = mirror / "GCA/000/001/405/29/ensembl/2023_03"
    geneset = root / "geneset"
    geneset.mkdir(parents=True)
    genome_dir = root / "genome"
    genome_dir.mkdir()
    attrs = 'gene_id "ENSG1"; gene_name "TEST"; transcript_id "ENST1"; transcript_version "2";'
    rows = (
        '1\tensembl\tgene\t3\t12\t.\t-\t.\tgene_id "ENSG1"; gene_name "TEST";\n'
        f'1\tensembl\ttranscript\t3\t12\t.\t-\t.\t{attrs}\n'
        f'1\tensembl\texon\t3\t6\t.\t-\t.\t{attrs} exon_id "e2"; exon_number "2";\n'
        f'1\tensembl\texon\t9\t12\t.\t-\t.\t{attrs} exon_id "e1"; exon_number "1";\n'
    )
    (geneset / "genes.gtf.gz").write_bytes(gzip.compress(rows.encode()))
    alt = 'ALT\tensembl\tgene\t1\t4\t.\t+\t.\tgene_id "alt"; gene_name "TEST";\n'
    (geneset / "genes-including_alt.gtf.gz").write_bytes(gzip.compress((rows + alt).encode()))
    (geneset / "cdna.fa.bgz").write_bytes(bgzf(b">ENST1.2\nCCCCGGTT\n"))
    (geneset / "pep.fa.bgz").write_bytes(bgzf(b">ENSP1.1\nMP\n"))
    (genome_dir / "softmasked.fa.bgz").write_bytes(bgzf(b">1\nAAAAccccGGGG\n>ALT\nACGT\n"))
    options = dict(include_alt=include_alt, genome_fasta=True, genome_fasta_mask="soft",
                   species="human", server=mirror.as_uri(),
                   cache_directory_path=str(tmp_path / "cache"))
    with EnsemblAnnotation("GCA_000001405.29", "2023_03", **options) as annotation:
        # Construction is offline and does not create the cache.
        assert not Path(options["cache_directory_path"]).exists()
        annotation.download()
        annotation.index()
        assert annotation.gene_by_id("ENSG1").strand == "-"
        assert annotation.gene_ids_of_gene_name("TEST") == (["ENSG1", "alt"] if include_alt else ["ENSG1"])
        transcript = annotation.transcript_by_id("ENST1")
        assert transcript.sequence == "CCCCGGTT"
        assert [exon.id for exon in transcript.exons] == ["e1", "e2"]
        assert annotation.protein_sequences.get("ENSP1.1") == "MP"
        assert annotation.sequence("1", 3, 6, strand="-", mask="raw") == "ggTT"
        encoded = annotation.to_json()
    # Local indexes and optional DNA survive a JSON round trip.
    with EnsemblAnnotation.from_json(encoded) as reopened:
        assert reopened.species.latin_name == "homo_sapiens"
        assert reopened.transcript_by_id("ENST1").sequence == "CCCCGGTT"
        assert reopened.sequence("1", 9, 12, strand="-") == "CCCC"


def test_cache_selections_are_distinct_and_dna_is_optional(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    options = [dict(), dict(include_alt=True), dict(provider="community"),
               dict(annotation_date="2023_04"), dict(assembly_accession="GCF_000001405.29")]
    directories = set()
    for changes in options:
        kwargs = dict(assembly_accession="GCA_000001405.29", annotation_date="2023_03",
                      reference_name="GRCh38")
        with EnsemblAnnotation(**{**kwargs, **changes}) as annotation:
            assert not annotation.requires_genome_fasta
            directories.add(annotation.download_cache.cache_directory_path)
    assert len(directories) == len(options)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("options", [
    {"assembly_accession": "GRCh38"}, {"assembly_accession": "GCA_000001405.0"},
    {"annotation_date": "2026-07"}, {"annotation_date": "2023_13"},
    {"provider": "../community"}, {"genome_fasta_mask": "invalid"},
])
def test_invalid_dataset_selection_rejected(options):
    kwargs = dict(assembly_accession="GCA_000001405.29", annotation_date="2023_03")
    with pytest.raises(ValueError):
        EnsemblAnnotation(**{**kwargs, **options})


def test_non_boolean_coverage_rejected():
    with pytest.raises(TypeError):
        EnsemblAnnotation("GCA_000001405.29", "2023_03", include_alt="False")

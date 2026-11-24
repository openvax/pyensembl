"""Tests for native T2T-CHM13v2.0 and pinned rapid-release assemblies (#541)."""

import gzip
import hashlib
import json
from dataclasses import replace

import pytest

import pyensembl
from pyensembl import (
    ContigAlias,
    PinnedAssembly,
    PinnedGenome,
    RCRS_ACCESSION,
    ReferenceAlias,
    available_releases,
    contig_aliases,
    find_pinned_reference_name,
    find_species_by_reference,
    genome_for_reference_name,
    is_rcrs_reference,
    mitochondrion_note,
    normalize_pinned_contig,
    normalize_reference_name,
    pinned_assemblies,
    pinned_contig_lengths,
    pinned_contig_spellings,
    pinned_reference_aliases,
    resolve_pinned_assembly,
    shell,
)
from pyensembl.annotation_integrity import AnnotationIntegrityError
from pyensembl.pinned_genomes import _all_pinned_assemblies


def _make_mini_rapid_assets(tmp_path):
    """Build tiny valid gzipped GTF/cDNA/pep files for hermetic PinnedGenome tests."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    gtf_text = (
        '1\tensembl\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "TEST1"; gene_biotype "protein_coding";\n'
        '1\tensembl\ttranscript\t100\t200\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; transcript_biotype "protein_coding";\n'
        '1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1"; exon_id "ENSE00000000001";\n'
        '1\tensembl\tCDS\t100\t108\t.\t+\t0\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; protein_id "ENSP00000000001";\n'
    )
    cdna_text = ">ENST00000000001\nATGGCTTAA\n"
    pep_text = ">ENSP00000000001\nMA\n"

    paths = {}
    specs = {}
    for role, filename, content in (
        ("gtf", "Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz", gtf_text),
        ("transcripts", "Homo_sapiens-GCA_009914755.4-2022_07-cdna.fa.gz", cdna_text),
        ("proteins", "Homo_sapiens-GCA_009914755.4-2022_07-pep.fa.gz", pep_text),
    ):
        path = tmp_path / filename
        raw = gzip.compress(content.encode("utf-8"), mtime=0)
        path.write_bytes(raw)
        paths[role] = path
        specs[role] = {
            "url": "https://example.invalid/%s" % filename,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "size": len(raw),
        }
    return paths, specs


def _patch_pinned_assets(monkeypatch, specs):
    real = resolve_pinned_assembly("t2t-chm13")
    patched = replace(real, assets=specs)
    monkeypatch.setattr(
        pyensembl.pinned_genomes,
        "_all_pinned_assemblies",
        lambda: (patched,),
    )
    return patched


def test_catalog_metadata_and_aliases():
    assemblies = pinned_assemblies()
    assert len(assemblies) >= 1
    assert any(a.name == "t2t-chm13" and a.reference_name == "T2T-CHM13v2.0" for a in assemblies)
    assert len(pinned_assemblies(species="human")) >= 1
    assert pinned_assemblies(species="mouse") == ()

    for alias in (
        "t2t-chm13",
        "T2T-CHM13v2.0",
        "t2t_chm13v2.0",
        "hs1",
        "HS1",
        "GCA_009914755.4",
        "GCF_009914755.1",
        "chm13v2.0",
        "chm13v2.0_maskedY",
        "chm13v2.0_noY",
        "chm13v2.0_maskedY_rCRS",
        "chm13v2.0_maskedY.rCRS",
        "T2T-CHM13v1.1",
        "GCA_009914755.3",
        "chm13.draft_v1.1",
    ):
        resolved = resolve_pinned_assembly(alias)
        assert isinstance(resolved, PinnedAssembly)
        assert resolved.reference_name == "T2T-CHM13v2.0"
        assert resolved.assembly_accession == "GCA_009914755.4"
        assert resolved.annotation_name == "ensembl-rapid"
        assert resolved.annotation_version == "2022_07"

    aliases = pinned_reference_aliases("T2T-CHM13v2.0")
    assert all(isinstance(row, ReferenceAlias) for row in aliases)
    assert any(row.name == "hs1" and row.assembly == "t2t-chm13" for row in aliases)
    assert any(row.name == "T2T-CHM13v1.0" and row.assembly is None for row in aliases)


@pytest.mark.parametrize(
    "unsupported",
    [
        "T2T-CHM13v1.0",
        "T2T-CHM13 v1.0",
        "t2t-chm13.v1.0",
        "GCA_009914755.2",
        "chm13.draft_v1.0",
        "chm13.draft_v0.9",
        "T2T-CHM13v0.7",
        "GCA_009914755.1",
        "chm13.draft_v0.7",
    ],
)
def test_other_versions_rejected_with_explanation(unsupported):
    with pytest.raises(ValueError, match="coordinates differ|draft"):
        resolve_pinned_assembly(unsupported)
    with pytest.raises(ValueError, match="coordinates differ|draft"):
        normalize_reference_name(unsupported)
    with pytest.raises(ValueError, match="coordinates differ|draft"):
        PinnedGenome(unsupported)


def test_find_pinned_reference_name_and_rcrs():
    assert find_pinned_reference_name("/refs/chm13v2.0_maskedY_rCRS.fa") == "chm13v2.0_maskedY_rCRS"
    assert find_pinned_reference_name("/refs/hs1/sample.vcf.gz") == "hs1"
    assert find_pinned_reference_name("T2T-CHM13 v1.0.fa") == "T2T-CHM13v1.0"
    assert find_pinned_reference_name("hs10.fa") is None
    assert find_pinned_reference_name("ths1.fa") is None

    assert is_rcrs_reference("chm13v2.0_maskedY_rCRS") is True
    assert is_rcrs_reference("chm13v2.0_maskedY.rCRS") is True
    assert is_rcrs_reference("hs1") is False
    assert is_rcrs_reference("T2T-CHM13v2.0") is False
    assert is_rcrs_reference("GCF_009914755.1") is False


def test_contig_aliases_and_mitochondrion_disambiguation():
    native = contig_aliases("T2T-CHM13v2.0")
    assert all(isinstance(c, ContigAlias) for c in native)
    assert len(native) == 25
    mt_row = next(c for c in native if c.name == "MT")
    assert mt_row.ucsc == "chrM"
    assert mt_row.genbank == "CP068254.1"
    assert mt_row.length == 16569

    rcrs = contig_aliases("chm13v2.0_maskedY_rCRS")
    assert all(c.name != "MT" for c in rcrs)
    rcrs_row = next(c for c in rcrs if c.name == RCRS_ACCESSION)
    assert rcrs_row.ucsc == "chrM"
    assert rcrs_row.genbank == "J01415.2"
    assert rcrs_row.refseq == RCRS_ACCESSION
    assert rcrs_row.length == 16569

    noy = contig_aliases("chm13v2.0_noY")
    assert all(c.name != "Y" for c in noy)
    assert len(noy) == 24

    refseq = contig_aliases("GCF_009914755.1")
    assert all(c.name != "MT" for c in refseq)
    assert len(refseq) == 24

    assert normalize_pinned_contig("chr22", "hs1") == "22"
    assert normalize_pinned_contig("CP068256.2", "T2T-CHM13v2.0") == "22"
    assert normalize_pinned_contig("NC_060946.1", "T2T-CHM13v2.0") == "22"
    assert normalize_pinned_contig("chrM", "T2T-CHM13v2.0", chrm_is_rcrs=False) == "MT"
    assert normalize_pinned_contig("chrM", "T2T-CHM13v2.0", chrm_is_rcrs=True) == RCRS_ACCESSION
    assert normalize_pinned_contig("M", "T2T-CHM13v2.0", chrm_is_rcrs=True) == RCRS_ACCESSION
    assert normalize_pinned_contig("chr1", "GRCh38") is None

    spellings = pinned_contig_spellings("hs1", chrm_is_rcrs=True)
    assert RCRS_ACCESSION in spellings
    assert {"chrM", "M", "J01415.2", RCRS_ACCESSION}.issubset(spellings[RCRS_ACCESSION])
    without_mt = pinned_contig_spellings("hs1", chrm_is_rcrs=False, include_mitochondrion=False)
    assert "MT" not in without_mt

    lengths = pinned_contig_lengths("hs1")
    assert lengths["1"] == 248387328
    assert lengths["MT"] == 16569

    assert "577" in mitochondrion_note("MT", "T2T-CHM13v2.0")
    assert RCRS_ACCESSION in mitochondrion_note(RCRS_ACCESSION, "hs1")
    assert mitochondrion_note("1", "T2T-CHM13v2.0") is None


def test_reference_name_and_available_releases_integration():
    assert normalize_reference_name("hs1") == "T2T-CHM13v2.0"
    assert normalize_reference_name("GCA_009914755.4") == "T2T-CHM13v2.0"
    assert find_species_by_reference("hs1").latin_name == "homo_sapiens"

    genome = genome_for_reference_name("hs1")
    assert isinstance(genome, PinnedGenome)
    assert genome.reference_name == "T2T-CHM13v2.0"
    assert genome.assembly_accession == "GCA_009914755.4"
    assert genome.species.latin_name == "homo_sapiens"
    assert genome_for_reference_name("T2T-CHM13v2.0") is genome

    for kwargs in (
        dict(assembly="hs1"),
        dict(assembly="T2T-CHM13v2.0"),
        dict(assembly_accession="GCA_009914755.4"),
        dict(species="human", assembly="t2t-chm13"),
    ):
        rows = available_releases(**kwargs)
        assert len(rows) == 1
        row = rows[0]
        assert row["kind"] == "pinned"
        assert row["release"] == "2022_07"
        assert row["assembly_accession"] == "GCA_009914755.4"
        assert row["provider"] == "ensembl-rapid"
        assert row["latest"] is True

    # Default human selection must remain GRCh38.
    human_rows = available_releases("human")
    assert all(r["assembly"] == "GRCh38" for r in human_rows)


def test_pinned_genome_hermetic_install_index_verify_and_roundtrip(tmp_path, monkeypatch):
    source_paths, specs = _make_mini_rapid_assets(tmp_path / "src")
    _patch_pinned_assets(monkeypatch, specs)

    cache_dir = tmp_path / "cache" / "T2T-CHM13v2.0" / "ensembl-rapid2022_07"
    genome = PinnedGenome("hs1", cache_directory_path=str(cache_dir))
    assert genome.installed() is False
    assert genome.inspect_data()["annotation"] == "missing"
    assert genome.genetic_code("1").table_id == 1
    assert genome.genetic_code(RCRS_ACCESSION).table_id == 2

    # Tampered source path fails SHA-256 verification before caching.
    bad_gtf = tmp_path / "bad.gtf.gz"
    bad_gtf.write_bytes(b"corrupted")
    with pytest.raises(AnnotationIntegrityError, match="Size mismatch|SHA-256 mismatch"):
        genome.download(source_paths={**source_paths, "gtf": bad_gtf})
    assert genome.installed() is False

    # Valid pinned install + index.
    genome.download(source_paths=source_paths)
    genome.index()
    assert genome.installed() is True

    # Nuclear rapid-release gene AND bundled rCRS mitochondrial genes are both indexed.
    assert genome.gene_by_id("ENSG00000000001").gene_name == "TEST1"
    assert genome.transcript_by_id("ENST00000000001").sequence == "ATGGCTTAA"
    assert genome.protein_sequence("ENSP00000000001") == "MA"
    rcrs_genes = genome.genes(contig=RCRS_ACCESSION)
    assert len(rcrs_genes) > 0
    assert any(g.gene_name == "MT-CO1" for g in rcrs_genes)
    mt_tx = next(t for t in genome.transcripts(contig=RCRS_ACCESSION) if t.is_protein_coding)
    assert mt_tx.genetic_code.table_id == 2

    # Offline verification passes for all assets + rCRS layers + merged GTF.
    reports = genome.verify_sources(check_remote=True)
    assert len(reports) >= 6
    assert all(r["status"] == "verified" for r in reports)
    assert all(r["upstream_status"] == "pinned" for r in reports)

    # Portable serialization round-trip.
    reloaded = PinnedGenome.from_dict(genome.to_dict())
    assert isinstance(reloaded, PinnedGenome)
    assert reloaded.reference_name == "T2T-CHM13v2.0"
    assert reloaded.assembly_accession == "GCA_009914755.4"
    assert reloaded.installed() is True
    assert reloaded.gene_by_id("ENSG00000000001").gene_name == "TEST1"

    # Mutating a cached file is caught by verify_sources().
    raw_pep = genome.download_cache.local_path(specs["proteins"]["url"])
    with open(raw_pep, "ab") as handle:
        handle.write(b"x")
    after = {r["role"]: r for r in genome.verify_sources()}
    assert after["proteins"]["status"] == "changed"


def test_cli_pinned_selection_available_list_inspect_verify(tmp_path, monkeypatch):
    source_paths, specs = _make_mini_rapid_assets(tmp_path / "src")
    _patch_pinned_assets(monkeypatch, specs)
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))

    # available table and --assembly JSON
    table = shell.format_available_species(use_color=False)
    assert "T2T-CHM13v2.0" in table
    assert "2022_07 (pinned)" in table

    avail_args = shell.parser.parse_args(["available", "--assembly", "hs1", "--json"])
    payload = json.loads(
        shell.run_available(
            species=avail_args.species,
            assembly=avail_args.assembly,
            assembly_accession=avail_args.assembly_accession,
            provider=avail_args.provider,
            reference_name=avail_args.reference_name,
            include_alt=avail_args.include_alt,
            refresh=avail_args.refresh,
            json_output=avail_args.json,
        )
    )
    assert len(payload) == 1
    assert payload[0]["assembly"] == "T2T-CHM13v2.0"
    assert payload[0]["kind"] == "pinned"

    # CLI genome selection via --reference-name / --assembly / --assembly-accession
    for argv in (
        ["install", "--reference-name", "T2T-CHM13v2.0"],
        ["install", "--reference-name", "hs1"],
        ["install", "--assembly", "hs1"],
        ["install", "--assembly-accession", "GCA_009914755.4"],
        ["install", "--reference-name", "T2T-CHM13v2.0", "--release", "2022_07"],
    ):
        genomes = shell.collect_selected_genomes(shell.parser.parse_args(argv))
        assert len(genomes) == 1
        assert isinstance(genomes[0], PinnedGenome)

    # Unsupported version and mismatched release raise ValueError
    with pytest.raises(ValueError, match="coordinates differ"):
        shell.collect_selected_genomes(shell.parser.parse_args(["install", "--reference-name", "T2T-CHM13v1.0"]))
    with pytest.raises(ValueError, match="not available"):
        shell.collect_selected_genomes(
            shell.parser.parse_args(["install", "--reference-name", "T2T-CHM13v2.0", "--release", "110"])
        )
    with pytest.raises(ValueError, match="belongs to homo_sapiens"):
        shell.collect_selected_genomes(
            shell.parser.parse_args(["install", "--reference-name", "hs1", "--species", "mouse"])
        )

    # Install into default PYENSEMBL_CACHE_DIR, then check list, inspect, verify
    genome = PinnedGenome("T2T-CHM13v2.0")
    genome.download(source_paths=source_paths)
    genome.index()

    listed = shell.format_installed_genomes(use_color=False)
    assert "T2T-CHM13v2.0" in listed
    assert "2022_07" in listed
    assert "indexed" in listed

    inspected = json.loads(shell.inspect_genomes([genome], json_output=True))
    assert inspected[0]["installed"] is True
    assert inspected[0]["reference_name"] == "T2T-CHM13v2.0"

    verify_text, ok = shell.run_verify([genome], check_remote=False, json_output=True)
    assert ok is True
    verify_payload = json.loads(verify_text)
    assert verify_payload[0]["release"] == "2022_07"
    assert all(s["status"] == "verified" for s in verify_payload[0]["sources"])

    # Keep _all_pinned_assemblies clean reference check
    assert _all_pinned_assemblies is not None

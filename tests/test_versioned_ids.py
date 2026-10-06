"""
Issue #445: a version suffix on an Ensembl ID is optional, but when given it
must match this annotation. Bare IDs match whatever version is installed.
"""

from os.path import join
from tempfile import TemporaryDirectory

import pytest

from pyensembl import Genome
from pyensembl.versioned_ids import match_version

from .common import eq_, grch38
from .data import custom_genome_cache_path, data_path
from .test_versioned_protein_fasta import _make_genome as _make_gencode_genome
from .test_versions import mouse_genome


def _ensembl_gtf(with_versions=True):
    """One coding transcript: gene v3, transcript v7, exon v2, protein v5."""

    def versions(**kwargs):
        if not with_versions:
            return ""
        return "".join(' %s_version "%d";' % item for item in kwargs.items())

    gene = 'gene_id "ENSGTEST00000020001";%s gene_name "VID1"; gene_biotype "protein_coding";' % versions(gene=3)
    transcript = '%s transcript_id "ENSTTEST00000020001";%s transcript_name "VID1-201"; transcript_biotype "protein_coding";' % (gene, versions(transcript=7))
    exon = '%s exon_number "1"; exon_id "ENSETEST00000020001";%s' % (transcript, versions(exon=2))
    cds = '%s protein_id "ENSPTEST00000020001";%s' % (exon, versions(protein=5))
    rows = [
        ("gene", 100, 800, ".", gene),
        ("transcript", 100, 800, ".", transcript),
        ("exon", 100, 800, ".", exon),
        ("CDS", 100, 798, "0", cds),
    ]
    return "".join(
        "1\ttest\t%s\t%d\t%d\t.\t+\t%s\t%s\n" % (feature, start, end, frame, attributes)
        for feature, start, end, frame, attributes in rows
    )


def _ensembl_genome(tmpdir, gtf_versions=True, fasta_versions=True):
    """Ensembl 83+ records versions in the GTF and FASTA headers, 77-82 only
    in the GTF, and earlier releases nowhere."""
    gtf = join(tmpdir, "e.gtf")
    with open(gtf, "w") as f:
        f.write(_ensembl_gtf(with_versions=gtf_versions))
    cdna = join(tmpdir, "e.cdna.fa")
    pep = join(tmpdir, "e.pep.fa")
    for path, header, sequence in [
        (cdna, "ENSTTEST00000020001.7", "ATGCCCAAATTT"),
        (pep, "ENSPTEST00000020001.5", "MPKF"),
    ]:
        if not fasta_versions:
            header = header.rpartition(".")[0]
        with open(path, "w") as f:
            f.write(">%s\n%s\n" % (header, sequence))
    genome = Genome(
        reference_name="GRCh38",
        annotation_name="_test_versioned_ids",
        gtf_path_or_url=gtf,
        transcript_fasta_paths_or_urls=[cdna],
        protein_fasta_paths_or_urls=[pep],
        cache_directory_path=tmpdir,
    )
    genome.index()
    return genome


@pytest.fixture(params=[True, False], ids=["fasta-versions", "gtf-only-versions"])
def ensembl(request):
    with TemporaryDirectory() as tmpdir:
        yield _ensembl_genome(tmpdir, fasta_versions=request.param)


def test_versioned_ids_return_the_same_objects(ensembl):
    assert ensembl.gene_by_id("ENSGTEST00000020001.3") is ensembl.gene_by_id("ENSGTEST00000020001")
    transcript = ensembl.transcript_by_id("ENSTTEST00000020001.7")
    assert transcript is ensembl.transcript_by_id("ENSTTEST00000020001")
    eq_(transcript.id, "ENSTTEST00000020001")
    eq_(transcript.versioned_id, "ENSTTEST00000020001.7")
    exon = ensembl.exon_by_id("ENSETEST00000020001.2")
    assert exon is ensembl.exon_by_id("ENSETEST00000020001")
    eq_(exon.versioned_id, "ENSETEST00000020001.2")


def test_versioned_ids_in_relationship_lookups(ensembl):
    eq_(ensembl.transcript_by_protein_id("ENSPTEST00000020001.5").id, "ENSTTEST00000020001")
    eq_(ensembl.gene_by_protein_id("ENSPTEST00000020001.5").id, "ENSGTEST00000020001")
    eq_(ensembl.transcript_id_of_protein_id("ENSPTEST00000020001.5"), "ENSTTEST00000020001")
    eq_(ensembl.gene_name_of_gene_id("ENSGTEST00000020001.3"), "VID1")
    eq_(ensembl.transcript_ids_of_gene_id("ENSGTEST00000020001.3"), ["ENSTTEST00000020001"])
    eq_(ensembl.exon_ids_of_transcript_id("ENSTTEST00000020001.7"), ["ENSETEST00000020001"])
    eq_(ensembl.transcript_ids_of_exon_id("ENSETEST00000020001.2"), ["ENSTTEST00000020001"])
    eq_(ensembl.locus_of_transcript_id("ENSTTEST00000020001.7").start, 100)


def test_versioned_ids_in_sequence_lookups(ensembl):
    eq_(ensembl.protein_sequence("ENSPTEST00000020001.5"), "MPKF")
    eq_(ensembl.protein_sequence("ENSPTEST00000020001"), "MPKF")
    eq_(ensembl.transcript_sequence("ENSTTEST00000020001.7"), "ATGCCCAAATTT")
    eq_(ensembl.transcript_sequence("ENSTTEST00000020001"), "ATGCCCAAATTT")


@pytest.mark.parametrize(
    "lookup, identifier, installed",
    [
        ("gene_by_id", "ENSGTEST00000020001.2", "ENSGTEST00000020001.3"),
        ("transcript_by_id", "ENSTTEST00000020001.8", "ENSTTEST00000020001.7"),
        ("exon_by_id", "ENSETEST00000020001.1", "ENSETEST00000020001.2"),
        ("gene_by_protein_id", "ENSPTEST00000020001.4", "ENSPTEST00000020001.5"),
        ("transcript_by_protein_id", "ENSPTEST00000020001.4", "ENSPTEST00000020001.5"),
        ("locus_of_gene_id", "ENSGTEST00000020001.4", "ENSGTEST00000020001.3"),
        ("protein_sequence", "ENSPTEST00000020001.4", "ENSPTEST00000020001.5"),
        ("transcript_sequence", "ENSTTEST00000020001.6", "ENSTTEST00000020001.7"),
    ],
)
def test_mismatched_version_names_the_installed_one(ensembl, lookup, identifier, installed):
    with pytest.raises(ValueError, match="%s is not in this annotation, which has %s" % (identifier, installed)):
        getattr(ensembl, lookup)(identifier)


def test_mismatch_is_not_cached_as_a_hit(ensembl):
    ensembl.transcript_by_id("ENSTTEST00000020001")
    with pytest.raises(ValueError, match="which has ENSTTEST00000020001.7"):
        ensembl.transcript_by_id("ENSTTEST00000020001.6")


def test_unknown_versioned_ids(ensembl):
    with pytest.raises(ValueError, match="not found"):
        ensembl.transcript_by_id("ENSTTEST00000099999.1")
    with pytest.raises(ValueError, match="not found"):
        ensembl.exon_by_id("ENSETEST00000099999.1")
    eq_(ensembl.protein_sequence("ENSPTEST00000099999.1"), None)


def test_annotation_without_versions_rejects_versioned_ids():
    with TemporaryDirectory() as tmpdir:
        genome = _ensembl_genome(tmpdir, gtf_versions=False, fasta_versions=False)
        eq_(genome.transcript_by_id("ENSTTEST00000020001").version, None)
        eq_(genome.protein_sequence("ENSPTEST00000020001"), "MPKF")
        message = "doesn't record versions for ENSTTEST00000020001"
        with pytest.raises(ValueError, match=message):
            genome.transcript_by_id("ENSTTEST00000020001.7")
        with pytest.raises(ValueError, match=message):
            genome.transcript_sequence("ENSTTEST00000020001.7")
        with pytest.raises(ValueError, match="doesn't record versions"):
            genome.protein_sequence("ENSPTEST00000020001.5")


def test_gencode_ids_accept_bare_and_check_versions():
    """GENCODE GTFs embed the version in each ID."""
    with TemporaryDirectory() as tmpdir:
        genome = _make_gencode_genome(tmpdir)
        genome.index()
        transcript = genome.transcript_by_id("ENSTTEST00000000001")
        eq_(transcript.id, "ENSTTEST00000000001.5")
        assert transcript is genome.transcript_by_id("ENSTTEST00000000001.5")
        eq_(genome.gene_by_protein_id("ENSPTEST00000000001").id, "ENSGTEST00000000001.4")
        eq_(genome.exon_by_id("ENSETEST00000000001").id, "ENSETEST00000000001.2")
        with pytest.raises(ValueError, match="which has ENSTTEST00000000001.5"):
            genome.transcript_by_id("ENSTTEST00000000001.4")
        with pytest.raises(ValueError, match="which has ENSPTEST00000000001.3"):
            genome.protein_sequence("ENSPTEST00000000001.2")


def test_tair_isoform_suffix_is_not_a_version():
    genome = Genome(
        reference_name="TAIR10",
        annotation_name="_test_versioned_ids_tair",
        cache_directory_path=custom_genome_cache_path("versioned-ids-tair"),
        gtf_path_or_url=data_path("arabidopsis.tair10.partial.gtf"),
        transcript_fasta_paths_or_urls=[data_path("arabidopsis.tair10.partial.cdna.fa")],
    )
    genome.index()
    transcript_id = genome.transcript_ids()[0]
    eq_(genome.transcript_by_id(transcript_id).id, transcript_id)
    stem, _, isoform = transcript_id.rpartition(".")
    other_isoform = "%s.%d" % (stem, int(isoform) + 100)
    with pytest.raises(ValueError, match="not found"):
        genome.transcript_by_id(other_isoform)
    eq_(genome.transcript_sequence(other_isoform), None)


def test_ensembl_81_checks_sequence_versions_against_the_gtf():
    """Ensembl 77-82 FASTA headers carry no versions; the GTF does."""
    mouse_genome.index()
    eq_(mouse_genome.protein_sequences.fasta_version("ENSMUSP00000099398"), None)
    sequence = mouse_genome.protein_sequence("ENSMUSP00000099398")
    eq_(mouse_genome.protein_sequence("ENSMUSP00000099398.3"), sequence)
    with pytest.raises(ValueError, match="which has ENSMUSP00000099398.3"):
        mouse_genome.protein_sequence("ENSMUSP00000099398.2")


def test_match_version_ambiguous_bare_id():
    with pytest.raises(ValueError, match="matches several versions"):
        match_version("ENSG1", {"ENSG1.1": 1, "ENSG1.2": 2})


def test_issue_445_examples_on_ensembl_93():
    eq_(grch38.protein_sequence("ENSP00000269305.4")[:10], "MEEPQSDPSV")
    for stale in ["ENSP00000269305.3", "ENSP00000269305.99"]:
        with pytest.raises(ValueError, match="which has ENSP00000269305.4"):
            grch38.protein_sequence(stale)
    eq_(grch38.transcript_by_id("ENST00000269305.8").name, "TP53-201")
    eq_(grch38.gene_by_id("ENSG00000141510.16").name, "TP53")
    eq_(grch38.transcript_by_protein_id("ENSP00000269305.4").id, "ENST00000269305")
    eq_(grch38.gene_by_protein_id("ENSP00000269305.4").id, "ENSG00000141510")

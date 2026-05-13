"""
Issue #335 (part 2): tolerate versioned protein / transcript IDs in FASTA
lookups for GENCODE-style genomes.

GENCODE embeds the assembly version directly in the ``protein_id`` and
``transcript_id`` GTF attributes (e.g. ``protein_id "ENSP00000123456.3"``),
while pyensembl's FASTA parser strips ``.N`` suffixes from ENS-prefix
headers. Before this fix the literal lookup with the versioned ID missed,
so ``Transcript.protein_sequence`` and ``Genome.transcript_sequence``
returned ``None`` for GENCODE FASTAs even when the sequence was present.
"""

import os

from pyensembl import Genome

from .common import TemporaryDirectory, eq_


GENCODE_STYLE_GTF = """\
1\ttest\ttranscript\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_name "FAKE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\texon\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; exon_number "1"; exon_id "ENSETEST00000000001.2"; gene_name "FAKE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tCDS\t100\t798\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; exon_number "1"; exon_id "ENSETEST00000000001.2"; protein_id "ENSPTEST00000000001.3"; gene_name "FAKE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tstart_codon\t100\t102\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_name "FAKE1"; protein_id "ENSPTEST00000000001.3"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tstop_codon\t799\t801\t.\t+\t0\tgene_id "ENSGTEST00000000001.4"; transcript_id "ENSTTEST00000000001.5"; gene_name "FAKE1"; protein_id "ENSPTEST00000000001.3"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
"""

# GENCODE-style FASTA: pipe-delimited header carrying the versioned ID.
PROTEIN_FASTA = (
    ">ENSPTEST00000000001.3|ENSTTEST00000000001.5|ENSGTEST00000000001.4|"
    "OTTHUMG00000000001.1|OTTHUMT00000000001.1|FAKE1-001|FAKE1|266\n"
    "MFAKEPROTEINSEQUENCE\n"
)
TRANSCRIPT_FASTA = (
    ">ENSTTEST00000000001.5|ENSGTEST00000000001.4|OTTHUMG00000000001.1|"
    "OTTHUMT00000000001.1|FAKE1-001|FAKE1|800|protein_coding|\n"
    "ATGCCCAAATTTGGGCCCAAATTT\n"
)


def _make_genome(tmpdir, with_transcript=True, with_protein=True):
    gtf_path = os.path.join(tmpdir, "gencode_style.gtf")
    with open(gtf_path, "w") as f:
        f.write(GENCODE_STYLE_GTF)
    transcript_fasta = None
    protein_fasta = None
    if with_transcript:
        transcript_fasta = os.path.join(tmpdir, "gencode_style.cdna.fa")
        with open(transcript_fasta, "w") as f:
            f.write(TRANSCRIPT_FASTA)
    if with_protein:
        protein_fasta = os.path.join(tmpdir, "gencode_style.pep.fa")
        with open(protein_fasta, "w") as f:
            f.write(PROTEIN_FASTA)
    return Genome(
        reference_name="GRCh38",
        annotation_name="_test_gencode_versioned_335",
        gtf_path_or_url=gtf_path,
        transcript_fasta_paths_or_urls=[transcript_fasta] if transcript_fasta else [],
        protein_fasta_paths_or_urls=[protein_fasta] if protein_fasta else [],
        cache_directory_path=tmpdir,
    )


def test_transcript_protein_sequence_with_gencode_versioned_id():
    """The CDS row stores ``protein_id "ENSPTEST00000000001.3"`` (with the
    GENCODE version suffix); the FASTA parser strips the suffix so the
    dict key is ``ENSPTEST00000000001``. ``Transcript.protein_sequence``
    must recover the sequence by stripping the version before lookup."""
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir)
        genome.index()
        transcript = genome.transcript_by_id("ENSTTEST00000000001.5")
        # protein_id is preserved verbatim from the GTF (with version)
        eq_(transcript.protein_id, "ENSPTEST00000000001.3")
        # but the sequence is still recovered
        eq_(transcript.protein_sequence, "MFAKEPROTEINSEQUENCE")


def test_transcript_sequence_with_gencode_versioned_id():
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir)
        genome.index()
        transcript = genome.transcript_by_id("ENSTTEST00000000001.5")
        eq_(transcript.sequence, "ATGCCCAAATTTGGGCCCAAATTT")


def test_genome_protein_sequence_accepts_versioned_id():
    """Genome.protein_sequence("ENSP....N") should resolve even when the
    FASTA dict keys are unversioned."""
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir, with_transcript=False)
        genome.index()
        eq_(
            genome.protein_sequence("ENSPTEST00000000001.3"),
            "MFAKEPROTEINSEQUENCE",
        )
        # unversioned form still works (existing behavior)
        eq_(
            genome.protein_sequence("ENSPTEST00000000001"),
            "MFAKEPROTEINSEQUENCE",
        )
        # bogus id returns None, not a KeyError or stripped-form hit
        eq_(genome.protein_sequence("ENSPNOTAREAL00000000.1"), None)


def test_genome_transcript_sequence_accepts_versioned_id():
    with TemporaryDirectory() as tmpdir:
        genome = _make_genome(tmpdir, with_protein=False)
        genome.index()
        eq_(
            genome.transcript_sequence("ENSTTEST00000000001.5"),
            "ATGCCCAAATTTGGGCCCAAATTT",
        )
        eq_(
            genome.transcript_sequence("ENSTTEST00000000001"),
            "ATGCCCAAATTTGGGCCCAAATTT",
        )

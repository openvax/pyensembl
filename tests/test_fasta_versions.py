"""
Issue #351: preserve FASTA-header versions in SequenceData instead of
stripping them at parse time.

Covers:
  * `_parse_header_id` retains versioned ENS IDs and handles GENCODE
    pipe-delimited headers.
  * `_parse_header_id` leaves non-ENS IDs (e.g. TAIR) alone.
  * `SequenceData` keys versioned IDs verbatim, builds a
    `_stripped_index` mapping bare ENS -> versioned, and exposes
    `fasta_version`.
  * `lookup_sequence_with_version_fallback` resolves both directions:
    versioned caller -> bare FASTA (Ensembl case) and bare caller ->
    versioned FASTA (GENCODE case).
"""
import os
import pickle
from os.path import join

from pyensembl import SequenceData
from pyensembl.common import dump_pickle
from pyensembl.fasta import _parse_header_id, _split_ens_version
from pyensembl.sequence_data import lookup_sequence_with_version_fallback

from .common import TemporaryDirectory, eq_


# -----------------------------
# _parse_header_id
# -----------------------------


def test_parse_header_id_keeps_ens_version_suffix():
    eq_(_parse_header_id(b">ENSP00000123456.3"), "ENSP00000123456.3")


def test_parse_header_id_bare_ens_unchanged():
    eq_(_parse_header_id(b">ENSP00000123456"), "ENSP00000123456")


def test_parse_header_id_splits_on_first_space():
    # Stock Ensembl FASTA header form
    eq_(
        _parse_header_id(b">ENSP00000123456.3 pep:protein_coding chromosome:GRCh38:1:1:1:1"),
        "ENSP00000123456.3",
    )


def test_parse_header_id_splits_on_first_pipe():
    # GENCODE FASTA header form
    eq_(
        _parse_header_id(
            b">ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|"
            b"OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|326"
        ),
        "ENSP00000493376.2",
    )


def test_parse_header_id_does_not_strip_tair_isoform_suffix():
    # TAIR .1 is an isoform, not a version — must not be touched.
    eq_(_parse_header_id(b">AT1G01010.1"), "AT1G01010.1")


# -----------------------------
# _split_ens_version
# -----------------------------


def test_split_ens_version_versioned():
    eq_(_split_ens_version("ENSP00000123456.3"), ("ENSP00000123456", 3))


def test_split_ens_version_bare():
    eq_(_split_ens_version("ENSP00000123456"), ("ENSP00000123456", None))


def test_split_ens_version_non_ens():
    # TAIR isoform suffix is not a version.
    eq_(_split_ens_version("AT1G01010.1"), ("AT1G01010.1", None))


def test_split_ens_version_non_integer_suffix():
    # Defensive: an unparseable suffix shouldn't blow up.
    eq_(_split_ens_version("ENSP00000123456.bogus"), ("ENSP00000123456.bogus", None))


# -----------------------------
# SequenceData
# -----------------------------


def _write_fasta(path, header_id, sequence):
    with open(path, "w") as f:
        f.write(">" + header_id + "\n" + sequence + "\n")


def test_sequence_data_keys_versioned_ens_id_verbatim():
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "v.fa")
        _write_fasta(fasta, "ENSPTEST00000001.3", "MFAKE")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # versioned key is the canonical one
        eq_(sd.get("ENSPTEST00000001.3"), "MFAKE")
        # bare ENS resolves via the stripped index
        eq_(sd._stripped_index["ENSPTEST00000001"], "ENSPTEST00000001.3")
        eq_(sd.fasta_version("ENSPTEST00000001.3"), 3)
        # fasta_version accepts the bare form too
        eq_(sd.fasta_version("ENSPTEST00000001"), 3)


def test_sequence_data_keys_bare_id_when_header_has_no_version():
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "bare.fa")
        _write_fasta(fasta, "ENSPTEST00000001", "MFAKE")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        eq_(sd.get("ENSPTEST00000001"), "MFAKE")
        # no version, no stripped_index entry, no fasta_version
        assert "ENSPTEST00000001" not in sd._stripped_index
        eq_(sd.fasta_version("ENSPTEST00000001"), None)


def test_sequence_data_tair_isoform_kept_as_is():
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "tair.fa")
        _write_fasta(fasta, "AT1G01010.1", "MTAIRSEQ")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # AT1G01010.1 is keyed as-is (no stripping at all)
        eq_(sd.get("AT1G01010.1"), "MTAIRSEQ")
        # AT1G01010 (without the .1) is NOT a valid lookup
        eq_(sd.get("AT1G01010"), None)
        # And fasta_version is None for it because TAIR isn't ENS-prefixed
        eq_(sd.fasta_version("AT1G01010.1"), None)


def test_sequence_data_stripped_index_keeps_highest_version_on_conflict():
    """If a FASTA contains two entries that bare-collide on an ENS ID
    with different version suffixes, the higher version is the canonical
    alias (version numbers are monotonically increasing)."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "dup.fa")
        with open(fasta, "w") as f:
            # write the older version first, then the newer
            f.write(">ENSPTEST00000001.1\n" + "MOLD\n")
            f.write(">ENSPTEST00000001.5\n" + "MNEW\n")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # both versioned forms reachable
        eq_(sd.get("ENSPTEST00000001.1"), "MOLD")
        eq_(sd.get("ENSPTEST00000001.5"), "MNEW")
        # bare form resolves to the highest version
        eq_(sd._stripped_index["ENSPTEST00000001"], "ENSPTEST00000001.5")


def test_old_bare_keyed_pickle_still_loads_under_new_code():
    """Caches written by older pyensembl (before this PR) keyed FASTA
    entries on the bare ENS ID. The new parser preserves versions, but
    a pre-existing pickle from the old code path must still load — its
    bare-keyed dict gets walked through `_add_to_fasta_dictionary` which
    leaves `_stripped_index` empty (the source dict has no versioned
    keys to index), and lookups via the version-fallback helper still
    resolve correctly because the version-stripped retry path handles
    the bare-cache case.
    """
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "legacy.fa")
        # write a versioned FASTA (so the new parser would key versioned)
        _write_fasta(fasta, "ENSPLEGACY00000001.4", "MLEGACY")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        # simulate an old (v1) bare-keyed pickle sitting on disk where
        # the new parser expects its cache file
        old_pickle_dict = {"ENSPLEGACY00000001": "MLEGACY"}
        dump_pickle(old_pickle_dict, sd.fasta_dictionary_pickle_paths[0])
        # index should pick up the existing pickle, NOT re-parse the FASTA
        sd.index()
        # bare lookup hits directly
        eq_(sd.get("ENSPLEGACY00000001"), "MLEGACY")
        # versioned lookup falls back via the strip-and-retry path
        eq_(
            lookup_sequence_with_version_fallback(sd, "ENSPLEGACY00000001.4"),
            "MLEGACY",
        )
        # _stripped_index stays empty for a bare-keyed cache
        eq_(sd._stripped_index, {})
        # fasta_version returns None for bare-only caches (the version
        # info was never captured at parse time)
        eq_(sd.fasta_version("ENSPLEGACY00000001.4"), None)
        eq_(sd.fasta_version("ENSPLEGACY00000001"), None)


def test_sequence_data_pickle_round_trip_rebuilds_stripped_index():
    """The pickle stores the FASTA dict; on reload the SequenceData
    rebuilds `_stripped_index` and `_versions` from the dict contents.
    Verify both versioned and bare lookups still work after a reload."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "rt.fa")
        _write_fasta(fasta, "ENSPTEST00000003.7", "MROUND")
        # first pass writes the pickle
        sd1 = SequenceData([fasta], cache_directory_path=tmpdir)
        sd1.index()
        # second SequenceData should pick up the pickle and still
        # populate _stripped_index from it
        sd2 = SequenceData([fasta], cache_directory_path=tmpdir)
        # trigger load
        eq_(sd2.get("ENSPTEST00000003.7"), "MROUND")
        eq_(sd2._stripped_index["ENSPTEST00000003"], "ENSPTEST00000003.7")
        eq_(sd2.fasta_version("ENSPTEST00000003"), 7)


# -----------------------------
# lookup_sequence_with_version_fallback
# -----------------------------


def test_lookup_direct_hit_versioned():
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "d.fa")
        _write_fasta(fasta, "ENSPTEST00000004.2", "MDIRECT")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        eq_(
            lookup_sequence_with_version_fallback(sd, "ENSPTEST00000004.2"),
            "MDIRECT",
        )


def test_lookup_versioned_caller_against_bare_fasta():
    """Caller has a versioned ID (e.g. GENCODE GTF protein_id with .N)
    but the FASTA on disk was the older bare style."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "b.fa")
        _write_fasta(fasta, "ENSPTEST00000005", "MBAREONLY")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        eq_(
            lookup_sequence_with_version_fallback(sd, "ENSPTEST00000005.9"),
            "MBAREONLY",
        )


def test_lookup_bare_caller_against_versioned_fasta():
    """Caller has a bare ID (e.g. Ensembl GTF protein_id w/o version)
    but the FASTA on disk was the GENCODE pipe-delimited form."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "v.fa")
        _write_fasta(fasta, "ENSPTEST00000006.2", "MVERSIONEDONLY")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # bare lookup must resolve via _stripped_index
        eq_(
            lookup_sequence_with_version_fallback(sd, "ENSPTEST00000006"),
            "MVERSIONEDONLY",
        )


def test_lookup_returns_none_for_unknown_id():
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "n.fa")
        _write_fasta(fasta, "ENSPTEST00000007.1", "MFAKE")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        eq_(lookup_sequence_with_version_fallback(sd, "ENSPNOTHERE.1"), None)
        eq_(lookup_sequence_with_version_fallback(sd, "ENSPNOTHERE"), None)
        eq_(lookup_sequence_with_version_fallback(sd, ""), None)
        eq_(lookup_sequence_with_version_fallback(sd, None), None)


def test_lookup_does_not_strip_tair_isoform_suffix():
    """Non-ENS .N suffixes are isoform identifiers, not versions, and
    must not be stripped during the version-fallback path."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "tair.fa")
        _write_fasta(fasta, "AT1G01010.1", "MTAIR")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # exact match works
        eq_(lookup_sequence_with_version_fallback(sd, "AT1G01010.1"), "MTAIR")
        # bare lookup must NOT silently succeed (the .1 is an isoform suffix)
        eq_(lookup_sequence_with_version_fallback(sd, "AT1G01010"), None)
        # AT1G01010.2 (asking for a different isoform) must NOT collide
        eq_(lookup_sequence_with_version_fallback(sd, "AT1G01010.2"), None)


def test_legacy_pickle_without_stripped_index_falls_back_cleanly():
    """If a SequenceData ends up loading a pickle that pre-dates this
    PR (no _stripped_index attribute), lookups should still resolve via
    the literal dict get. Belt-and-braces against an upgrade path."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "legacy.fa")
        _write_fasta(fasta, "ENSPLEGACY00000001.4", "MLEGACY")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # simulate an old pickle that only knows the bare key
        sd._fasta_dictionary = {"ENSPLEGACY00000001": "MLEGACY"}
        sd._stripped_index = None  # legacy state
        sd._versions = None
        eq_(
            lookup_sequence_with_version_fallback(sd, "ENSPLEGACY00000001.4"),
            "MLEGACY",
        )


# -----------------------------
# Pickle contents
# -----------------------------


def test_pickled_fasta_dictionary_uses_versioned_keys():
    """The pickle is the source of truth for what gets cached. The
    versioned form should make it onto disk so the stripped_index can
    be rebuilt from it after a process restart."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "p.fa")
        _write_fasta(fasta, "ENSPTEST00000008.5", "MPICKLED")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        with open(sd.fasta_dictionary_pickle_paths[0], "rb") as f:
            data = pickle.load(f)
        assert "ENSPTEST00000008.5" in data
        assert data["ENSPTEST00000008.5"] == "MPICKLED"


# -----------------------------
# Transcript.fasta_version / Protein.fasta_version
# -----------------------------


_VERSIONED_GTF = """\
1\ttest\ttranscript\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000010001"; gene_version "1"; transcript_id "ENSTTEST00000010001"; transcript_version "7"; gene_name "FV1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\texon\t100\t800\t.\t+\t.\tgene_id "ENSGTEST00000010001"; gene_version "1"; transcript_id "ENSTTEST00000010001"; transcript_version "7"; exon_number "1"; exon_id "ENSETEST00000010001"; exon_version "1"; gene_name "FV1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tCDS\t100\t798\t.\t+\t0\tgene_id "ENSGTEST00000010001"; gene_version "1"; transcript_id "ENSTTEST00000010001"; transcript_version "7"; exon_number "1"; exon_id "ENSETEST00000010001"; exon_version "1"; protein_id "ENSPTEST00000010001"; protein_version "5"; gene_name "FV1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tstart_codon\t100\t102\t.\t+\t0\tgene_id "ENSGTEST00000010001"; transcript_id "ENSTTEST00000010001"; protein_id "ENSPTEST00000010001"; gene_name "FV1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1\ttest\tstop_codon\t799\t801\t.\t+\t0\tgene_id "ENSGTEST00000010001"; transcript_id "ENSTTEST00000010001"; protein_id "ENSPTEST00000010001"; gene_name "FV1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
"""


def _build_versioned_genome(tmpdir, transcript_fasta=True, protein_fasta=True):
    """Build a tiny Genome whose GTF, cDNA FASTA, and protein FASTA all
    carry the same versioned ENS IDs."""
    from pyensembl import Genome

    gtf_path = join(tmpdir, "v.gtf")
    with open(gtf_path, "w") as f:
        f.write(_VERSIONED_GTF)
    cdna_path = None
    pep_path = None
    if transcript_fasta:
        cdna_path = join(tmpdir, "v.cdna.fa")
        # Header records version 7 — same as the GTF's transcript_version
        _write_fasta(cdna_path, "ENSTTEST00000010001.7", "ATGCCCAAATTTGGGCCCAAATTT")
    if protein_fasta:
        pep_path = join(tmpdir, "v.pep.fa")
        # Header records version 5 — same as the GTF's protein_version
        _write_fasta(pep_path, "ENSPTEST00000010001.5", "MFAKEPROTEIN")
    g = Genome(
        reference_name="GRCh38",
        annotation_name="_test_fasta_versions",
        gtf_path_or_url=gtf_path,
        transcript_fasta_paths_or_urls=[cdna_path] if cdna_path else [],
        protein_fasta_paths_or_urls=[pep_path] if pep_path else [],
        cache_directory_path=tmpdir,
    )
    g.index()
    return g


def test_transcript_fasta_version_returns_header_version():
    """`Transcript.fasta_version` returns the int the cDNA FASTA header
    carried (7), independent of the GTF's `transcript_version` (also 7
    in this fixture but distinct in concept)."""
    with TemporaryDirectory() as tmpdir:
        g = _build_versioned_genome(tmpdir)
        t = g.transcript_by_id("ENSTTEST00000010001")
        eq_(t.transcript_version, 7)  # from the GTF
        eq_(t.fasta_version, 7)  # from the FASTA header


def test_transcript_fasta_version_none_when_no_transcript_fasta():
    with TemporaryDirectory() as tmpdir:
        g = _build_versioned_genome(tmpdir, transcript_fasta=False)
        t = g.transcript_by_id("ENSTTEST00000010001")
        eq_(t.fasta_version, None)


def test_protein_fasta_version_returns_header_version():
    """`Transcript.protein.fasta_version` returns the int the protein
    FASTA header carried (5), independent of the GTF's `protein_version`."""
    with TemporaryDirectory() as tmpdir:
        g = _build_versioned_genome(tmpdir)
        t = g.transcript_by_id("ENSTTEST00000010001")
        eq_(t.protein.protein_version, 5)  # from the GTF
        eq_(t.protein.fasta_version, 5)  # from the FASTA header


def test_protein_fasta_version_none_when_no_protein_fasta():
    with TemporaryDirectory() as tmpdir:
        g = _build_versioned_genome(tmpdir, protein_fasta=False)
        t = g.transcript_by_id("ENSTTEST00000010001")
        # protein object is still constructed (protein_id exists in CDS row)
        # but fasta_version returns None because there's no protein FASTA
        assert t.protein is not None
        eq_(t.protein.fasta_version, None)


def test_protein_fasta_version_disagrees_with_gtf_when_files_diverge():
    """Mixing a fresh GTF (protein_version=5) with a stale FASTA
    (header version 3) — the two accessors disagree and downstream
    tools can detect the mismatch."""
    with TemporaryDirectory() as tmpdir:
        from pyensembl import Genome
        gtf_path = join(tmpdir, "v.gtf")
        with open(gtf_path, "w") as f:
            f.write(_VERSIONED_GTF)  # claims protein version 5
        pep_path = join(tmpdir, "v.pep.fa")
        # FASTA header records version 3 — different from the GTF's 5
        _write_fasta(pep_path, "ENSPTEST00000010001.3", "MFAKEPROTEIN")
        g = Genome(
            reference_name="GRCh38",
            annotation_name="_test_fasta_versions_disagree",
            gtf_path_or_url=gtf_path,
            protein_fasta_paths_or_urls=[pep_path],
            cache_directory_path=tmpdir,
        )
        g.index()
        t = g.transcript_by_id("ENSTTEST00000010001")
        eq_(t.protein.protein_version, 5)
        eq_(t.protein.fasta_version, 3)


def test_protein_constructed_without_genome_returns_none_fasta_version():
    """Direct `Protein(...)` construction outside a Genome shouldn't blow
    up if the caller asks for `fasta_version` — just returns None."""
    from pyensembl import Protein

    p = Protein(protein_id="ENSPSTANDALONE.1", protein_version=1)
    eq_(p.fasta_version, None)

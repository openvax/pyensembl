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
from pyensembl.fasta import _parse_header_id, _split_ens_version
from pyensembl.sequence_data import (
    FASTA_PICKLE_SCHEMA_VERSION,
    lookup_sequence_with_version_fallback,
)

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


def test_sequence_data_pickle_filename_includes_schema_version():
    """Pickles must be namespaced so the upgrade from the v1 layout
    (bare-keyed dict) to the v2 layout (versioned-keyed dict with
    stripped_index) doesn't load a stale cache."""
    with TemporaryDirectory() as tmpdir:
        fasta = join(tmpdir, "vc.fa")
        _write_fasta(fasta, "ENSPTEST00000002.4", "MFAKE")
        sd = SequenceData([fasta], cache_directory_path=tmpdir)
        sd.index()
        # the actual pickle path must encode the schema version so an
        # older cached pickle without the index is ignored, not loaded.
        for pickle_path in sd.fasta_dictionary_pickle_paths:
            assert FASTA_PICKLE_SCHEMA_VERSION in pickle_path, pickle_path
            assert os.path.exists(pickle_path)


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

"""
Test SequenceData object to make sure it's correctly parsing FASTA files
and that we're able to clear and regenrate its cached representation of
a FASTA dictionary
"""
from os.path import exists

from nose.tools import assert_raises
from pyensembl import SequenceData

from .common import TemporaryDirectory
from .data import data_path


FASTA_PATH = data_path("mouse.ensembl.81.partial.ENSMUSG00000017167.fa")


def test_sequence_type():
    with TemporaryDirectory() as tmpdir:
        seqs_dna = SequenceData(
            FASTA_PATH,
            cache_directory_path=tmpdir)
        seq = seqs_dna.get("ENSMUST00000138942")
        assert seq is not None, \
            "Failed to find sequence for ENSMUST00000138942"
        assert isinstance(seq, str), \
            "Wrong sequence type, expected %s but got %s" % (str, type(seq))

def test_check_ensembl_id():
    with TemporaryDirectory() as tmpdir:
        seqs = SequenceData(
            FASTA_PATH,
            require_ensembl_ids=True,
            cache_directory_path=tmpdir)
        with assert_raises(ValueError):
            seqs.get("WeirdID")

def test_missing_sequence():
    with TemporaryDirectory() as tmpdir:
        seqs = SequenceData(FASTA_PATH, cache_directory_path=tmpdir)
        seq = seqs.get("NotInFasta")
        assert seq is None, "Should get None back for missing sequence"

def test_clear_cache():
    with TemporaryDirectory() as tmpdir:
        seqs = SequenceData(FASTA_PATH, cache_directory_path=tmpdir)
        assert not seqs._fasta_dictionary, \
            "Expected _fasta_dictionary to load lazily"

        seqs._load_or_create_fasta_dictionary_pickle()
        assert len(seqs._fasta_dictionary) > 0, \
            "FASTA dictionary didn't get created"

        seqs.clear_cache()
        assert not seqs._fasta_dictionary, \
            "Expected FASTA dictionary to be empty after clear_cache()"
        assert not exists(seqs.fasta_dictionary_pickle_path), \
            "Cached pickle file should have been deleted"

        seqs._load_or_create_fasta_dictionary_pickle()
        assert exists(seqs.fasta_dictionary_pickle_path), \
            "Cached pickle file should have been created"

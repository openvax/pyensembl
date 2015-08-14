"""
Test SequenceData object to make sure it's correctly parsing FASTA files
and that we're able to clear and regenrate its cached representation of
a FASTA dictionary
"""
from os.path import exists

from nose.tools import raises
from pyensembl import SequenceData

from skbio import DNASequence
from .data import data_path

FASTA_PATH = data_path("mouse.ensembl.81.partial.ENSMUSG00000017167.fa")

def test_reverse_sequence():
    seqs = SequenceData(FASTA_PATH, sequence_type=DNASequence)
    seq = seqs.get("ENSMUST00000138942")
    assert len(seq) > 0
    # make sure returned sequence supports complement and reverse methods:
    assert len(seq.complement()) == len(seq)

def test_complement_sequence():
    seqs = SequenceData(FASTA_PATH, sequence_type=DNASequence)
    seq = seqs.get("ENSMUST00000138942")
    assert len(seq) > 0
    # make sure returned sequence supports complement and reverse methods:
    assert len(seq.reverse()) == len(seq)

@raises(ValueError)
def test_check_ensembl_id():
    seqs = SequenceData(FASTA_PATH, require_ensembl_ids=True)
    seqs.get("WeirdID")

def test_missing_sequence():
    seqs = SequenceData(FASTA_PATH)
    seq = seqs.get("NotInFasta")
    assert seq is None, "Should get None back for missing sequence"

def test_clear_cache():
    seqs = SequenceData(FASTA_PATH)
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




"""
Test SequenceData object to make sure it's correctly parsing FASTA files
and that we're able to clear and regenrate its cached representation of
a FASTA dictionary
"""
from os.path import exists
from tempfile import TemporaryDirectory

import pytest

from pyensembl import SequenceData

from .data import data_path


FASTA_PATH = data_path("mouse.ensembl.81.partial.ENSMUSG00000017167.fa")


@pytest.mark.parametrize("path_indices", [[], [0], [0, 1], [0, 0, 1]])
def test_sequence_data_as_dictionary_key(tmp_path, path_indices):
    paths = [tmp_path / "first.fa", tmp_path / "second.fa"]
    for path in paths:
        path.write_text(">sequence\nACGT\n")
    selected = [str(paths[index]) for index in path_indices]
    first = SequenceData(selected)
    reordered = SequenceData(list(reversed(selected)))

    assert first == reordered
    assert hash(first) == hash(reordered)
    assert {first: "cached"}[reordered] == "cached"
    assert len({first, reordered}) == 1


def test_sequence_data_path_multiplicity_distinguishes_keys(tmp_path):
    path = tmp_path / "sequence.fa"
    path.write_text(">sequence\nACGT\n")
    single = SequenceData(str(path))
    repeated = SequenceData([str(path), str(path)])

    assert single != repeated
    cached = {single: "once", repeated: "twice"}
    assert len(cached) == 2
    assert cached[SequenceData([str(path)])] == "once"
    assert cached[SequenceData([str(path), str(path)])] == "twice"


def test_sequence_type():
    with TemporaryDirectory() as tmpdir:
        seqs_dna = SequenceData([FASTA_PATH], cache_directory_path=tmpdir)
        seq = seqs_dna.get("ENSMUST00000138942")
        assert seq is not None, "Failed to find sequence for ENSMUST00000138942"
        assert isinstance(seq, str), "Wrong sequence type, expected %s but got %s" % (
            str,
            type(seq),
        )


def test_missing_sequence():
    with TemporaryDirectory() as tmpdir:
        seqs = SequenceData([FASTA_PATH], cache_directory_path=tmpdir)
        seq = seqs.get("NotInFasta")
        assert seq is None, "Should get None back for missing sequence"


def test_clear_cache():
    with TemporaryDirectory() as tmpdir:
        seqs = SequenceData([FASTA_PATH], cache_directory_path=tmpdir)
        assert not seqs._fasta_dictionary, "Expected _fasta_dictionary to load lazily"

        seqs._load_or_create_fasta_dictionary_pickle()
        assert len(seqs._fasta_dictionary) > 0, "FASTA dictionary didn't get created"

        seqs.clear_cache()
        assert (
            not seqs._fasta_dictionary
        ), "Expected FASTA dictionary to be empty after clear_cache()"
        for pickle_path in seqs.fasta_dictionary_pickle_paths:
            assert exists(
                pickle_path
            ), "Cached pickle file should be preserved by clear_cache()"

        seqs.delete_index_files()
        for pickle_path in seqs.fasta_dictionary_pickle_paths:
            assert not exists(
                pickle_path
            ), "Cached pickle file should be deleted by delete_index_files()"

        seqs._load_or_create_fasta_dictionary_pickle()
        for pickle_path in seqs.fasta_dictionary_pickle_paths:
            assert exists(pickle_path), "Cached pickle file should have been created"

from mock import Mock, patch
from nose.tools import assert_raises

from pyensembl import EnsemblRelease

# pylint: disable=no-value-for-parameter
# pylint and mocking don't go very well together. pylint complains,
# for example, about _test_fai_index not being called with the
# mock_pyfaidx parameter.

@patch('pyensembl.reference_transcripts.ReferenceTranscripts.local_fasta_path')
@patch('pyensembl.reference_transcripts.exists')
@patch('pyensembl.reference_transcripts.pyfaidx')
def _test_fai_index(mock_pyfaidx, mock_exists, mock_local_fasta_path,
                    fai_exists, force_index):
    """
    Depending on whether we already have a .fai file created by
    pyfaidx, and whether we're using the force flag for
    EnsemblRelease(...)._index, we may or may not expect write_fai
    to be called.

    Return True if it was called.
    """
    data = EnsemblRelease(54)

    # Skip the GTF stuff, since we're testing pyfaidx
    data.db._connect_if_exists = Mock(return_value=False)
    data.db._create_database = Mock()

    mock_exists.return_value = fai_exists
    data._index(force=force_index)

    return mock_pyfaidx.Fasta().faidx.write_fai.called


def test_fai_exists_index():
    called = _test_fai_index(fai_exists=True, force_index=False)
    assert not called, 'Expected no new .fai file'


def test_fai_not_exists_index():
    called = _test_fai_index(fai_exists=False, force_index=False)
    assert called, 'Expected a new .fai file'


def test_fai_exists_force_index():
    called = _test_fai_index(fai_exists=True, force_index=True)
    assert called, 'Expected a new .fai file'


def test_fai_not_exists_force_index():
    called = _test_fai_index(fai_exists=False, force_index=True)
    assert called, 'Expected a new .fai file'


@patch('pyensembl.reference_transcripts.ReferenceTranscripts.index')
def _test_db_index(mock_index, db_exists):
    """
    Return True if the GTF database gets created, which should
    be different depending on whether the database already existed.

    Note: we need to mock the reference transcript indexing, as we're
    testing GTF indexing.
    """
    data = EnsemblRelease(54)
    data.db._connect_if_exists = Mock(return_value=db_exists)
    data.db._create_database = Mock()

    data._index(force=False)

    return data.db._create_database.called


def test_db_not_exists_index():
    called = _test_db_index(db_exists=False)
    assert called, "Expected a new database"


def test_db_exists_index():
    called = _test_db_index(db_exists=True)
    assert not called, "Expected no new database"


def test_auto_download_off():
    data = EnsemblRelease(54)
    data.db._connect_if_exists = Mock(return_value=False)
    assert_raises(ValueError, data.gene_names_at_locus,
        contig=6, position=29945884)

from __future__ import absolute_import

from mock import Mock, patch
from nose.tools import assert_raises

from pyensembl import EnsemblRelease

# pylint: disable=no-value-for-parameter
# pylint and mocking don't go very well together. pylint complains,
# for example, about _test_db_index not being called with the
# mock_index parameter.

@patch('pyensembl.sequence_data.SequenceData.index')
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
    data.index(force=False)

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

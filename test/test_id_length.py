from __future__ import absolute_import

from .common import major_releases

from nose.tools import nottest

@nottest
def check_id_length(method_name):
    for release in major_releases:
        method = getattr(release, method_name)
        # only load chromosome Y to speed up tests
        idents = method(contig="Y")
        assert len(idents) > 0, "No values returned by %s" % method_name
        assert all(len(ident) == 15 for ident in idents), \
           "Invalid IDs for %s: %s" % (
                method_name,
                [ident for ident in idents if len(ident) != 15])

def test_gene_id_length():
    check_id_length('gene_ids')

def test_transcript_id_length():
    check_id_length('transcript_ids')

def test_protein_id_length():
    check_id_length('protein_ids')

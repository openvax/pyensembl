"""
Tests for methods which return collections of transcript IDs that aren't
converting from some type of name or ID.
"""

from test_common import test_ensembl_releases
from pyensembl import EnsemblRelease

# subset of transcript IDs for HLA-A
HLA_A_TRANSCRIPT_IDS = [
        'ENST00000396634',
        'ENST00000376809',
        'ENST00000376806',
        'ENST00000376802',
        'ENST00000496081',
        'ENST00000495183',
        'ENST00000461903',
        'ENST00000479320',
]


def test_transcript_ids_ensembl77_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    ensembl = EnsemblRelease(77, auto_download=True)
    transcript_ids = ensembl.transcript_ids_at_locus(6, 29941260, 29945884)
    for transcript_id in HLA_A_TRANSCRIPT_IDS:
        assert transcript_id in transcript_ids, \
            "Transcript %s of HLA-A not found overlapping locus" % transcript_id

KNOWN_TRANSCRIPT_IDS = HLA_A_TRANSCRIPT_IDS + [
    'ENST00000398417', # transcript ID of SMAD4-001
    'ENST00000334701', # transcript ID of HSP90AA1-001
    'ENST00000599837', # transcript ID of CTAG1A-002
]

# TODO: add release 54 after transcript IDs for older GTFs are filled in
# See https://github.com/hammerlab/pyensembl/issues/20
@test_ensembl_releases(75, 77)
def test_all_transcript_ids(ensembl):
    transcript_ids = set(ensembl.transcript_ids())
    for transcript_id in KNOWN_TRANSCRIPT_IDS:
        assert transcript_id in transcript_ids, \
            "Missing transcript ID %s from %s" % (transcript_id, ensembl)


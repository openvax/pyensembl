from pyensembl import EnsemblRelease

from nose.tools import ok_

ensembl = EnsemblRelease(77, auto_download=True)

def test_transcript_id_of_protein_id_release77():
    transcript_id = ensembl.transcript_id_of_protein_id(
        "ENSP00000485678")
    ok_('ENST00000623976', transcript_id)


from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(77, auto_download=True)
# chr6:29,945,884  is a position for HLA-A
# based on:
# http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884

def test_transcript_ids_hla_a():
    transcript_ids = ensembl.transcript_ids_at_locus(6, 29941260, 29945884)
    expected = [
        'ENST00000396634',
        'ENST00000376809',
        'ENST00000376806',
        'ENST00000620168',
        'ENST00000376802',
        'ENST00000496081',
        'ENST00000495183',
        'ENST00000461903',
        'ENST00000479320',
    ]
    for x in expected:
        assert x in transcript_ids, "Expected transcript: %s" % x

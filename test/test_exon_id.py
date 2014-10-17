from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(77)

TP53_EXON_IDS_RELEASE_77 = [
    'ENSE00002337729', 'ENSE00002419584',
    'ENSE00003625790', 'ENSE00003518480',
    'ENSE00003723991', 'ENSE00003712342',
    'ENSE00001657961', 'ENSE00003725258',
    'ENSE00003740946', 'ENSE00002204316',
    'ENSE00002064269', 'ENSE00003750554',
    'ENSE00003634848', 'ENSE00003492844',
    'ENSE00003735852', 'ENSE00003545950',
    'ENSE00003605891', 'ENSE00002051192',
    'ENSE00002084733', 'ENSE00003726882',
    'ENSE00001146308', 'ENSE00002667911',
    'ENSE00003752869', 'ENSE00003739898',
    'ENSE00003753508', 'ENSE00002034209',
    'ENSE00002030826', 'ENSE00001596491',
    'ENSE00002037735', 'ENSE00003736616',
    'ENSE00002672443', 'ENSE00002226620',
    'ENSE00003715195', 'ENSE00003750794',
    'ENSE00003745267', 'ENSE00003746220',
    'ENSE00003656695', 'ENSE00003669712',
    'ENSE00002051873', 'ENSE00002048269',
    'ENSE00002670535', 'ENSE00002677565',
    'ENSE00003532881', 'ENSE00003520683',
    'ENSE00002076714', 'ENSE00002062958',
    'ENSE00002073243', 'ENSE00003670707',
    'ENSE00002065802', 'ENSE00002362269'
]

def test_exon_ids_of_gene_id():
    exon_ids = ensembl.exon_ids_of_gene_id('ENSG00000141510')
    assert len(exon_ids) == len(TP53_EXON_IDS_RELEASE_77), \
        "Wrong number of exons, expected %d but got %d (n_distinct=%d)" % (
            len(TP53_EXON_IDS_RELEASE_77),
            len(exon_ids),
            len(set(exon_ids))
        )
    assert all(exon_id in TP53_EXON_IDS_RELEASE_77 for exon_id in exon_ids)

def test_exon_ids_of_gene_name():
    exon_ids = ensembl.exon_ids_of_gene_name("TP53")
    assert len(exon_ids) == len(TP53_EXON_IDS_RELEASE_77), \
        "Wrong number of exons, expected %d but got %d (n_distinct=%d)" % (
            len(TP53_EXON_IDS_RELEASE_77),
            len(exon_ids),
            len(set(exon_ids))
        )
    assert all(exon_id in TP53_EXON_IDS_RELEASE_77 for exon_id in exon_ids)

TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77 = [
    'ENSE00002064269',
    'ENSE00003723991',
    'ENSE00003712342',
    'ENSE00003725258',
    'ENSE00003740946',
    'ENSE00003750554',
    'ENSE00003634848',
    'ENSE00003492844'
]

def test_exon_ids_of_transcript_name():
    exon_ids = ensembl.exon_ids_of_transcript_name("TP53-026")
    assert len(exon_ids) == len(TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77), \
        "Expected %d exons, got %d"  % (
            len(TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77),
            len(exon_ids)
        )
    assert all(
            exon_id in TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77
            for exon_id in exon_ids
    )

def exon_ids_of_transcript_id():
    exon_ids = ensembl.exon_ids_of_transcript_id("ENST00000610623")
    assert len(exon_ids) == len(TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77), \
        "Expected %d exons, got %d"  % (
            len(TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77),
            len(exon_ids)
        )
    assert all(
            exon_id in TP53_TRANSCRIPT_26_EXON_IDS_RELEASE_77
            for exon_id in exon_ids
    )
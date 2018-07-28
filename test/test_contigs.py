from pyensembl import genome_for_reference_name

grch38 = genome_for_reference_name("GRCh38")

def test_contig_names():
    contig_names = set(grch38.contigs())
    for chrom in list(range(1, 23)) + ["X", "Y", "MT"]:
        assert str(chrom) in contig_names, (chrom, contig_names)

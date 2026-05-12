"""
Issue #190: Genome.merged_gene_intervals(contig) returns sorted
non-overlapping (start, end) tuples after merging adjacent/overlapping
gene loci on the contig.
"""

from pyensembl.common import merge_intervals

from .common import TemporaryDirectory, eq_
from pyensembl import Genome


def test_merge_intervals_handles_overlap_adjacency_and_disjoint():
    # Overlapping
    eq_(merge_intervals([(1, 10), (5, 15)]), [(1, 15)])
    # Adjacent (end+1 == next start) merges
    eq_(merge_intervals([(1, 10), (11, 20)]), [(1, 20)])
    # One-base gap stays separate
    eq_(merge_intervals([(1, 10), (12, 20)]), [(1, 10), (12, 20)])
    # Unsorted input is reordered
    eq_(merge_intervals([(50, 60), (1, 10), (5, 15)]), [(1, 15), (50, 60)])
    # Containment
    eq_(merge_intervals([(1, 100), (10, 20)]), [(1, 100)])
    # Empty
    eq_(merge_intervals([]), [])


def test_merged_gene_intervals_against_synthetic_gtf():
    """Build a tiny GTF with three overlapping and one disjoint gene on
    chr1 plus a single gene on chr2, then assert the merged intervals."""
    with TemporaryDirectory() as tmpdir:
        import os
        gtf_path = os.path.join(tmpdir, "merge.gtf")
        with open(gtf_path, "w") as f:
            # chr1: gene1 [100,200], gene2 [150,250] overlap → [100,250]
            #       gene3 [251,300] adjacent → merges into [100,300]
            #       gene4 [500,600] disjoint
            # chr1 -: gene5 on minus strand for strand-filter test
            # chr2: gene6 standalone
            for gene_id, contig, start, end, strand in [
                ("G1", "1", 100, 200, "+"),
                ("G2", "1", 150, 250, "+"),
                ("G3", "1", 251, 300, "+"),
                ("G4", "1", 500, 600, "+"),
                ("G5", "1", 700, 800, "-"),
                ("G6", "2", 1000, 1100, "+"),
            ]:
                f.write(
                    "%s\ttest\tgene\t%d\t%d\t.\t%s\t.\t"
                    'gene_id "%s"; gene_name "%s";\n'
                    % (contig, start, end, strand, gene_id, gene_id)
                )
                f.write(
                    "%s\ttest\ttranscript\t%d\t%d\t.\t%s\t.\t"
                    'gene_id "%s"; transcript_id "%s_t"; gene_name "%s";\n'
                    % (contig, start, end, strand, gene_id, gene_id, gene_id)
                )
                f.write(
                    "%s\ttest\texon\t%d\t%d\t.\t%s\t.\t"
                    'gene_id "%s"; transcript_id "%s_t"; exon_number "1"; gene_name "%s";\n'
                    % (contig, start, end, strand, gene_id, gene_id, gene_id)
                )

        genome = Genome(
            reference_name="GRCh38",
            annotation_name="merge_intervals_test",
            gtf_path_or_url=gtf_path,
            cache_directory_path=tmpdir,
        )
        genome.index()

        # chr1 (any strand): overlapping triplet collapses, gene4 stays,
        # minus-strand gene5 is included when strand is not filtered.
        eq_(
            genome.merged_gene_intervals("1"),
            [(100, 300), (500, 600), (700, 800)],
        )
        # Strand-restricted: gene5 is excluded
        eq_(
            genome.merged_gene_intervals("1", strand="+"),
            [(100, 300), (500, 600)],
        )
        # Only the minus-strand gene
        eq_(
            genome.merged_gene_intervals("1", strand="-"),
            [(700, 800)],
        )
        # chr2: single gene
        eq_(genome.merged_gene_intervals("2"), [(1000, 1100)])
        # Empty contig returns []
        eq_(genome.merged_gene_intervals("3"), [])

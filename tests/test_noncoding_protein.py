"""An annotation of only noncoding RNAs has no canonical protein table (#459)."""

import pytest

from pyensembl import Genome


@pytest.mark.parametrize("include_cds", [False, True])
def test_protein_accessors_without_protein_id(tmp_path, include_cds):
    gtf = tmp_path / "noncoding.gtf"
    attrs = ('gene_id "g"; transcript_id "t"; gene_name "G"; '
             'gene_biotype "lncRNA"; transcript_biotype "lncRNA";')
    rows = ["1\ttest\t%s\t100\t108\t.\t+\t.\t%s\n" % (feature, attrs)
            for feature in ("gene", "transcript")]
    rows.append('1\ttest\texon\t100\t108\t.\t+\t.\t%s exon_id "e"; exon_number "1";\n' % attrs)
    if include_cds:
        rows.append("1\ttest\tCDS\t100\t105\t.\t+\t0\t%s\n" % attrs)
    gtf.write_text("".join(rows))
    with Genome(reference_name="synthetic", annotation_name="noncoding",
                gtf_path_or_url=str(gtf), cache_directory_path=str(tmp_path)) as genome:
        genome.index()
        transcript = genome.transcript_by_id("t")
        assert transcript.protein_id is None
        assert transcript.protein is None
        assert transcript.protein_sequence is None
        assert not transcript.is_protein_coding

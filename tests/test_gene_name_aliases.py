"""Alias lookup is explicit, ambiguity-preserving, and annotation-scoped."""

import gzip

import pytest

from pyensembl import EnsemblRelease, GeneNameAliases, Genome


# Deliberately add a second gene to the real p53 alias to exercise ambiguity.
HGNC = (
    'symbol\talias_symbol\tprev_symbol\tensembl_gene_id\tstatus\n'
    'TP53\tp53|LFS1\tTRP53\tENSG00000141510\tApproved\n'
    'BRAF\tp53\tBRAF1\tENSG00000157764\tApproved\n'
    'MISSING\tp53\t\tENSG00000999999\tApproved\n'
    'WITHDRAWN\tp53\t\tENSG00000100000\tEntry Withdrawn\n'
    'NO_ID\tNO_GENE\t\t\tApproved\n'
)


@pytest.fixture
def genome(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000141510.7"; gene_name "TP53";\n'
        '1\ttest\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG00000157764"; gene_name "BRAF";\n'
        '1\ttest\tgene\t40\t50\t.\t+\t.\tgene_id "AT1G01010.1"; gene_name "NAC001";\n'
    )
    with Genome("synthetic", "aliases", gtf_path_or_url=str(gtf),
                cache_directory_path=str(tmp_path / "cache")) as result:
        result.index()
        yield result


@pytest.mark.parametrize("compressed", [False, True])
def test_hgnc_alias_lookup_preserves_ambiguity_and_versions(genome, tmp_path, compressed):
    path = tmp_path / ("hgnc.tsv.gz" if compressed else "hgnc.tsv")
    path.write_bytes(gzip.compress(HGNC.encode()) if compressed else HGNC.encode())
    aliases = GeneNameAliases.from_hgnc(path)
    assert aliases.species == "homo_sapiens"
    assert aliases.source == str(path)
    assert "WITHDRAWN" not in aliases
    assert "NO_GENE" not in aliases
    assert genome.gene_ids_of_gene_name("p53", aliases=aliases) == [
        "ENSG00000141510.7", "ENSG00000157764"
    ]
    assert [gene.name for gene in genome.genes_by_name("TRP53", aliases=aliases)] == ["TP53"]
    assert genome.gene_ids_of_gene_name("TP53", aliases=aliases) == ["ENSG00000141510.7"]
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("P53", aliases=aliases)
    with pytest.raises(ValueError):
        genome.genes_by_name("p53")
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("MISSING", aliases=aliases)


def test_exact_and_alias_hits_are_combined_without_duplicates(genome):
    aliases = {"TP53": ["ENSG00000157764", "ENSG00000141510", "ENSG00000157764"]}
    assert genome.gene_ids_of_gene_name("TP53", aliases=aliases) == [
        "ENSG00000141510.7", "ENSG00000157764"
    ]
    assert genome.gene_ids_of_gene_name("TP53") == ["ENSG00000141510.7"]


def test_custom_aliases_preserve_non_ensembl_identifier_suffixes(genome):
    assert genome.gene_ids_of_gene_name("NAC", aliases={"NAC": "AT1G01010.1"}) == ["AT1G01010.1"]
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("NAC", aliases={"NAC": "AT1G01010"})


def test_aliases_do_not_require_gene_names_in_the_gtf(tmp_path):
    source = tmp_path / "no-names.gtf"
    source.write_text('1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n')
    with Genome("synthetic", "no-names", gtf_path_or_url=str(source),
                cache_directory_path=str(tmp_path / "cache")) as genome:
        genome.index()
        assert genome.genes_by_name("named", aliases={"named": ["g1"]})[0].id == "g1"


def test_human_aliases_rejected_for_known_nonhuman_releases(tmp_path):
    path = tmp_path / "hgnc.tsv"
    path.write_text(HGNC)
    with pytest.raises(ValueError, match="homo_sapiens cannot be used with mus_musculus"):
        EnsemblRelease(93, species="mouse").genes_by_name(
            "p53", aliases=GeneNameAliases.from_hgnc(path)
        )


def test_mapping_copies_inputs():
    original = {"alias": ["g1"]}
    aliases = GeneNameAliases(original)
    original["alias"].append("g2")
    assert aliases["alias"] == ("g1",)
    with pytest.raises(TypeError):
        aliases["alias"] = ["g3"]


@pytest.mark.parametrize("aliases", [{"": "g1"}, {"alias": None}, {"alias": [1]}, {"alias": ""}])
def test_invalid_mapping_rejected(aliases):
    with pytest.raises(ValueError):
        GeneNameAliases(aliases)


@pytest.mark.parametrize("text", ["", "symbol\talias_symbol\n", HGNC.splitlines()[0] + "\nTP53\n",
                                  HGNC.replace("ENSG00000141510", "not-an-ensembl-id")])
def test_invalid_hgnc_data_rejected(tmp_path, text):
    path = tmp_path / "invalid.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        GeneNameAliases.from_hgnc(path)

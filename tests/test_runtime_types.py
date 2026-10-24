"""Preserve public validation and cached queries without the typechecks package."""

import numpy as np
import pytest

from pyensembl import Genome, normalization
from .test_legacy_cache_compatibility import snapshot


class Integer(int):
    pass


class String(str):
    pass


@pytest.fixture(autouse=True)
def fresh_normalization_cache(monkeypatch):
    monkeypatch.setattr(normalization, "NORMALIZE_CHROMOSOME_CACHE", {})


@pytest.mark.parametrize("warm", [False, True])
@pytest.mark.parametrize("invalid", [1.0, np.int64(1), None, b"1", [], {}])
def test_chromosome_type_rejection_is_independent_of_cache(invalid, warm):
    if warm:
        assert normalization.normalize_chromosome(1) == "1"
    with pytest.raises(TypeError):
        normalization.normalize_chromosome(invalid)


@pytest.mark.parametrize("value,expected", [
    (1, "1"), (Integer(1), "1"), (True, "TRUE"), (False, "FALSE"),
    (String("x"), "X"), ("chrx", "chrX"), ("chrUn_gl000212", "chrUn_gl000212"),
])
def test_normalization_preserves_accepted_values_and_memoization(value, expected):
    first = normalization.normalize_chromosome(value)
    assert first == expected
    assert normalization.normalize_chromosome(value) is first


@pytest.mark.parametrize("value", [0, "0", ""])
def test_invalid_chromosome_values_remain_value_errors(value):
    with pytest.raises(ValueError):
        normalization.normalize_chromosome(value)


@pytest.mark.parametrize("value,expected", [("+", "+"), ("-1", "-"), (1, "+"), (-1, "-"), (True, "+"), (1.0, "+")])
def test_strand_value_semantics_are_unchanged(value, expected):
    assert normalization.normalize_strand(value) == expected


@pytest.fixture
def database(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache))
    gtf = tmp_path / "tiny.gtf"
    gtf.write_text('1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g"; gene_name "G";\n')
    genome = Genome(reference_name="test", annotation_name="types", gtf_path_or_url=str(gtf))
    genome.index()
    genome.close()
    before = snapshot(cache)
    try:
        yield genome.db
    finally:
        genome.close()
        assert snapshot(cache) == before


@pytest.mark.parametrize("field,value,error", [
    ("column_name", "", ValueError), ("column_name", 1, TypeError),
    ("column_name", None, TypeError), ("column_name", b"gene_id", TypeError),
    ("position", "1", TypeError), ("position", None, TypeError),
    ("position", 1.0, TypeError), ("position", np.int64(1), TypeError),
    ("end", "2", TypeError), ("end", 2.0, TypeError),
    ("end", np.int64(2), TypeError),
])
def test_locus_query_preserves_validation_errors(database, field, value, error):
    options = dict(column_name="gene_id", feature="gene", contig="1", position=1)
    options[field] = value
    with pytest.raises(error, match=field):
        database.column_values_at_locus(**options)


@pytest.mark.parametrize("coordinate", [1, Integer(1), True])
def test_locus_query_keeps_integer_subclasses_and_typed_cells(database, coordinate):
    assert database.column_values_at_locus(String("gene_id"), "gene", Integer(1), coordinate) == ["g"]
    assert database.column_values_at_locus("start", "gene", "1", coordinate, end=coordinate) == [1]
    assert database.distinct_column_values_at_locus("gene_name", "gene", "1", coordinate) == ["G"]


def test_locus_query_rejects_invalid_contig_after_cached_integer(database):
    assert database.column_values_at_locus("gene_id", "gene", 1, 1) == ["g"]
    with pytest.raises(TypeError):
        database.column_values_at_locus("gene_id", "gene", 1.0, 1)

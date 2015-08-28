from __future__ import absolute_import
from os.path import exists

from pyensembl import ensembl_grch37 as ensembl
import pandas as pd
from nose.tools import eq_

from .common import test_ensembl_releases

@test_ensembl_releases()
def gtf_path_endswith_gtf_gz(ensembl):
    path = ensembl.gtf.gtf_path
    assert exists(path)
    assert path.endswith(".gtf.gz")

def test_dataframe_column_at_locus():
    series = ensembl.gtf.dataframe_column_at_locus(
        column_name='gene_id',
        contig='1',
        start=11869)
    assert series is not None
    assert len(series) > 0
    eq_('ENSG00000223972', series.unique()[0])

def test_dataframe_at_locus():
    df = ensembl.gtf.dataframe_at_locus(contig='1',
                                        start=11869)
    assert df is not None
    assert len(df) > 0
    eq_('ENSG00000223972', df['gene_id'].unique()[0])

def test_slice_column():
    column_name_series = pd.Series(["a", "b", "c"])
    start_series = pd.Series([1000, 2000, 3000])
    end_series = pd.Series([4000, 2300, 3500])

    series = ensembl.gtf._slice_column(
        column_name_series, start_series, end_series,
        900, 1000)
    eq_({'a'}, set(series))

    series = ensembl.gtf._slice_column(
        column_name_series, start_series, end_series,
        1000, 5000)
    eq_({'a', 'b', 'c'}, set(series))

    series = ensembl.gtf._slice_column(
        column_name_series, start_series, end_series,
        1900, 2000)
    eq_({'a', 'b'}, set(series))

    series = ensembl.gtf._slice_column(
        column_name_series, start_series, end_series,
        1900, 3000)
    eq_({'a', 'b', 'c'}, set(series))

    series = ensembl.gtf._slice_column(
        column_name_series, start_series, end_series,
        3000, 6000)
    eq_({'a', 'c'}, set(series))

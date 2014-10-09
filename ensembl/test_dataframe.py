from ensembl.human_data import HumanData

def test_df_length():
    ensembl = HumanData(75)
    df = ensembl.dataframe()
    assert df is not None
    assert len(df) == 2828312
    assert len(df.columns) == 9

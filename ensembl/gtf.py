from os.path import exists

import pandas as pd

import cython

GTF_COLS = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]

def _read_gtf(filename):
    """
    Download GTF annotation data for given release (if not already present),
    parse it into a DataFrame, leave the attributes in a single string column
    """
    assert exists(filename)
    assert filename.endswith(".gtf") or filename.endswith(".gtf.gz"), \
        "Must be a GTF file (%s)" % filename
    compression = "gzip" if filename.endswith(".gz") else None
    df = pd.read_csv(
        filename,
        comment='#',
        sep='\t',
        names = GTF_COLS,
        na_filter=False,
        compression = compression)
    df['seqname'] = df['seqname'].map(str)
    return df

def _attribute_dictionaries(df):
    """
    Given a DataFrame with attributes in semi-colon
    separated string, split them apart into per-entry dictionaries
    """
    # this crazy tripe-nested comprehensions
    # is unfortunately the best way to parse
    # ~2 million attribute strings while still using
    # vaguely Pythonic code
    attr_dicts = [
        {
            key : value.replace("\"", "").replace(";", "")
            for (key,value) in (
                pair_string.split(" ", 2)
                for pair_string in attr_string.split("; ")
            )
        }
        for attr_string in df.attribute
    ]
    return attr_dicts

def _extend_with_attributes(df):
    n = len(df)
    extra_columns = {}
    column_order = []
    for i, attr_string in enumerate(df.attribute):
        pairs = (
            kv.strip().split(" ", 2)
            for kv in attr_string.split(";")
            if len(kv) > 0
        )
        for k,v in pairs:
            if k not in extra_columns:
                extra_columns[k] = [""] * n
                column_order.append(k)
            extra_columns[k][i] = v
    # make a copy of the DataFrame without attribute strings.
    df = df.drop("attribute", axis=1)

    print "Adding attribute columns: %s" % (column_order,)
    for k in column_order:
        df[k] = extra_columns[k]
        # delete from dictionary since these objects are big
        # and we might be runniung low on memory
        del extra_columns[k]
        # remove quotes around values
        df[k] = df[k].str.replace("\"", "")
    return df

def load_gtf_as_dataframe(filename):
    print "Reading GTF %s into DataFrame" % filename
    df = _read_gtf(filename)
    print "Extracting attributes for %d entries in GTF DataFrame" % len(df)
    df = _extend_with_attributes(df)
    return df

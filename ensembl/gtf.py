from os.path import exists
import pandas as pd

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

def load_gtf_as_dataframe(filename):
    """
    Download GTF annotation data for given release (if not already present),
    parse it into a dataframe.
    """
    assert exists(filename)
    assert filename.endswith(".gtf") or filename.endswith(".gtf.gz"), \
        "Must be a GTF file (%s)" % filename
    compression = "gzip" if filename.endswith(".gz") else None
    return pd.read_csv(
        filename,
        comment='#',
        sep='\t',
        names = GTF_COLS,
        compression = compression)

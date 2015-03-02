# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Parse an Ensembl GTF file into a dataframe, expanding all of its attributes
into their own columns,


Columns of a GTF file:

    seqname   - name of the chromosome or scaffold; chromosome names
                without a 'chr'
    source    - name of the program that generated this feature, or
                the data source (database or project name)
    feature   - feature type name. Current allowed features are
                {gene, transcript, exon, CDS, Selenocysteine, start_codon,
                stop_codon and UTR}
    start     - start position of the feature, with sequence numbering
                starting at 1.
    end       - end position of the feature, with sequence numbering
                starting at 1.
    score     - a floating point value indiciating the score of a feature
    strand    - defined as + (forward) or - (reverse).
    frame     - one of '0', '1' or '2'. Frame indicates the number of base pairs
                before you encounter a full codon. '0' indicates the feature
                begins with a whole codon. '1' indicates there is an extra
                base (the 3rd base of the prior codon) at the start of this feature.
                '2' indicates there are two extra bases (2nd and 3rd base of the
                prior exon) before the first codon. All values are given with
                relation to the 5' end.
    attribute - a semicolon-separated list of tag-value pairs (separated by a space),
                providing additional information about each feature. A key can be
                repeated multiple times.

(from ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/README)
"""

from __future__ import print_function, division, absolute_import
import logging
from os.path import exists

import numpy as np
import pandas as pd


GTF_COLS = [
    'seqname',
    # Different versions of GTF use this column as
    # of: (1) gene biotype (2) transcript biotype or (3) the annotation source
    #
    # See: https://www.biostars.org/p/120306/#120321
    'second_column',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]

def pandas_series_is_biotype(series):
    """
    Hackishly infer whether a Pandas series is either from
    gene_biotype or transcript_biotype annotations by checking
    whether the 'protein_coding' biotype annotation is among its values.
    """
    return 'protein_coding' in series.values

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
        names=GTF_COLS,
        na_filter=False,
        compression=compression)

    df['seqname'] = df['seqname'].map(str)

    # very old GTF files use the second column to store the gene biotype
    # others use it to store the transcript biotype and
    # anything beyond release 77+ will use it to store
    # the source of the annotation (e.g. "havana")
    if pandas_series_is_biotype(df['second_column']):
        # patch this later to either 'transcript_biotype' or 'gene_biotype'
        # depending on what other annotations are present
        column_name = 'biotype'
    else:
        column_name = 'source'
    df.rename(columns={'second_column' : column_name}, inplace=True)

    return df

def _attribute_dictionaries(df):
    """
    Given a DataFrame with attributes in semi-colon
    separated string, split them apart into per-entry dictionaries
    """
    # this crazy triply-nested comprehensions
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
        # TODO: tokenize attribute strings properly
        # for now, catch mistaken semicolons by replacing "xyz;" with "xyz"
        # Required to do this since the GTF for Ensembl 78 has
        # gene_name = "PRAMEF6;"
        # transcript_name = "PRAMEF6;-201"
        attr_string = attr_string.replace(';\"', '\"').replace(";-", "-")

        pairs = (
            kv.strip().split(" ", 2)
            for kv in attr_string.split(";")
            # simplest entry must be at least three characters:
            # (key name, space, value)
            if len(kv) > 2
        )
        for k,v in pairs:
            if k not in extra_columns:
                extra_columns[k] = [""] * n
                column_order.append(k)
            extra_columns[k][i] = v

    # make a copy of the DataFrame without attribute strings.
    df = df.drop("attribute", axis=1)

    logging.info("Adding attribute columns: %s", column_order)
    for k in column_order:
        assert k not in df, "Column '%s' appears in GTF twice" % k
        df[k] = extra_columns[k]
        # delete from dictionary since these objects are big
        # and we might be runniung low on memory
        del extra_columns[k]
        # remove quotes around values
        df[k] = df[k].str.replace("\"", "")
    return df

# In addition to the required 8 columns and names and IDs genes & transcripts,
# there might also be annotations like 'transcript_biotype' but these aren't
# available for all species/releases.

REQUIRED_ATTRIBUTE_COLUMNS = [
    'gene_name',
    'gene_id',
    'gene_biotype',
    'transcript_name',
    'transcript_id',
]

def _dataframe_from_groups(groups, feature, extra_column_names=[]):
    """
    Helper function used to construct a missing feature such as 'transcript'
    or 'gene'. For example, a sufficiently old release might only have
    'exon' entries, but these were tagged with which transcript_id and gene_id
    they're associated with. Grouping by transcript_id lets you reconstruct
    the feature='transcript' entries which are normally present in later
    releases.
    """
    start = groups.start.min()
    end = groups.end.max()
    strand = groups.strand.first()
    seqname = groups.seqname.first()
    gene_name = groups.gene_name.first()
    gene_biotype = groups.gene_biotype.first()
    def pick_protein_id(candidates):
        for c in candidates:
            if c is not None and len(c) > 0:
                return c
        return None
    protein_id = groups.protein_id.apply(pick_protein_id)
    columns = [seqname, gene_biotype, start, end, strand, gene_name, protein_id]
    for column_name in extra_column_names:
        column = groups[column_name].first()
        columns.append(column)
    df = pd.concat(columns, axis=1).reset_index()

    # score seems to be always this value, not sure why it's in the GTFs
    df['score'] = '.'

    # frame values only make sense for CDS entries, but need this column
    # so we concatenate these rows with the rest of the Ensembl entries
    df['frame'] = '.'

    df['feature'] = feature
    return df

def reconstruct_gene_rows(df):
    gene_id_groups = df.groupby(['gene_id'])
    genes_df = _dataframe_from_groups(gene_id_groups, feature='gene')
    return pd.concat([df, genes_df], ignore_index=True)

def reconstruct_transcript_rows(df):
    transcript_id_groups = df.groupby(['transcript_id'])
    transcripts_df = _dataframe_from_groups(
        transcript_id_groups,
        feature='transcript',
        extra_column_names=['transcript_name']
    )
    return pd.concat([df, transcripts_df], ignore_index=True)

def reconstruct_exon_id_column(df, inplace=True):
    """
    Construct missing exon_id column for older GTFs by concatenating
    transcript_id and exon_number

        e.g. ENST00000400674.exon2

    While not useful for joining against external data, these exon_ids are
    needed for methods like ensembl_release.exon_ids_at_location.

    Parameters
    ----------

    df : DataFrame
        Must have columns 'transcript_id' and 'exon_number'

    """

    assert 'exon_id' not in df
    assert 'transcript_id' in df
    assert 'exon_number' in df

    if not inplace:
        df = df.copy()

    # missing values in dataframes indicated by NaN
    df['exon_id'] = np.nan

    # assign exon IDs to any entry with an exon number
    # this includes features = {'exon', 'CDS', 'start_codon', 'stop_codon'}
    exon_id_mask = ~df.exon_number.isnull()

    transcript_ids = df['transcript_id'][exon_id_mask]

    # convert floating value to ints (to get rid of decimal) and then str
    exon_numbers = df['exon_number'][exon_id_mask].astype(int).astype(str)

    df['exon_id'][exon_id_mask] = transcript_ids + '.exon' + exon_numbers
    return df

def load_gtf_as_dataframe(filename):
    logging.info("Reading GTF %s into DataFrame", filename)
    df = _read_gtf(filename)
    logging.info(
        "Extracting attributes for %d entries in GTF DataFrame", len(df))
    df = _extend_with_attributes(df)

    # due to the annoying ambiguity of the second GTF column,
    # figure out if an older GTF's biotype is actually the gene_biotype
    # or transcript_biotype

    if 'biotype' in df.columns:
        assert 'transcript_biotype' not in df.columns, \
            "Inferred 2nd column as biotype but also found transcript_biotype"

        # Initially we could only figure out if the 2nd column was either
        # the source of the annotation or some kind of biotype (either
        # a gene_biotype or transcript_biotype).
        # Now we disambiguate between the two biotypes by checking if
        # gene_biotype is already present in another column. If it is,
        # the 2nd column is the transcript_biotype (otherwise, it's the
        # gene_biotype)
        if 'gene_biotype' in df.columns:
            rename_to = 'transcript_biotype'
        else:
            rename_to = 'gene_biotype'
        df.rename(columns={'biotype' : rename_to}, inplace=True)

    for column_name in REQUIRED_ATTRIBUTE_COLUMNS:
        assert column_name in df.columns, \
            "Missing required column '%s', available: %s" % (
                column_name,
                list(sorted(df.columns)))

    # older Ensembl releases only had features:
    #   - exon
    #   - CDS
    #   - start_codon
    #   - stop_codon
    #
    # Might have to manually reconstruct gene & transcript entries
    # by grouping the gene_id and transcript_id columns of existing features

    distinct_features = df.feature.unique()

    if 'gene' not in distinct_features:
        logging.info("Creating entries for feature='gene'")
        df = reconstruct_gene_rows(df)

    if 'transcript' not in distinct_features:
        logging.info("Creating entries for feature='transcript'")
        df = reconstruct_transcript_rows(df)

    if 'exon_id' not in df:
        logging.info("Creating 'exon_id' column")
        df = reconstruct_exon_id_column(df)

    return df

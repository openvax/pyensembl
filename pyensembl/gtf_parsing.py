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

from __future__ import print_function, division, absolute_import
import logging
from os.path import exists
import resource
import gzip

import numpy as np
import pandas as pd
from six.moves import intern

from .locus import normalize_chromosome


"""
Parse an Ensembl GTF file into a dataframe, expanding all of its attributes
into their own columns,


Columns of a GTF file:

    seqname   - name of the chromosome or scaffold; chromosome names
                without a 'chr' in Ensembl (but sometimes with a 'chr'
                elsewhere)
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


def memory_usage():
    """
    Returns number of megabytes of memory currently being used by this Python
    process
    """
    resources = resource.getrusage(resource.RUSAGE_SELF)
    resident_bytes = resources.ru_maxrss
    return resident_bytes / (1024 * 1024)


def _read_gtf(filename, chunksize=10**5):
    """
    Download GTF annotation data for given release (if not already present),
    parse it into a DataFrame, leave the attributes in a single string column
    """
    assert exists(filename)
    assert filename.endswith(".gtf") or filename.endswith(".gtf.gz"), \
        "Must be a GTF file (%s)" % filename
    # compression = "gzip" if filename.endswith(".gz") else None

    logging.info("Memory usage before GTF parsing: %0.4f MB" % memory_usage())

    if filename.endswith("gz") or filename.endswith("gzip"):
        def open_gtf(filename):
            return gzip.open(filename, mode="rt")
    else:
        def open_gtf(filename):
            return open(filename, buffering=chunksize)

    expected_number_fields = len(GTF_COLS)

    seqname_values = []
    second_column_values = []
    feature_values = []
    start_values = []
    end_values = []
    score_values = []
    strand_values = []
    frame_values = []
    attribute_values = []

    for i, line in enumerate(open_gtf(filename)):

        if line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == expected_number_fields, \
            "Wrong number of fields %d in %s" % (len(fields, filename))
        seq, second, feature, start, end, score, strand, frame, attr = fields
        seqname_values.append(normalize_chromosome(seq))
        second_column_values.append(intern(str(second)))
        feature_values.append(intern(str(feature)))
        start_values.append(int(start))
        end_values.append(int(end))
        if score == ".":
            score_values.append(np.nan)
        else:
            score_values.append(np.float(score))
        strand_values.append(intern(str(strand)))
        frame_values.append(intern(str(frame)))
        attribute_values.append(attr)

    logging.info("Memory usage after GTF parsing: %0.4f MB" % memory_usage())

    df = pd.DataFrame({})
    df["seqname"] = np.array(seqname_values, dtype="S1")
    df["second_column"] = second_column_values
    df["feature"] = feature_values
    df["start"] = np.array(start_values, dtype="int64")
    df["end"] = np.array(end_values, dtype="int64")
    df["score"] = np.array(score_values, dtype="float32")
    df["strand"] = np.array(strand_values, dtype="S1")
    df["frame"] = np.array(frame_values, dtype="S1")
    df["attribute"] = attribute_values
    logging.info("Memory usage after DataFrame construction: %0.4f MB" % (
        memory_usage(),))

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
    df[column_name] = df["second_column"]
    del df["second_column"]
    return df

def _extend_with_attributes(df):

    logging.info(
        "Memory usage before expanding GTF attributes: %0.4f MB" % (
            memory_usage(),))

    n = len(df)

    attribute_strings = df["attribute"]
    # remove attribute strings from dataframe we're going to return
    del df["attribute"]

    extra_columns = {}
    column_order = []

    for i, attribute_string in enumerate(attribute_strings):

        # Catch mistaken semicolons by replacing "xyz;" with "xyz"
        # Required to do this since the GTF for Ensembl 78 has
        # gene_name = "PRAMEF6;"
        # transcript_name = "PRAMEF6;-201"
        if ';\"' in attribute_string:
            attribute_string = attribute_string.replace(';\"', '\"')
        if ";-" in attribute_string:
            attribute_string = attribute_string.replace(";-", "-")

        # Split the semi-colon separated attributes in the last column of a GTF
        # into a list of (key, value) pairs.
        for kv in attribute_string.split(";"):
            # need at least 3 chars for minimal entry like 'k v'
            if len(kv) < 3 or " " not in kv:
                continue
            # We're slicing the first two elements out of split() because
            # Ensembl release 79 added values like:
            #   transcript_support_level "1 (assigned to previous version 5)";
            # ...which gets mangled by splitting on spaces.
            #
            # TODO: implement a proper parser!
            column_name, value = kv.strip().split(" ", 2)[:2]

            # 1) interning keys such as "gene_name" since they reocccur
            # millions of times
            # 2) remove quotes around values
            if '\"' in value:
                value = value.replace('\"', "")
            column_name = intern(str(column_name))
            value = intern(str(value))
            column = extra_columns.get(column_name)
            if column is None:
                column = [""] * n
                extra_columns[column_name] = column
                column_order.append(column_name)
            column[i] = value
    print("Extracted GTF attributes: %s" % column_order)
    for column_name in column_order:
        column = extra_columns[column_name]
        # special case for single character strings, since these are
        # commonly occurring values
        df[column_name] = column
        del extra_columns[column_name]
    logging.info(
        "Memory usage after expanding GTF attributes: %0.4f MB" % (
            memory_usage(),))
    print("DataFrame dtypes:")
    print(df.dtypes)
    return df

# In addition to the required 8 columns and IDs of genes & transcripts,
# there might also be annotations like 'transcript_biotype' but these aren't
# available for all species/releases.

# TODO: gene and transcript names, as well as biotypes, should have
# the option to be required.
REQUIRED_ATTRIBUTE_COLUMNS = [
    'gene_id',
    'transcript_id',
]

def _dataframe_from_groups(groups, feature):
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

    # Include these columns when they're available. We also want to
    # include gene_id if we're grouping by transcript_id and vice versa.
    conditional_column_names = [
        "gene_name", "gene_biotype", "transcript_name",
        "transcript_biotype", "gene_id", "transcript_id"
    ]

    def pick_protein_id(candidates):
        for c in candidates:
            if c is not None and len(c) > 0:
                return c
        return None

    columns = [seqname, start, end, strand]

    if "protein_id" in groups.first().columns:
        protein_id = groups.protein_id.apply(pick_protein_id)
        columns.append(protein_id)

    for conditional_column_name in conditional_column_names:
        if conditional_column_name in groups.first().columns:
            column = groups[conditional_column_name].first()
            columns.append(column)

    df = pd.concat(columns, axis=1).reset_index()

    # score seems to be always '.' in Ensembl GTF files value,
    # not sure why it's even given
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
        feature='transcript'
    )
    return pd.concat([df, transcripts_df], ignore_index=True)

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
        df.rename(columns={'biotype': rename_to}, inplace=True)

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
    # (And this also applies to other non-Ensembl GTF files.)
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
    return df

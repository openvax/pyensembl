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
from os.path import split, abspath, join, exists, splitext
import pandas as pd

from typechecks import require_string
from gtfparse import read_gtf_as_dataframe, create_missing_features

from .locus import normalize_chromosome, normalize_strand
from .memory_cache import MemoryCache

class GTF(object):
    """
    Parse a GTF gene annotation file from a given local path.
    Represent its contents as a Pandas DataFrame (optionally filtered
    by locus, column, contig, &c).
    """
    def __init__(self, gtf_path, cache_directory_path=None):

        if not (gtf_path.endswith(".gtf") or gtf_path.endswith(".gtf.gz")):
            raise ValueError("Wrong extension for GTF file: %s" % (gtf_path,))
        self.gtf_path = abspath(gtf_path)

        if not exists(self.gtf_path):
            raise ValueError("GTF file %s does not exist" % (self.gtf_path,))

        self.gtf_directory_path, self.gtf_filename = split(self.gtf_path)
        self.gtf_base_filename = splitext(self.gtf_filename)[0]

        # if cache directory isn't given then put cached files
        # alongside the GTF
        if cache_directory_path:
            self.cache_directory_path = cache_directory_path
        else:
            self.cache_directory_path = self.gtf_directory_path

        self.memory_cache = MemoryCache()

        # lazily load DataFrame of all GTF entries if necessary
        # for database construction
        self._dataframes = {}

    def __eq__(self, other):
        return (
            isinstance(other, GTF) and
            other.gtf_path == self.gtf_path)

    def __hash__(self):
        return hash(self.gtf_path)

    def clear_cache(self):
        # clear any dataframes we constructed
        self._dataframes.clear()

        # clear cached dataframes loaded from CSV
        self.memory_cache.clear_cached_objects()

    def data_subset_path(
            self,
            contig=None,
            feature=None,
            column=None,
            strand=None,
            distinct=False,
            extension=".csv"):
        """
        Path to cached file for storing materialized views of the genomic data.
        Typically this is a CSV file, the filename reflects which filters have
        been applied to the entries of the database.

        Parameters:

        contig : str, optional
            Path for subset of data restricted to given contig

        feature : str, optional
            Path for subset of data restrict to given feature

        column : str, optional
            Restrict to single column

        strand : str, optional
            Positive ("+") or negative ("-") DNA strand. Default = either.

        distinct : bool, optional
            Only keep unique values (default=False)
        """
        csv_filename = self.gtf_base_filename + ".expanded"
        if contig:
            contig = normalize_chromosome(contig)
            csv_filename += ".contig.%s" % (contig,)
        if feature:
            csv_filename += ".feature.%s" % (feature,)
        if column:
            csv_filename += ".column.%s" % (column,)
        if strand:
            if strand == "+":
                strand_string = "positive"
            elif strand == "-":
                strand_string = "negative"
            else:
                raise ValueError("Invalid strand value: %s" % strand)
            csv_filename += ".strand.%s" % strand_string
        if distinct:
            csv_filename += ".distinct"
        csv_filename += extension
        return join(self.cache_directory_path, csv_filename)

    def _load_full_dataframe_cached(self):
        """
        Loads full dataframe from cached CSV or constructs it from GTF
        """
        return self.memory_cache.cached_dataframe(
            csv_path=self.data_subset_path(),
            compute_fn=self._load_full_dataframe_from_gtf)

    def _load_full_dataframe_from_gtf(self):
        """
        Parse this genome source's GTF file and load it as a Pandas DataFrame
        """
        df = read_gtf_as_dataframe(
            self.gtf_path,
            column_converters={
                "seqname": normalize_chromosome,
                "strand": normalize_strand,
            },
            infer_biotype_column=True)

        features = set(df["feature"])
        column_names = set(df.keys())

        # older Ensembl releases don't have "gene" or "transcript"
        # features, so fill in those rows if they're missing
        if "gene" not in features:
            # if we have to reconstruct gene feature rows then
            # fill in values for 'gene_name' and 'gene_biotype'
            # but only if they're actually present in the GTF
            df = create_missing_features(
                dataframe=df,
                unique_keys={"gene": "gene_id"},
                extra_columns={
                    "gene": {
                        "gene_name",
                        "gene_biotype"
                    }.intersection(column_names),
                },
                missing_value="")

        if "transcript" not in features:
            df = create_missing_features(
                dataframe=df,
                unique_keys={"transcript": "transcript_id"},
                extra_columns={
                    "transcript": {
                        "gene_id",
                        "gene_name",
                        "gene_biotype",
                        "transcript_name",
                        "transcript_biotype",
                        "protein_id",
                    }.intersection(column_names)
                },
                missing_value="")
        return df

    def dataframe(
            self,
            contig=None,
            feature=None,
            strand=None,
            save_to_disk=False):
        """
        Load genome entries as a DataFrame, optionally restricted to
        particular contig or feature type.
        """
        if contig:
            contig = normalize_chromosome(contig)

        if strand:
            strand = normalize_strand(strand)

        if feature is not None:
            require_string(feature, "feature")

        key = (contig, feature, strand)

        if key not in self._dataframes:
            def _construct_df():
                full_df = self._load_full_dataframe_cached()

                assert len(full_df) > 0, \
                    "Dataframe representation of genomic database empty!"

                # rename since we're going to be filtering the entries but
                # may still want to access the full dataset
                df = full_df
                if contig:
                    df = df[df["seqname"] == contig]
                    if len(df) == 0:
                        raise ValueError("Contig not found: %s" % (contig,))

                if feature:
                    df = df[df["feature"] == feature]
                    if len(df) == 0:
                        # check to make sure feature was somewhere in
                        # the full dataset before returning an empty dataframe
                        features = full_df["feature"].unique()
                        if feature not in features:
                            raise ValueError(
                                "Feature not found: %s" % (feature,))
                if strand:
                    df = df[df["strand"] == strand]

                return df
            if save_to_disk:
                csv_path = self.data_subset_path(
                    contig=contig,
                    feature=feature,
                    strand=strand,
                    distinct=False)
                df = self.memory_cache.cached_dataframe(
                    csv_path=csv_path,
                    compute_fn=_construct_df)
            else:
                df = _construct_df()
            self._dataframes[key] = df

        return self._dataframes[key]

    def dataframe_column_at_locus(
            self,
            column_name,
            contig,
            start,
            end=None,
            offset=None,
            strand=None):
        """
        Subset of entries which overlap an inclusive range of loci
        """
        if end is None and offset is None:
            end = start
        elif offset is None:
            end = start + offset - 1

        df_contig = self.dataframe(contig=contig, strand=strand)

        assert column_name in df_contig, \
            "Unknown genome property: %s" % column_name

        return self._slice_column(
            df_contig[column_name],
            df_contig.start,
            df_contig.end,
            start,
            end)

    @staticmethod
    def _slice_column(column_name_series, start_series, end_series,
                      start, end):
        """
        Given a Series of data, and two Series' representing start
        and end positions for the end, keep all data for which
        the start and end positions overlap the start and end
        parameters (inclusively).

        For example:
        column_name_series: ENSG00000223972, ENSG00000000003, etc.
        start_series: 11869, 100635558, etc.
        end_series: 14409, 100635252, etc.
        start = 100000000
        end = 110000000

        Returns only ENSG00000000003, since it overlaps
        [100000000, 110000000].
        """
        df = pd.DataFrame({
            'name': column_name_series,
            'start': start_series,
            'end': end_series
        })
        df = GTF._slice(df, 'start', 'end', start, end)
        return df.name

    @staticmethod
    def _slice(df, df_start_col_name, df_end_col_name, start, end):
        """
        Given a DataFrame, along with the names of columns in the
        DataFrame representing start and end positions, keep all
        data for which the start and end positions overlap the start
        and end parameters (inclusively).
        """
        # No overlap because the whole thing is before start
        mask_left = df[df_end_col_name] < start
        # No overlap because the whole thing is after end
        mask_right = df[df_start_col_name] > end
        df = df[~mask_left & ~mask_right]
        return df

    def dataframe_at_locus(
            self,
            contig,
            start,
            end=None,
            offset=None,
            strand=None):
        """
        Subset of entries which overlap an inclusive range of
        chromosomal positions
        """
        if end is None and offset is None:
            end = start
        elif offset is None:
            end = start + offset - 1

        df_contig = self.dataframe(contig=contig, strand=strand)

        # find genes whose start/end boundaries overlap with the position
        return GTF._slice(df_contig, df_contig.start.name,
                          df_contig.end.name, start, end)

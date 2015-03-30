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
from os.path import split, join
import pandas as pd

import datacache
from typechecks import require_string

from .gtf_parsing import load_gtf_as_dataframe
from .common import CACHE_SUBDIR
from .locus import normalize_chromosome, normalize_strand
from .compute_cache import cached_dataframe, clear_cached_objects
from .url_templates import ENSEMBL_FTP_SERVER, gtf_url


class GTF(object):
    """
    Download and parse a GTF gene annotation file from a given URL.
    Represent its contents as a Pandas DataFrame (optionally filtered
    by locus, column, contig, &c).
    """
    def __init__(
            self,
            release,
            species="homo_sapiens",
            server=ENSEMBL_FTP_SERVER,
            auto_download=False,
            decompress=True):
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.release = release
        self.decompress = decompress
        self.url = gtf_url(
            ensembl_release=release,
            species=species,
            server=server)

        self.remote_filename = split(self.url)[1]

        assert self.remote_filename.endswith(".gtf.gz"), \
            "Expected remote GTF file %s to end with '.gtf.gz'" % (
                self.remote_filename,)

        self.auto_download = auto_download

        # lazily load DataFrame of all GTF entries if necessary
        # for database construction
        self._dataframes = {}

    def __eq__(self, other):
        return (
            isinstance(other, GTF) and
            other.url == self.url)

    def __hash__(self):
        return hash(self.url)

    def clear_cache(self):
        # clear any dataframes we constructed
        self._dataframes.clear()

        # clear cached dataframes loaded from CSV
        clear_cached_objects()

    def base_filename(self):
        """
        Trim extensions such as ".gtf" or ".gtf.gz",
        leaving only the base filename which should be used
        to construct other derived filenames for cached data.
        """
        assert ".gtf" in self.remote_filename, \
            "GTF filename must contain .gtf extension: %s" % (
                self.remote_filename,)

        return self.remote_filename.split(".gtf")[0]

    def local_filename(self):
        """Filename used for local copy of GTF"""
        if self.decompress:
            return self.base_filename() + ".gtf"
        else:
            return self.remote_filename

    def local_copy_exists(self):
        """Has a local copy of the GTF file been downloaded?"""
        return self.cache.exists(
            self.url,
            self.local_filename(),
            decompress=self.decompress)

    def local_gtf_path(self):
        """
        Returns local path to GTF file for given release of Ensembl,
        download from the Ensembl FTP server if not already cached and
        auto download is enabled.
        """
        if self.local_copy_exists() or self.auto_download:
            # Does a download if the cache is empty
            return self.cache.fetch(
                self.url,
                self.local_filename(),
                decompress=self.decompress)
        raise ValueError('Ensembl annotation data is not currently '
                         'installed for release %s. Run '
                         '"pyensembl install --release %s" or call '
                         '"EnsemblRelease(%s).install()"' %
                         ((self.release,) * 3))

    def local_dir(self):
        return split(self.local_gtf_path())[0]

    def local_data_file_path(
            self,
            contig=None,
            feature=None,
            column=None,
            strand=None,
            distinct=False,
            extension=".csv"):
        """
        Path to local file for storing materialized views of the Ensembl data.
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
        base = self.local_filename()
        dirpath = self.local_dir()
        csv_filename = base + ".expanded"
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
        return join(dirpath, csv_filename)

    def _load_full_dataframe(self):
        """
        Loads full dataframe from cached CSV or constructs it from GTF
        """
        csv_path = self.local_data_file_path()
        return cached_dataframe(csv_path, self._load_full_dataframe_from_gtf)

    def _load_full_dataframe_from_gtf(self):
        """
        Parse this release's GTF file and load it as a Pandas DataFrame
        """
        path = self.local_gtf_path()
        return load_gtf_as_dataframe(path)

    def dataframe(self, contig=None, feature=None, strand=None):
        """
        Load Ensembl entries as a DataFrame, optionally restricted to
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
            csv_path = self.local_data_file_path(
                contig=contig,
                feature=feature,
                strand=strand,
                distinct=False)

            def local_loader_fn():
                # pylint: disable=no-member
                # pylint has trouble with df.seqname and similar
                # statements in this function.

                full_df = self._load_full_dataframe()
                assert len(full_df) > 0, \
                    "Dataframe representation of Ensembl database empty!"

                # rename since we're going to be filtering the entries but
                # may still want to access the full dataset
                df = full_df
                if contig:
                    df = df[df.seqname == contig]
                    if len(df) == 0:
                        raise ValueError("Contig not found: %s" % (contig,))

                if feature:
                    df = df[df.feature == feature]
                    if len(df) == 0:
                        # check to make sure feature was somewhere in
                        # the full dataset before returning an empty dataframe
                        features = full_df.feature.unique()
                        if feature not in features:
                            raise ValueError(
                                "Feature not found: %s" % (feature,))
                if strand:
                    df = df[df.strand == strand]

                return df

            self._dataframes[key] = cached_dataframe(csv_path, local_loader_fn)

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
            "Unknown Ensembl property: %s" % column_name

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

    def download(self, force=False):
        """
        Download the GTF file if one does not exist. If `force` is
        True, overwrites any existing file.
        """
        if not self.local_copy_exists() or force:
            self.cache.fetch(
                self.url,
                self.local_filename(),
                decompress=self.decompress,
                force=force)

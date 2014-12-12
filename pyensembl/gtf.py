from os.path import split, join
from types import NoneType

from gtf_parsing import load_gtf_as_dataframe
from common import CACHE_SUBDIR
from locus import normalize_chromosome, normalize_strand
import memory_cache
from url_templates import ENSEMBL_FTP_SERVER, gtf_url_parts

import datacache

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
            server=ENSEMBL_FTP_SERVER):

        self.cache = datacache.Cache(CACHE_SUBDIR)

        self.species = species
        self.release = release
        self.server = server

        gtf_url_dir, gtf_filename = gtf_url_parts(
            ensembl_release=release,
            species=self.species,
            server=server)
        self.filename = gtf_filename
        self.url = join(gtf_url_dir, gtf_filename)

        # lazily load DataFrame of all GTF entries if necessary
        # for database construction
        self._dataframes = {}

    def clear_cache():
        # clear any dataframes we constructed
        self._dataframes.clear()

        # clear cached dataframes loaded from CSV
        memory_cache.clear_cached_objects()

    def base_filename(self):
        """
        Trim extensions such as ".gtf" or ".gtf.gz",
        leaving only the base filename which should be used
        to construct other derived filenames for cached data.
        """
        assert ".gtf" in self.filename, \
            "GTF filename must contain .gtf extension: %s" % self.filename
        return self.filename.split(".gtf")[0]

    def local_gtf_path(self):
        """
        Returns local path to GTF file for given release of Ensembl,
        download from the Ensembl FTP server if not already cached.
        """
        return self.cache.fetch(self.url, self.filename, decompress=False)

    def local_dir(self):
        return split(self.local_gtf_path())[0]

    def local_csv_path(self, contig=None, feature=None, column=None):
        """
        Path to CSV which the annotation data with expanded columns
        for optional attributes.

        Parameters:

        contig : str, optional
            Path for subset of data restricted to given contig

        feature : str, optional
            Path for subset of data restrict to given feature

        column : str, optional
            Restrict to single column
        """
        base = self.base_filename()
        dirpath = self.local_dir()
        csv_filename = base + ".expanded"
        if contig:
            contig = normalize_chromosome(contig)
            csv_filename += ".contig.%s" % (contig,)
        if feature:
            csv_filename += ".feature.%s" % (feature,)
        if column:
            csv_filename += ".column.%s" % (column,)
        csv_filename += ".csv"
        return join(dirpath, csv_filename)

    def _load_full_dataframe(self):
        """
        Loads full dataframe from cached CSV or constructs it from GTF
        """
        csv_path = self.local_csv_path()
        return memory_cache.load_csv(
            csv_path, self._load_full_dataframe_from_gtf)


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

        if not isinstance(feature, (NoneType, str, unicode)):
            raise TypeError(
                    "Expected feature to be string, got %s : %s" % (
                        feature, type(feature)))

        key = (contig, feature, strand)

        if key not in self._dataframes:
            csv_path = self.local_csv_path(contig=contig, feature=feature)

            def local_loader_fn():
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

            self._dataframes[key] = memory_cache.load_csv(
                csv_path, local_loader_fn)

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
            end = position
        elif offset is None:
            end = position + offset - 1

        df_contig = self.dataframe(contig=contig, strand=strand)

        assert column_name in df_contig, \
            "Unknown Ensembl property: %s" % column_name

        return self._slice_column(
            df_contig[column_name],
            df_contig.start,
            df_contig.end,
            start,
            end)

    def dataframe_at_locus(self, contig, position, end=None, offset=None):
        """
        Subset of entries which overlap an inclusive range of
        chromosomal positions
        """
        if end is None and offset is None:
            end = position
        elif offset is None:
            end = position + offset - 1

        df_contig = self.dataframe(contig=contig)

        # find genes whose start/end boundaries overlap with the position
        overlap_start = df_contig.start <= end
        overlap_end = df_contig.end >= position
        overlap = overlap_start & overlap_end
        return df_contig[overlap]

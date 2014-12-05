"""
Contains the EnsemblRelease class, with its millions of accessors and wrappers
around the Ensembl annotation database.
"""

from glob import glob
import logging
from os.path import join, exists, split
from os import remove
import sqlite3
from types import NoneType

from exon import Exon
from gene import Gene
from gtf import load_gtf_as_dataframe
from locus import normalize_chromosome, normalize_strand
import memory_cache
from release_info import which_human_reference, check_release_version
from transcript import Transcript


import datacache
import numpy as np
import pandas as pd


# directory which contains GTF files, missing the release number
URL_DIR_TEMPLATE = 'ftp://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/'
FILENAME_TEMPLATE = "Homo_sapiens.%s.%d.gtf.gz"
CACHE_SUBDIR = "ensembl"

class EnsemblRelease(object):

    def __init__(self, release, lazy_loading=True):
        self.release = check_release_version(release)
        self.gtf_url_dir = URL_DIR_TEMPLATE % self.release
        self.reference_name =  which_human_reference(self.release)
        self.gtf_filename = FILENAME_TEMPLATE  % (
            self.reference_name, self.release
        )
        self.gtf_url = join(self.gtf_url_dir, self.gtf_filename)
        self._local_gtf_path = None

        # lazily load DataFrame of all GTF entries if necessary
        # for database construction
        self._df_cache = {}

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

        if lazy_loading:
            # lazily construct sqlite3 database
            self._db = None
        else:
            self._db = self._connect_or_create_database()

    def __str__(self):
        return "EnsemblRelease(release=%d, gtf_url='%s')" % (
            self.release, self.gtf_url)

    def __repr__(self):
        return str(self)

    def base_gtf_filename(self):
        """
        Trim extensions such as ".gtf" or ".gtf.gz",
        leaving only the base filename which should be used
        to construct other derived filenames for cached data.
        """
        assert ".gtf" in self.gtf_filename, \
            "GTF filename must contain .gtf extension: %s" % self.gtf_filename
        parts = self.gtf_filename.split(".gtf")
        return parts[0]

    def _delete_cached_files(self):
        base = self.base_gtf_filename()
        dirpath = self.local_gtf_dir()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)

    def clear_cached_data(self):
        self._df_cache = {}
        self._delete_cached_files()
        # TODO: combine caching of parsed CSV files and
        # caching of partial reads from the database
        memory_cache.clear_cached_objects()


    def local_gtf_path(self):
        """
        Returns local path to GTF file for given release of Ensembl,
        download from the Ensembl FTP server if not already cached.
        """
        if self._local_gtf_path is None:
            self._local_gtf_path = datacache.fetch_file(
                self.gtf_url,
                filename=self.gtf_filename,
                decompress=False,
                subdir=CACHE_SUBDIR)
        assert self._local_gtf_path is not None
        return self._local_gtf_path

    def local_gtf_dir(self):
        return split(self.local_gtf_path())[0]

    def local_csv_path(self, contig=None, feature=None):
        """
        Path to CSV which the annotation data with expanded columns
        for optional attributes.

        Parameters:

        contig : str, optional
            Path for subset of data restricted to given contig

        feature : str, optional
            Path for subset of data restrict to given feature
        """
        base = self.base_gtf_filename()
        dirpath = self.local_gtf_dir()
        csv_filename = base + ".expanded"
        if contig:
            contig = normalize_chromosome(contig)
            csv_filename += ".contig.%s" % (contig,)
        if feature:
            csv_filename += ".feature.%s" % (feature,)
        csv_filename += ".csv"
        return join(dirpath, csv_filename)

    def local_db_filename(self):
        base = self.base_gtf_filename()
        return base + ".db"

    def local_db_path(self):
        dirpath = self.local_gtf_dir()
        filename = self.local_db_filename()
        return join(dirpath, filename)

    def _load_full_dataframe_from_gtf(self):
        """
        Parse this release's GTF file and load it as a Pandas DataFrame
        """
        path = self.local_gtf_path()
        return load_gtf_as_dataframe(path)

    def _load_full_dataframe(self):
        """
        Loads full dataframe from cached CSV or constructs it from GTF
        """
        csv_path = self.local_csv_path()
        return memory_cache.load_csv(
            csv_path, self._load_full_dataframe_from_gtf)

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

        if key not in self._df_cache:
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

            self._df_cache[key] = memory_cache.load_csv(
                csv_path, local_loader_fn)

        return self._df_cache[key]

    def _database_indices(self, available_columns):
        """
        If a release is missing some column we want to index on,
        we have to drop any indices which use that column.
        """
        candidate_column_groups = [
            ['seqname', 'start', 'end'],
            ['seqname', 'start', 'end', 'strand'],
            ['seqname'],
            ['gene_name'],
            ['gene_id'],
            ['transcript_id'],
            ['transcript_name'],
            ['exon_id'],
        ]
        indices = []

        # Since queries are often restricted by feature type
        # we should include that column in combination with all
        # other indices we anticipate might improve performance
        for column_group in candidate_column_groups:
            skip = False
            for column_name in column_group:
                # some columns, such as 'exon_id',
                # are not available in all releases of Ensembl
                if column_name not in available_columns:
                    self.logger.info(
                        "Skipping database index for %s",
                        "+".join(column_group)
                    )
                    skip = True
            if skip:
                continue
            indices.append(column_group)
            indices.append(column_group + ['feature'])
        indices.append(['feature'])
        return indices

    def _create_database(self):
        print "Creating database: %s" % self.local_db_path()
        filename = self.local_db_filename()
        df = self.dataframe()

        available_columns = set(df.columns)
        indices = self._database_indices(available_columns)

        db = datacache.db_from_dataframe(
            db_filename=filename,
            table_name="ensembl",
            df=df,
            subdir=CACHE_SUBDIR,
            overwrite=False,
            indices=indices)
        return db

    def _connect_or_create_database(self):
        """
        If database already exists, open a connection.
        Otherwise, create it.
        """
        db_path = self.local_db_path()
        if exists(db_path):
            db = sqlite3.connect(db_path)
            # maybe file got created but not filled
            if datacache.db.db_table_exists(db, 'ensembl'):
                return db
        return self._create_database()

    def db(self):
        """
        Return the sqlite3 database for this Ensembl release
        (download and/or construct it necessary).
        As a side effect, stores the database connection in self._db
        """
        if self._db is None:
            self._db = self._connect_or_create_database()
        return self._db

    def _slice_column(self, values, starts, ends, start, end):
        """
        Return subset of values array which overlap with the given
        start/end positions.
        """
        overlap_start = starts <= end
        overlap_end = ends >= start
        # find genes whose start/end boundaries overlap with the position
        return values[overlap_start & overlap_end]

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

    def _db_columns(self, db):
        table_info = db.execute("PRAGMA table_info(ensembl);").fetchall()
        return [info[1] for info in table_info]

    def _db_column_exists(self, db, column_name):
        return column_name in self._db_columns(db)

    def column_at_locus(
            self,
            column_name,
            contig,
            position,
            end=None,
            strand=None):
        """
        Get the non-null values of a column from the database
        at a particular range of loci
        """

        if not isinstance(column_name, (str, unicode)):
            raise TypeError(
                "Expected column_name to be str, got %s : %s" % (
                    column_name, type(column_name)))

        contig = normalize_chromosome(contig)

        if not isinstance(position, (int, long)):
            raise TypeError(
                "Expected position to be integer, got %s : %s" % (
                    position, type(position)))

        if end is None:
            end = position

        if not isinstance(end, (int, long)):
            raise TypeError(
                "Expected end to be integer, got %s : %s" % (end, type(end)))

        db = self.db()

        if not self._db_column_exists(db, column_name):
            raise ValueError("Unknown Ensembl property: %s" % (column_name,))

        query = """
            SELECT %s
            FROM ensembl
            WHERE seqname=?
            AND start <= ?
            AND end >= ?
        """  % (column_name,)

        query_params = [contig, end, position]

        if strand:
            query += " AND strand = ?"
            query_params.append(strand)

        # self.logger.info("Running query: %s" % query)
        results = db.execute(query, query_params).fetchall()
        # each result is a tuple, so pull out its first element
        return [result[0] for result in results if result[0] is not None]


    def _property_values_at_locus(
            self, property_name, contig, position, end=None, strand=None):
        col = self.column_at_locus(
            property_name,
            contig,
            position,
            end=end,
            strand=strand)
        return list(sorted({c for c in col if c}))

    def genes_at_locus(self, contig, position, end=None, strand=None):
        gene_ids = self.gene_ids_at_locus(
            contig, position, end=end, strand=strand)
        return [self.gene_by_id(gene_id) for gene_id in gene_ids]

    def transcripts_at_locus(self, contig, position, end=None, strand=None):
        transcript_ids = self.transcript_ids_at_locus(
            contig, position, end=end, strand=strand)
        return [
            self.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    def exons_at_locus(self, contig, position, end=None, strand=None):
        exon_ids = self.exon_ids_at_locus(
            contig, position, end=end, strand=strand)
        return [
            self.exon_by_id(exon_id)
            for exon_id in exon_ids
        ]

    def gene_ids_at_locus(self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
            'gene_id', contig, position, end=end, strand=strand)

    def gene_names_at_locus(self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
             'gene_name', contig, position, end=end, strand=strand)

    def exon_ids_at_locus(self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
            'exon_id', contig, position, end=end, strand=strand)

    def transcript_ids_at_locus(self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
            'transcript_id', contig, position, end=end, strand=strand)

    def transcript_names_at_locus(
            self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
            'transcript_name', contig, position, end=end, strand=strand)

    def protein_ids_at_locus(self, contig, position, end=None, strand=None):
        return self._property_values_at_locus(
            'protein_id', contig, position, end=end, strand=strand)

    def run_sql_query(self, sql, required=False):
        """
        Given an arbitrary SQL query, run it against the Ensembl database
        and return the results.
        """
        db = self.db()
        cursor = db.execute(sql)
        results = cursor.fetchall()
        if required and len(results) == 0:
            raise ValueError(
                "No results found in Ensembl for query:\n%s" % (sql,))
        return results

    def query(
            self,
            select_column_names,
            filter_column,
            filter_value,
            feature,
            distinct=False,
            required=False):
        """
        Construct a SQL query and run against the ensembl sqlite3 database,
        filtered both by the feature type and a user-provided column/value.
        """
        query = """
            select %s%s
            from ensembl
            where
                %s = '%s' and
                feature='%s'
        """ % ("distinct " if distinct else "",
               ", ".join(select_column_names),
               filter_column,
               filter_value,
               feature)
        return self.run_sql_query(query, required=required)

    def _query_feature(self, column, feature, distinct=True):
        """
        Run a SQL query against the ensembl sqlite3 database, filtered
        only on the feature type.
        """
        query = "select %s%s from ensembl where feature='%s'" % (
               "distinct " if distinct else "",
               column,
               feature
        )
        return self.run_sql_query(query)

    def _query_distinct_on_contig(self, column_name, feature, contig):
        contig = normalize_chromosome(contig)
        results = self.query(
            select_column_names=[column_name],
            filter_column="seqname",
            filter_value=contig,
            feature=feature,
            distinct=True)
        return [result[0] for result in results if result is not None]

    ###################################################
    #
    #         Locations: (contig, start, stop)
    #
    ###################################################

    def _get_locations(self, property_name, property_value, feature):
        return self.query(
            select_column_names=["seqname", "start", "end", "strand"],
            filter_column=property_name,
            filter_value=property_value,
            feature=feature,
            distinct=True,
            required=True)

    def locations_of_gene_name(self, gene_name):
        """
        Given a gene name returns list of tuples with fields:
            (chromosome, start, stop, strand)
        You can get multiple results since a gene might have multiple copies
        in the genome.
        """
        return self._get_locations('gene_name', gene_name, 'gene')

    def _get_unique_location(self, property_name, property_value, feature):
        locations = self._get_locations(property_name, property_value, feature)
        if len(locations) == 0:
            raise ValueError("%s not found: %s" % (feature, gene_id,))
        elif len(locations) > 1:
            raise ValueError("%s has multiple loci: %s" % (feature, gene_id))
        return locations[0]

    def location_of_gene_id(self, gene_id):
        """
        Given a gene ID returns (chromosome, start, stop, strand)
        """
        return self._get_unique_location('gene_id', gene_id, 'gene')

    def location_of_transcript_id(self, transcript_id):
        return self._get_unique_location(
            'transcript_id',
            transcript_id,
            feature='transcript')

    def location_of_exon_id(self, exon_id):
        """
        Given an exon ID returns (chromosome, start, stop)
        """
        return self._get_unique_location('exon_id', exon_id, feature='exon')

    ###################################################
    #
    #             Gene Info Objects
    #
    ###################################################

    def gene_by_id(self, gene_id):
        return Gene(gene_id, self.db())

    def genes_by_name(self, gene_name):
        gene_ids = self.gene_ids_of_gene_name(gene_name)
        return [self.gene_by_id(gene_id) for gene_id in gene_ids]

    def gene_by_protein_id(self, protein_id):
        gene_id = self.gene_id_of_protein_id(protein_id)
        return self.gene_by_id(gene_id)

    ###################################################
    #
    #             Gene Names
    #
    ###################################################

    def _query_gene_name(self, property_name, property_value, feature_type):
        results = self.query(
            select_column_names=["gene_name"],
            filter_column=property_name,
            filter_value=property_value,
            feature=feature_type,
            distinct=True,
            required=True)
        return str(results[0][0])

    def gene_names(self):
        return self._query_feature('gene_name', 'gene')

    def gene_names_on_contig(self, contig):
        """
        Return all genes on given chromosome/contig
        """
        return self._query_distinct_on_contig(
            column_name='gene_name',
            feature='gene',
            contig=contig)

    def gene_name_of_gene_id(self, gene_id):
        return self._query_gene_name("gene_id", gene_id, 'gene')

    def gene_name_of_transcript_id(self, transcript_id):
        return self._query_gene_name(
            "transcript_id", transcript_id, 'transcript')

    def gene_name_of_transcript_name(self, transcript_name):
        return self._query_gene_name(
            "transcript_name", transcript_name, 'transcript')

    def gene_name_of_exon_id(self, transcript_id):
        return self._query_gene_name("exon_id", exon_id, 'exon')

    ###################################################
    #
    #             Gene IDs
    #
    ###################################################

    def gene_ids(self):
        return self._query_feature('gene_id', 'gene')

    def gene_ids_on_contig(self, contig):
        """
        What are all the gene IDs on a given chromosome/contig?
        """
        return self._query_distinct_on_contig(
            column_name='gene_id',
            feature='gene',
            contig=contig)

    def gene_ids_of_gene_name(self, gene_name):
        """
        What are the Ensembl gene IDs associated with a given gene name?
        (due to copy events, there might be multiple genes per name)
        """
        results = self.query(
            select_column_names=["gene_id"],
            filter_column="gene_name",
            filter_value=gene_name,
            feature="gene",
            required=True)
        results = [
            str(result_tuple[0])
            for result_tuple in results
            if result_tuple[0]
        ]

        if len(results) == 0:
            raise ValueError("Gene name not found: %s" % gene_name)

        return results


    ###################################################
    #
    #             Transcript Info Objects
    #
    ###################################################

    def transcript_by_id(self, transcript_id):
        return Transcript(transcript_id, self.db())

    def transcripts_by_name(self, transcript_name):
        transcript_ids = self.transcript_ids_of_transcript_name(transcript_name)
        return [
            self.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    def transcript_by_name(self, transcript_name):
        """
        Get the single transcript associated with a particular name,
        raise an exception if there are zero or multiple transcripts.
        """
        transcripts = self.transcripts_by_name(transcript_name)

        if len(transcripts) == 0:
            raise ValueError(
                "No transcripts found with name = %s" % transcript_name)
        elif len(transcripts) > 1:
            raise ValueError(
                "Multiple transcripts found with name = %s (IDs = %s)" %(
                    transcript_name,
                    [transcript.id for transcript in transcripts]))
        return transcripts[0]

    def transcript_by_protein_id(self, protein_id):
        transcript_id = self.transcript_id_of_protein_id(protein_id)
        return self.transcript_by_id(transcript_id)

    ###################################################
    #
    #            Transcript Names
    #
    ###################################################

    def _query_transcript_names(self, property_name, value):
        results = self.query(
            select_column_names=['transcript_name'],
            filter_column=property_name,
            filter_value=value,
            feature='transcript',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def transcript_names(self):
        return self._query_feature('transcript_name', 'transcript')

    def transcript_names_on_contig(self, contig):
        """
        What are all the transcript names on a given chromosome/contig?
        """
        return self._query_distinct_on_contig(
            column_name='transcript_name',
            feature='transcript',
            contig=contig)

    def transcript_name_of_gene_name(self, gene_name):
        return _query_transcript_names('gene_name', gene_name)[0]

    def transcript_name_of_transcript_id(self, transcript_id):
        return _query_transcript_names('transcript_id', transcript_id)[0]

    ###################################################
    #
    #            Transcript IDs
    #
    ###################################################

    def _query_transcript_ids(self, property_name, value):
        results = self.query(
            select_column_names=['transcript_id'],
            filter_column=property_name,
            filter_value=value,
            feature='transcript',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def transcript_ids_of_gene_id(self, gene_id):
        return self._query_transcript_ids('gene_id', gene_id)

    def transcript_ids_of_gene_name(self, gene_name):
        return self._query_transcript_ids('gene_name', gene_name)

    def transcript_ids_of_transcript_name(self, transcript_name):
        return self._query_transcript_ids('transcript_name', transcript_name)

    def transcript_ids_of_exon_id(self, exon_id):
        return self._query_transcript_ids('exon_id', exon_id)

    ###################################################
    #
    #             Exon Info Objects
    #
    ###################################################

    def exon_by_id(self, exon_id):
        return Exon(exon_id, self.db())

    def exon_by_transcript_and_number(self, transcript_id, exon_number):
        transcript = self.transcript_by_id(transcript_id)
        if len(transcript.exons) > exon_number:
            raise ValueError(
                "Invalid exon number for transcript %s" % transcript_id)

        # exon numbers in Ensembl are 1-based, need to subtract 1 to get
        # a list index
        return transcript.exons[exon_number - 1]

    ###################################################
    #
    #                Exon IDs
    #
    ###################################################

    def _query_exon_ids(self, property_name, value):
        results = self.query(
            select_column_names=['exon_id'],
            filter_column=property_name,
            filter_value=value,
            feature='exon',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def exon_ids(self):
        return self._query_feature('exon_id', 'exon')

    def exon_ids_of_gene_id(self, gene_id):
        return self._query_exon_ids('gene_id', gene_id)

    def exon_ids_of_gene_name(self, gene_name):
        return self._query_exon_ids('gene_name', gene_name)

    def exon_ids_of_transcript_name(self, transcript_name):
        return self._query_exon_ids('transcript_name', transcript_name)

    def exon_ids_of_transcript_id(self, transcript_id):
        return self._query_exon_ids('transcript_id', transcript_id)


    ###################################################
    #
    #                Start Codons
    #
    ###################################################

    def _query_start_codon_location(self, property_name, value):
        return self.query(
            select_column_names=['seqname', 'start', 'end', 'strand'],
            filter_column=property_name,
            filter_value=value,
            feature='start_codon',
            distinct=True,
            required=True)

    def start_codon_of_transcript_id(self, transcript_id):
        return self._query_start_codon_location('transcript_id', transcript_id)

    def start_codon_of_transcript_name(self, transcript_name):
        return self._query_start_codon_location(
            'transcript_name', transcript_name)

    ###################################################
    #
    #                Stop Codons
    #
    ###################################################

    def _query_stop_codon_location(self, property_name, value):
        return self.query(
            select_column_names=['seqname', 'start', 'end', 'strand'],
            filter_column=property_name,
            filter_value=value,
            feature='start_codon',
            distinct=True,
            required=True)

    def stop_codon_of_transcript_id(self, transcript_id):
        return self._query_stop_codon_location('transcript_id', transcript_id)

    def stop_codon_of_transcript_name(self, transcript_name):
        return self._query_stop_codon_location(
            'transcript_name', transcript_name)





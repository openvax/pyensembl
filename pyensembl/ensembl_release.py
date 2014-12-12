"""
Contains the EnsemblRelease class, with its millions of accessors and wrappers
around the Ensembl annotation database.
"""

from glob import glob
import logging
from os.path import join, exists, split
from os import remove
import sqlite3

from common import CACHE_SUBDIR
from exon import Exon
from gene import Gene
from gtf import GTF
from locus import normalize_chromosome, normalize_strand
from reference_transcripts import ReferenceTranscripts
from release_info import check_release_number
from transcript import Transcript
from url_templates import ENSEMBL_FTP_SERVER

import datacache
import numpy as np
import pandas as pd


class EnsemblRelease(object):

    def __init__(self, release, lazy_loading=True, server=ENSEMBL_FTP_SERVER):
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.release = check_release_number(release)
        self.species = "homo_sapiens"
        self.server = server
        self.gtf = GTF(self.release, self.species, server)
        self.reference = ReferenceTranscripts(
            self.release, self.species, server)

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

        if lazy_loading:
            # lazily construct sqlite3 database
            self._db = None
        else:
            self._db = self._connect_or_create_database()

    def __str__(self):
        return "EnsemblRelease(release=%d, gtf_url='%s', fasta_url='%s')" % (
            self.release, self.gtf.url, self.reference.url)

    def __repr__(self):
        return str(self)

    def _delete_cached_files(self):
        """
        Any files which start with the same name as our GTF file
        is assumed to be some view of the this release's data and thus
        safe to delete.
        """
        base = self.gtf.base_filename()
        dirpath = self.local_gtf_dir()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)

    def clear_cache(self):
        self.gtf.clear_cache()
        self.reference.clear_cache()
        self._delete_cached_files()

    def local_db_filename(self):
        base = self.gtf.base_filename()
        return base + ".db"

    def local_db_path(self):
        dirpath = self.gtf.local_dir()
        filename = self.local_db_filename()
        return join(dirpath, filename)

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
        df = self.gtf.dataframe()

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

    @property
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



    def _db_columns(self, db):
        table_info = db.execute("PRAGMA table_info(ensembl);").fetchall()
        return [info[1] for info in table_info]

    def _db_column_exists(self, db, column_name):
        return column_name in self._db_columns(db)

    def column_values_at_locus(
            self,
            column_name,
            feature,
            contig,
            position,
            end=None,
            strand=None,
            distinct=False,
            sorted=False):
        """
        Get the non-null values of a column from the database
        at a particular range of loci
        """

        # TODO: combine with the query method, since they overlap
        # significantly

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

        if not self._db_column_exists(self.db, column_name):
            raise ValueError("Unknown Ensembl property: %s" % (column_name,))

        if distinct:
            distinct_string = "DISTINCT "
        else:
            distinct_string = ""

        query = """
            SELECT %s%s
            FROM ensembl
            WHERE feature = ?
            AND seqname=?
            AND start <= ?
            AND end >= ?

        """  % (distinct_string, column_name,)

        query_params = [feature, contig, end, position]

        if strand:
            query += " AND strand = ?"
            query_params.append(strand)

        tuples = self.db.execute(query, query_params).fetchall()
        # each result is a tuple, so pull out its first element
        results = [t[0] for t in tuples if t[0] is not None]
        if sorted:
            results.sort()
        return results


    def distinct_column_values_at_locus(
            self,
            column,
            feature,
            contig,
            position,
            end=None,
            strand=None):
        """
        Gather all the distinct values for a property/column at some specified
        locus.

        Parameters
        ----------
        column : str
            Which property are we getting the values of.

        feature : str
            Which type of entry (e.g. transcript, exon, gene) is the property
            associated with?

        contig : str
            Chromosome or unplaced contig name

        position : int
            Chromosomal position

        end : int, optional
            End position of a range, if unspecified assume we're only looking
            at the single given position.

        strand : str, optional
            Either the positive ('+') or negative strand ('-'). If unspecified
            then check for values on either strand.
        """
        return self.column_values_at_locus(
            column,
            feature,
            contig,
            position,
            end=end,
            strand=strand,
            distinct=True,
            sorted=True)

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
        return self.distinct_column_values_at_locus(
            column='gene_id',
            feature='gene',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def gene_names_at_locus(self, contig, position, end=None, strand=None):
        return self.distinct_column_values_at_locus(
             column='gene_name',
             feature='gene',
             contig=contig,
             position=position,
             end=end,
             strand=strand)

    def exon_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.distinct_column_values_at_locus(
            column='exon_id',
            feature='exon',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.distinct_column_values_at_locus(
            column='transcript_id',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_names_at_locus(
            self, contig, position, end=None, strand=None):
        return self.distinct_column_values_at_locus(
            column='transcript_name',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def protein_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.distinct_column_values_at_locus(
            column='protein_id',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def run_sql_query(self, sql, required=False, query_params=[]):
        """
        Given an arbitrary SQL query, run it against the Ensembl database
        and return the results.
        """
        cursor = self.db.execute(sql, query_params)
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
            SELECT %s%s
            FROM ensembl
            WHERE %s = ?
            AND feature= ?
        """ % ("distinct " if distinct else "",
               ", ".join(select_column_names),
               filter_column)
        query_params = [filter_value, feature]
        return self.run_sql_query(
            query, required=required, query_params=query_params)

    def _query_feature_values(
            self,
            column,
            feature,
            distinct=True,
            contig=None):
        """
        Run a SQL query against the ensembl sqlite3 database, filtered
        only on the feature type.
        """
        query = """
            SELECT %s%s
            FROM ensembl
            WHERE feature=?
        """ % ("DISTINCT " if distinct else "", column)
        query_params = [feature]

        if contig:
            contig = normalize_chromosome(contig)
            query += " AND seqname = ?"
            query_params.append(contig)

        rows = self.run_sql_query(query, query_params = query_params)
        return [row[0] for row in rows if row is not None]

    def _all_feature_values(self, column, feature, distinct=True, contig=None):
        """
        Cached lookup of all values for a particular feature property

        Parameters
        ----------

        column : str
            Name of property (e.g. exon_id)

        feature : str
            Type of entry (e.g. exon)

        distinct : bool, optional
            Keep only unique values

        contig : str, optional
            Restrict query to particular contig

        Returns a single column Pandas DataFrame constructed from query results.
        """
        csv_path = self.local_csv_path(
            feature=feature,
            column=column,
            contig=contig)
        def run_query():
            values = self._query_feature_values(
                column=column,
                feature=feature,
                distinct=distinct,
                contig=contig)
            return pd.DataFrame({column : values})
        return memory_cache.load_csv(csv_path, expensive_action=run_query)

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
        return Gene(gene_id, self.db)

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
        return self._all_feature_values('gene_name', 'gene')

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
        return self._all_feature_values('gene_id', 'gene')

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

    def transcripts(self, contig=None):
        """
        Construct Transcript object for every transcript entry in
        the Ensembl database. Optionally restrict to a particular
        chromosome using the `contig` argument.
        """
        # DataFrame with single column 'transcript_id'
        transcript_ids_df = self.transcript_ids(contig=contig)
        db = self.db
        return [
            Transcript(transcript_id, db)
            for transcript_id in transcript_ids_df['transcript_id']
        ]

    def transcript_by_id(self, transcript_id):
        return Transcript(transcript_id, self.db)

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

    def transcript_names(self, contig=None):
        return self._all_feature_values(
            column='transcript_name',
            feature='transcript',
            contig=contig)

    def transcript_names_on_contig(self, contig):
        """
        What are all the transcript names on a given chromosome/contig?
        """
        return self._query_distinct_on_contig(
            column_name='transcript_name',
            feature='transcript',
            contig=contig)

    def transcript_names_of_gene_name(self, gene_name):
        return self._query_transcript_names('gene_name', gene_name)

    def transcript_name_of_transcript_id(self, transcript_id):
        transcript_names = self._query_transcript_names(
            'transcript_id', transcript_id)
        if len(transcript_names) == 0:
            raise ValueError(
                "No transcript names for transcript ID = %s" % transcript_id)
        assert len(transcript_names) == 1, \
            "Multiple transcript names for transcript ID = %s" % transcript_id
        return transcript_names[0]

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

    def transcript_ids(self, contig=None):
        return self._all_feature_values(
            column='transcript_id',
            feature='transcript',
            contig=contig)

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

    def exons(self, contig=None):
        """
        Create exon object for all exons in the database, optionally
        restrict to a particular chromosome using the `contig` argument.
        """
        # DataFrame with single column called 'exon_id'
        exon_ids_df = self.exon_ids(contig=contig)
        db = self.db
        return [Exon(exon_id, db) for exon_id in exon_ids_df['exon_id']]

    def exon_by_id(self, exon_id):
        return Exon(exon_id, self.db)

    def exon_by_transcript_id_and_number(self, transcript_id, exon_number):
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

    def exon_ids(self, contig=None):
        return self._all_feature_values(
            column='exon_id',
            feature='exon',
            contig=contig)

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






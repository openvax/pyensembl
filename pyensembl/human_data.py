from glob import glob
import logging
from os.path import join, exists, split
from os import remove
import sqlite3

from gtf import load_gtf_as_dataframe
from locus import normalize_chromosome
from memory_cache import load_csv, clear_cached_objects

import datacache
import numpy as np
import pandas as pd


MIN_ENSEMBL_RELEASE = 48
MAX_ENSEMBL_RELEASE = 77

def _check_release(release):
    """
    Convert a user-provided release number into
    an integer, check to make sure it's in the
    valid range of Ensembl releases
    """
    try:
        release = int(release)
    except:
       assert False, "%s is not a valid Ensembl release" % release
    assert release >= MIN_ENSEMBL_RELEASE
    assert release <= MAX_ENSEMBL_RELEASE
    return release



# mapping from Ensembl release to which reference assembly it uses
_human_references = {}

# Ensembl release 48-54 use NCBI36 as a reference
for i in xrange(48,55):
    _human_references[i] = 'NCBI36'

# Ensembl releases 55-75 use CRCh37 as a reference
for i in xrange(55,76):
    _human_references[i] = 'GRCh37'

# Ensembl releases 76 and 77 use GRCh38
for i in xrange(76,78):
    _human_references[i] = 'GRCh38'

def _which_human_reference(release):
    release = _check_release(release)
    assert release in _human_references, \
        "No reference found for release %d" % release
    return _human_references[release]


# directory which contains GTF files, missing the release number
URL_DIR_TEMPLATE = 'ftp://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/'
FILENAME_TEMPLATE = "Homo_sapiens.%s.%d.gtf.gz"
CACHE_SUBDIR = "ensembl"

class EnsemblRelease(object):

    def __init__(self, release):
        self.release = _check_release(release)
        self.gtf_url_dir = URL_DIR_TEMPLATE % self.release
        self.reference_name =  _which_human_reference(self.release)
        self.gtf_filename = FILENAME_TEMPLATE  % (
            self.reference_name, self.release
        )
        self.gtf_url = join(self.gtf_url_dir, self.gtf_filename)
        self._local_gtf_path = None

        # lazily load DataFrame of all GTF entries
        self._df = None

        # lazily construct sqlite3 database
        self._db = None

        # dictionary mapping each contig name to a dictionary of columns
        self._contig_columns = {}

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

    def _clear_contig_columns(self):
        for d in self._contig_columns:
            d.clear()
        self._contig_columns.clear()

    def _load_contig_column_from_db(self, contig_name, col_name):
        db = self.db()
        assert self._db_col_exists(db, col_name), \
            "Column %s not in Ensembl database" % col_name
        query = "select %s from ensembl where seqname='%s' order by ROWID;" % (
            col_name, contig_name)
        df = pd.read_sql(query, db)
        return df[col_name]

    def _contig_column(self, contig_name, col_name):
        if contig_name not in self._contig_columns:
            self._contig_columns[contig_name] = {}

        d = self._contig_columns[contig_name]

        # if this contig/column pair isn't already cached, then
        # extract it from the database
        if col_name in d:
            col = d[col_name]
        else:
            col = self._load_contig_column_from_db(contig_name, col_name)
            if not col.dtype.hasobject:
              col = np.array(col)
            d[col_name] = col
        return col

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

    def delete_cached_files(self):
        base = self.base_gtf_filename()
        dirpath = self.local_gtf_dir()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)
        # TODO: combine caching of parsed CSV files and
        # caching of partial reads from the database
        clear_cached_objects()
        self._clear_contig_columns()

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
        assert self._local_gtf_path
        return self._local_gtf_path

    def local_gtf_dir(self):
        return split(self.local_gtf_path())[0]

    def local_csv_path(self):
        """
        Path to CSV which the annotation data with expanded columns
        for optional attributes
        """
        base = self.base_gtf_filename()
        dirpath = self.local_gtf_dir()
        csv_filename = base + ".expanded.csv"
        return join(dirpath, csv_filename)

    def local_db_filename(self):
        base = self.base_gtf_filename()
        return base + ".db"

    def local_db_path(self):
        dirpath = self.local_gtf_dir()
        filename = self.local_db_filename()
        return join(dirpath, filename)

    def local_csv_path_for_contig(self, contig_name):
        """
        Path to CSV file containing subset of Ensembl data
        restricted to given contig_name
        """
        dirpath = self.local_gtf_dir()
        base = self.base_gtf_filename()
        contig_name = normalize_chromosome(contig_name)
        return join(dirpath, base + ".contig.%s.csv" % contig_name)


    def _load_dataframe_from_gtf(self):
        """
        Parse this release's GTF file and load it as a Pandas DataFrame
        """
        path = self.local_gtf_path()
        return load_gtf_as_dataframe(path)


    def dataframe(self):
        if self._df is None:
            csv_path = self.local_csv_path()
            self._df = load_csv(csv_path, self._load_dataframe_from_gtf)
        return self._df

    def _create_database(self):
        df = self.dataframe()
        filename = self.local_db_filename()
        print "Creating database: %s" % self.local_db_path()
        db = datacache.db_from_dataframe(
            db_filename=filename,
            table_name="ensembl",
            df=df,
            subdir=CACHE_SUBDIR,
            overwrite=False,
            indices = [
                ['seqname', 'start', 'end'],
                ['seqname', 'start', 'end', 'strand'],
                ['seqname'],
                ['gene_name'],
                ['gene_id'],
                ['transcript_id'],
                ['exon_id'],
                ['feature'],
            ])
        return db


    def db(self):
        if self._db is None:
            db_path = self.local_db_path()
            if exists(db_path):
                db = sqlite3.connect(db_path)
                # maybe file got created but not filled
                if not datacache.db.db_table_exists(db, 'ensembl'):
                    db = self._create_database()
            else:
                db = self._create_database()
            self._db = db
        return self._db

    def dataframe_for_contig(self, contig_name):
        """
        Load a subset of the Ensembl data for a specific contig
        """
        contig_name = normalize_chromosome(contig_name)
        contig_csv_path = self.local_csv_path_for_contig(contig_name)
        def create_dataframe():
            df = self.dataframe()
            mask = df.seqname == contig_name
            subset = df[mask]
            assert len(subset) > 0, "Contig not found: %s" % contig_name
            return subset
        return load_csv(contig_csv_path, create_dataframe)

    def dataframe_at_loci(self, contig_name, start, end):
        """
        Subset of entries which overlap an inclusive range of loci
        """
        contig_name = normalize_chromosome(contig_name)
        df_chr = self.dataframe_for_contig(contig_name)

        # find genes whose start/end boundaries overlap with the position
        overlap_start = df_chr.start <= end
        overlap_end = df_chr.end >= start
        overlap = overlap_start & overlap_end
        df_overlap = df_chr[overlap]
        return df_overlap

    def _slice_column(self, col, starts, ends, start, end):
        overlap_start = starts <= end
        overlap_end = ends >= start
        # find genes whose start/end boundaries overlap with the position
        return col[overlap_start & overlap_end]

    def dataframe_column_at_loci(self, contig_name, start, end, col_name):
        """
        Subset of entries which overlap an inclusive range of loci
        """
        contig_name = normalize_chromosome(contig_name)
        df_chr = self.dataframe_for_contig(contig_name)
        assert col_name in df_chr, "Unknown Ensembl property: %s" % col_name
        return self._slice_column(
            df_chr[col_name], df_chr.start, df_chr.end, start, end)

    def dataframe_at_locus(self, contig_name, position):
        return self.dataframe_at_loci(
            contig_name=contig_name,
            start=position,
            stop=position)

    def _db_cols(self, db):
        table_info = db.execute("PRAGMA table_info(ensembl);").fetchall()
        return [info[1] for info in table_info]

    def _db_col_exists(self, db, col_name):
        return col_name in self._db_cols(db)

    def db_column_loci(self, contig_name, start, end, col_name):
        """
        Get the non-null values of a column from the database
        at a particular range of loci
        """
        db = self.db()

        assert self._db_col_exists(db, col_name), \
            "Unknown Ensembl property: %s" % col_name
        query = """
            select distinct %s
            from ensembl
            where
                seqname='%s'
                and start <= %d
                and end >= %d
        """  % (col_name, contig_name, end, start)
        # self.logger.info("Running query: %s" % query)
        results = db.execute(query).fetchall()
        # each result is a tuple, so pull out its first element
        return [result[0] for result in results if result[0] is not None]

    def column_at_loci(self, contig_name, start, end, property_name):
        contig_name = normalize_chromosome(contig_name)
        starts = self._contig_column(contig_name, 'start')
        ends = self._contig_column(contig_name, 'end')
        col = self._contig_column(contig_name, property_name)
        return self._slice_column(col, starts, ends, start, end)

    def _property_values_at_loci(self, contig_name, start, end, property_name):
        col = self.column_at_loci(contig_name, start, end, property_name)
        return list(sorted(set(col)))

    def _property_values_at_locus(self, contig_name, position, property_name):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=position,
            end=position,
            property_name=property_name)

    def gene_ids_at_locus(self, contig_name, position):
        return self._property_values_at_locus(contig_name, position, 'gene_id')

    def gene_names_at_locus(self, contig_name, position):
        return self._property_values_at_locus(
            contig_name, position, 'gene_name')

    def exon_ids_at_locus(self, contig_name, position):
        return self._property_values_at_locus(
            contig_name, position, 'exon_id'
        )

    def transcript_ids_at_locus(self, contig_name, position):
        return self._property_values_at_locus(
            contig_name, position, 'transcript_id')

    def transcript_names_at_locus(self, contig_name, position):
        return self._property_values_at_locus(
            contig_name, position, 'transcript_name')

    def gene_ids_at_loci(self, contig_name, start, end):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=start,
            end=end,
            property_name='gene_id')

    def gene_names_at_loci(self, contig_name, start, end):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=start,
            end=end,
            property_name='gene_name')

    def exon_ids_at_loci(self, contig_name, start, end):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=start,
            end=end,
            property_name='exon_id')

    def transcript_ids_at_loci(self, contig_name, start, end):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=start,
            end=end,
            property_name='transcript_id')

    def transcript_names_at_loci(self, contig_name, start, end):
        return self._property_values_at_loci(
            contig_name=contig_name,
            start=start,
            end=end,
            property_name='transcript_name')

    def _run_query(self, sql, required=False):
        db = self.db()
        cursor = db.execute(sql)
        return cursor.fetchall()

    def _query(
            self,
            select_column_names,
            filter_column,
            filter_value,
            feature,
            distinct=False,
            required=False):
        """
        Run a SQL query against the ensembl sqlite3 database, filtered
        both by the feature type and a user-provided column/value.
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
        results = self._run_query(query, required=required)
        if required:
            assert len(results) > 0, \
                "%s with %s = %s not found" % (
                    feature,
                    filter_column,
                    filter_value
                )
        return results

    def _query_feature(self, cols, feature, distinct=False):
        """
        Run a SQL query against the ensembl sqlite3 database, filtered
        only on the feature type.
        """
        query = "select %s%s from ensembl where feature='%s'" % (
               "distinct " if distinct else "",
               ", ".join(cols),
               feature
        )
        return self._run_query(query)

    def _query_distinct_on_contig(self, col_name, feature, contig_name):
        contig_name = normalize_chromosome(contig_name)
        results = self._query(
            select_column_names=[col_name],
            filter_column="seqname",
            filter_value=contig_name,
            feature=feature,
            distinct=True)
        return [result[0] for result in results if result is not None]

    ###################################################
    #
    #         Locations: (contig, start, stop)
    #
    ###################################################

    def _get_location(self, property_name, property_value, feature_type):
        """
        We're currently assuming every gene/transcript/exon has a
        single location.
        This assumption doesn't hold for >1000 gene names!
        TODO: Return multiple locations.
        """
        results = self._query(
            select_column_names=["seqname", "start", "end", "strand"],
            filter_column=property_name,
            filter_value=property_value,
            feature=feature_type,
            distinct=True,
            required=True)
        return results[0]

    def location_of_gene_name(self, gene_name):
        """
        Given a gene name returns (chromosome, start, stop)
        """
        return self._get_location('gene_name', gene_name, 'gene')

    def location_of_gene_id(self, gene_id):
        """
        Given a gene ID returns (chromosome, start, stop)
        """
        return self._get_location('gene_id', gene_id, 'gene')

    def location_of_transcript_id(self, transcript_id):
        return self._get_location('transcript_id', transcript_id, 'transcript')

    def location_of_exon_id(self, exon_id):
        """
        Given an exon ID returns (chromosome, start, stop)
        """
        return self._get_location('exon_id', exon_id, 'exon')

    ###################################################
    #
    #             Gene Names
    #
    ###################################################

    def _query_gene_name(self, property_name, property_value, feature_type):
        results = self._query(
            select_column_names=["gene_name"],
            filter_col=property_name,
            filter_value=property_value,
            feature=feature_type,
            distinct=True,
            required=True)
        return str(results[0][0])

    def gene_names(self):
        return self._query_feature('gene_name', distinct=True)

    def gene_names_on_contig(self, contig_name):
        """
        Return all genes on given chromosome/contig
        """
        return self._query_distinct_on_contig('gene_name', 'gene', contig_name)

    def gene_name_of_gene_id(self, gene_id):
        return self._query_gene_name("gene_id", gene_id, 'gene')

    def gene_name_of_transcript_id(self, transcript_id):
        return self._query_gene_name("transcript_id", transcript_id, 'transcript')

    def gene_name_of_transcript_name(self, transcript_id):
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
        return self._query_feature('gene_id', distinct=True)

    def gene_ids_on_contig(self, contig_name):
        """
        What are all the gene IDs on a given chromosome/contig?
        """
        return self._query_distinct_on_contig('gene_id', 'gene', contig_name)

    def gene_id_of_gene_name(self, gene_name):
        """
        What's the Ensembl gene ID of the given gene name?
        """
        results = self._query(
            select_column_names=["gene_id"],
            filter_column="gene_name",
            filter_value=gene_name,
            feature="gene",
            required=True)
        assert results, "Gene name not found: %s" % gene_name
        return str(results[0][0])


    ###################################################
    #
    #            Transcript IDs
    #
    ###################################################

    def _query_transcript_ids(self, property_name, value):
        results = self._query(
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
    #                Exon IDs
    #
    ###################################################

    def _query_exon_ids(self, property_name, value):
        results = self._query(
            select_column_names=['exon_id'],
            filter_column=property_name,
            filter_value=value,
            feature='exon',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def exon_ids_of_gene_id(self, gene_id):
        return self._query_exon_ids('gene_id', gene_id)

    def exon_ids_of_gene_name(self, gene_name):
        return self._query_exon_ids('gene_name', gene_name)

    def exon_ids_of_transcript_name(self, transcript_name):
        return self._query_exon_ids('transcript_name', transcript_name)

    def exon_ids_of_transcript_id(self, transcript_id):
        return self._query_exon_ids('transcript_id', transcript_id)

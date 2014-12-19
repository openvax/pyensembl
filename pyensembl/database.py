from os.path import join, exists
import sqlite3

from locus import normalize_chromosome, normalize_strand

import datacache

class Database(object):
    """
    Wrapper around sqlite3 database so that the rest of the
    library doesn't have to worry about constructing the .db file or
    writing SQL queries directly.
    """

    def __init__(self, gtf):
        self.gtf = gtf
        self._connection = None

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
    def connection(self):
        """
        Return the sqlite3 database for this Ensembl release
        (download and/or construct it necessary).
        As a side effect, stores the database connection in self._db
        """
        if self._connection is None:
            self._connection = self._connect_or_create_database()
        return self._connection

    def columns(self):
        sql = "PRAGMA table_info(ensembl);"
        table_info = self.connection.execute(sql).fetchall()
        return [info[1] for info in table_info]

    def column_exists(self, column_name):
        return column_name in self.columns()

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

        if not self.column_exists(column_name):
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

        tuples = self.connection.execute(query, query_params).fetchall()
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

    def run_sql_query(self, sql, required=False, query_params=[]):
        """
        Given an arbitrary SQL query, run it against the Ensembl database
        and return the results.
        """
        cursor = self.connection.execute(sql, query_params)
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
        sql = """
            SELECT %s%s
            FROM ensembl
            WHERE %s = ?
            AND feature= ?
        """ % ("distinct " if distinct else "",
               ", ".join(select_column_names),
               filter_column)
        query_params = [filter_value, feature]
        return self.run_sql_query(
            sql, required=required, query_params=query_params)

    def query_one(
            self,
            select_column_names,
            filter_column,
            filter_value,
            feature,
            distinct=False,
            required=False):
        results = self.query(
            select_column_names,
            filter_column,
            filter_value,
            feature,
            distinct=distinct,
            required=required)

        if len(results) == 0:
            raise ValueError("%s not found: %s" % (filter_column, filter_value))
        elif len(results) > 1:
            raise ValueError(
                "Found multiple entries with %s=%s (%s)" % (
                    filter_column, filter_value, results))
        return results[0]

    def query_feature_values(
            self,
            column,
            feature,
            distinct=True,
            contig=None,
            strand=None):
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

        if strand:
            strand = normalize_strand(strand)
            query += " AND strand = ?"
            query_params.append(strand)

        rows = self.run_sql_query(query, query_params = query_params)
        return [row[0] for row in rows if row is not None]

    def query_distinct_on_contig(self, column_name, feature, contig):
        return self.query_feature_values(
            column_name=column_name,
            feature=feature,
            contig=contig,
            distinct=True)

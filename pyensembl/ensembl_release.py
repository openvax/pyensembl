"""
Contains the EnsemblRelease class, with its millions of accessors and wrappers
around the Ensembl annotation database.
"""
from __future__ import print_function, division, absolute_import

from glob import glob
import logging
from os.path import join
from os import remove

from .compute_cache import cached_object
from .database import Database
from .exon import Exon
from .gene import Gene
from .gtf import GTF
from .reference_transcripts import ReferenceTranscripts
from .release_info import check_release_number, MAX_ENSEMBL_RELEASE
from .transcript import Transcript
from .url_templates import ENSEMBL_FTP_SERVER


class EnsemblRelease(object):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database and provides a wide
    variety of helper methods for accessing this data.
    """

    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 server=ENSEMBL_FTP_SERVER,
                 auto_download=False):
        self.release = check_release_number(release)
        self.species = "homo_sapiens"
        self.server = server
        self.auto_download = auto_download
        self.gtf = GTF(self.release, self.species, server,
                       auto_download=auto_download)
        self.db = Database(gtf=self.gtf, auto_download=auto_download)
        self.reference = ReferenceTranscripts(
            ensembl_release=self.release,
            species=self.species,
            server=server,
            auto_download=auto_download)

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

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
        base = self.gtf.local_filename()
        dirpath = self.gtf.local_gtf_path()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)

    def clear_cache(self):
        self.gtf.clear_cache()
        self.reference.clear_cache()
        self._delete_cached_files()

    def all_feature_values(
            self,
            column,
            feature,
            distinct=True,
            contig=None,
            strand=None):
        """
        Cached lookup of all values for a particular feature property from
        the database, caches repeated queries in memory and
        stores them as a CSV.

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

        strand : str, optional
            Restrict results to '+' or '-' strands

        Returns a list constructed from query results.
        """
        # since we're constructing a list, rather than a DataFrame,
        # we're going to store it using pickling rather than the Pandas
        # CSV serializer. Change the default extension from ".csv" to ".pickle"
        pickle_path = self.gtf.local_data_file_path(
            feature=feature,
            column=column,
            contig=contig,
            strand=strand,
            distinct=distinct,
            extension=".pickle")
        def run_query():
            results = self.db.query_feature_values(
                column=column,
                feature=feature,
                distinct=distinct,
                contig=contig,
                strand=strand)
            assert isinstance(results, list), \
                "Expected list from Database.query_feature_values, got %s" % (
                    type(results))
            return results
        return cached_object(pickle_path, compute_fn=run_query)

    def install(self):
        """
        Explicitely download and index any data for this release that 
        is not yet downloaded and/or indexed.
        """
        self._download(force=False)
        self._index(force=False)

    def index(self):
        """
        Explicitely index all data for this release, regardless of
        whether it is already indexed. Raises an error if data is
        not downloaded.
        """
        self._index(force=True)

    def _index(self, force):
        if self.db.create(force=force):
            logging.info("Annotation data for release %s has just "
                         "been indexed" % self.release)
        else:
            logging.info("Annotation data for release %s is already "
                         "indexed" % self.release)
        if self.reference.index(force=force):
            logging.info("Transcript sequence data for release %s "
                         "has just been indexed" % self.release)
        else:
            logging.info("Transcript sequence data for release %s is "
                         "already indexed" % self.release)

    def download(self):
        """
        Explicitely download all data for this release, regardless of
        whether it is already downloaded.
        """
        self._download(force=True)

    def _download(self, force):
        if self.gtf.download(force=force):
            logging.info("Annotation data for release %s has just "
                         "been downloaded" % self.release)
        else:
            logging.info("Annotation data for release %s is already "
                         "downloaded" % self.release)
        if self.reference.download_transcript_sequences(force=force):
            logging.info("Transcript sequence data for release %s "
                         "has just been downloaded" % self.release)
        else:
            logging.info("Transcript sequence data for release %s is "
                         "already downloaded" % self.release)

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
        return self.db.distinct_column_values_at_locus(
            column='gene_id',
            feature='gene',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def gene_names_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
             column='gene_name',
             feature='gene',
             contig=contig,
             position=position,
             end=end,
             strand=strand)

    def exon_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column='exon_id',
            feature='exon',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column='transcript_id',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_names_at_locus(
            self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column='transcript_name',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def protein_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column='protein_id',
            feature='transcript',
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    ###################################################
    #
    #         Methods which return Locus objects
    #         containing (contig, start, stop, strand)
    #         of various genomic entities
    #
    ###################################################

    def locus_of_gene_id(self, gene_id):
        """
        Given a gene ID returns Locus with: chromosome, start, stop, strand
        """
        return self.db.query_locus(
            filter_column='gene_id',
            filter_value=gene_id,
            feature='gene')

    def loci_of_gene_names(self, gene_name):
        """
        Given a gene name returns list of Locus objects with fields:
            chromosome, start, stop, strand
        You can get multiple results since a gene might have multiple copies
        in the genome.
        """
        return self.db.query_loci('gene_name', gene_name, 'gene')

    def locus_of_transcript_id(self, transcript_id):
        return self.db.query_locus(
            filter_column='transcript_id',
            filter_value=transcript_id,
            feature='transcript')

    def locus_of_exon_id(self, exon_id):
        """
        Given an exon ID returns Locus
        """
        return self.db.query_locus('exon_id', exon_id, feature='exon')

    ###################################################
    #
    #             Gene Info Objects
    #
    ###################################################

    def genes(self, contig=None, strand=None):
        """
        Returns all Gene objects in Ensembl. Can be restricted to a
        particular contig/chromosome and strand by the following arguments:

        Parameters
        ----------
        contig : str
            Only return genes on the given contig.

        strand : str
            Only return genes on this strand.
        """
        gene_ids = self.gene_ids(contig=contig, strand=strand)
        return [self.gene_by_id(gene_id) for gene_id in gene_ids]

    def gene_by_id(self, gene_id):
        """
        Construct a Gene object for the given gene ID.
        """
        return Gene(gene_id, self.db, self.reference)

    def genes_by_name(self, gene_name):
        """
        Get all the unqiue genes with the given name (there might be multiple
        due to copies in the genome), return a list containing a Gene object
        for each distinct ID.
        """
        gene_ids = self.gene_ids_of_gene_name(gene_name)
        return [self.gene_by_id(gene_id) for gene_id in gene_ids]

    def gene_by_protein_id(self, protein_id):
        """
        Get the gene ID associated with the given protein ID,
        return its Gene object
        """
        gene_id = self.gene_id_of_protein_id(protein_id)
        return self.gene_by_id(gene_id)

    ###################################################
    #
    #             Gene Names
    #
    ###################################################

    def _query_gene_name(self, property_name, property_value, feature_type):
        results = self.db.query(
            select_column_names=["gene_name"],
            filter_column=property_name,
            filter_value=property_value,
            feature=feature_type,
            distinct=True,
            required=True)
        return str(results[0][0])

    def gene_names(self, contig=None, strand=None):
        """
        Return all genes in the database,
        optionally restrict to a chromosome and/or strand.
        """
        return self.all_feature_values(
            column='gene_name',
            feature='gene',
            contig=contig,
            strand=strand)

    def gene_name_of_gene_id(self, gene_id):
        return self._query_gene_name("gene_id", gene_id, 'gene')

    def gene_name_of_transcript_id(self, transcript_id):
        return self._query_gene_name(
            "transcript_id", transcript_id, 'transcript')

    def gene_name_of_transcript_name(self, transcript_name):
        return self._query_gene_name(
            "transcript_name", transcript_name, 'transcript')

    def gene_name_of_exon_id(self, exon_id):
        return self._query_gene_name("exon_id", exon_id, 'exon')

    ###################################################
    #
    #             Gene IDs
    #
    ###################################################

    def _query_gene_ids(self, property_name, value, feature='gene'):
        results = self.db.query(
            select_column_names=['gene_id'],
            filter_column=property_name,
            filter_value=value,
            feature=feature,
            distinct=True,
            required=True)
        return [str(result[0]) for result in results if result[0]]

    def gene_ids(self, contig=None, strand=None):
        """
        What are all the gene IDs
        (optionally restrict to a given chromosome/contig and/or strand)
        """
        return self.all_feature_values(
            column='gene_id',
            feature='gene',
            contig=contig,
            strand=strand)

    def gene_ids_of_gene_name(self, gene_name):
        """
        What are the Ensembl gene IDs associated with a given gene name?
        (due to copy events, there might be multiple genes per name)
        """
        results = self._query_gene_ids('gene_name', gene_name)
        if len(results) == 0:
            raise ValueError("Gene name not found: %s" % gene_name)
        return results

    def gene_id_of_protein_id(self, protein_id):
        """
        What is the Ensembl gene ID associated with a given protein ID?
        """
        results = self._query_gene_ids('protein_id', protein_id,
                                       feature='CDS')
        if len(results) == 0:
            raise ValueError("Protein ID not found: %s" % protein_id)
        assert len(results) == 1, \
            ("Should have only one gene ID for a given protein ID, "
             "but found %d: %s" % (len(results), results))
        return results[0]

    ###################################################
    #
    #             Transcript Info Objects
    #
    ###################################################

    def transcripts(self, contig=None, strand=None):
        """
        Construct Transcript object for every transcript entry in
        the Ensembl database. Optionally restrict to a particular
        chromosome using the `contig` argument.
        """
        transcript_ids = self.transcript_ids(contig=contig, strand=strand)
        return [
            Transcript(transcript_id, self.db, self.reference)
            for transcript_id in transcript_ids
        ]

    def transcript_by_id(self, transcript_id):
        return Transcript(transcript_id, self.db, self.reference)

    def transcripts_by_name(self, transcript_name):
        transcript_ids = self.transcript_ids_of_transcript_name(transcript_name)
        return [
            self.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    def transcript_by_protein_id(self, protein_id):
        transcript_id = self.transcript_id_of_protein_id(protein_id)
        return self.transcript_by_id(transcript_id)

    ###################################################
    #
    #            Transcript Names
    #
    ###################################################

    def _query_transcript_names(self, property_name, value):
        results = self.db.query(
            select_column_names=['transcript_name'],
            filter_column=property_name,
            filter_value=value,
            feature='transcript',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def transcript_names(self, contig=None, strand=None):
        """
        What are all the transcript names in the database
        (optionally, restrict to a given chromosome and/or strand)
        """
        return self.all_feature_values(
            column='transcript_name',
            feature='transcript',
            contig=contig,
            strand=strand)

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

    def _query_transcript_ids(self, property_name, value,
                              feature='transcript'):
        results = self.db.query(
            select_column_names=['transcript_id'],
            filter_column=property_name,
            filter_value=value,
            feature=feature,
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def transcript_ids(self, contig=None, strand=None):
        return self.all_feature_values(
            column='transcript_id',
            feature='transcript',
            contig=contig,
            strand=strand)

    def transcript_ids_of_gene_id(self, gene_id):
        return self._query_transcript_ids('gene_id', gene_id)

    def transcript_ids_of_gene_name(self, gene_name):
        return self._query_transcript_ids('gene_name', gene_name)

    def transcript_ids_of_transcript_name(self, transcript_name):
        return self._query_transcript_ids('transcript_name', transcript_name)

    def transcript_ids_of_exon_id(self, exon_id):
        return self._query_transcript_ids('exon_id', exon_id)

    def transcript_id_of_protein_id(self, protein_id):
        """
        What is the Ensembl transcript ID associated with a given
        protein ID?
        """
        results = self._query_transcript_ids('protein_id', protein_id,
                                             feature='CDS')
        if len(results) == 0:
            raise ValueError("Protein ID not found: %s" % protein_id)
        assert len(results) == 1, \
            ("Should have only one transcript ID for a given protein ID, "
             "but found %d: %s" % (len(results), results))
        return results[0]

    ###################################################
    #
    #             Exon Info Objects
    #
    ###################################################

    def exons(self, contig=None, strand=None):
        """
        Create exon object for all exons in the database, optionally
        restrict to a particular chromosome using the `contig` argument.
        """
        # DataFrame with single column called 'exon_id'
        exon_ids = self.exon_ids(contig=contig, strand=strand)
        return [Exon(exon_id, self.db) for exon_id in exon_ids]

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
        results = self.db.query(
            select_column_names=['exon_id'],
            filter_column=property_name,
            filter_value=value,
            feature='exon',
            distinct=True,
            required=True)
        return [result[0] for result in results]

    def exon_ids(self, contig=None, strand=None):
        return self.all_feature_values(
            column='exon_id',
            feature='exon',
            contig=contig,
            strand=strand)

    def exon_ids_of_gene_id(self, gene_id):
        return self._query_exon_ids('gene_id', gene_id)

    def exon_ids_of_gene_name(self, gene_name):
        return self._query_exon_ids('gene_name', gene_name)

    def exon_ids_of_transcript_name(self, transcript_name):
        return self._query_exon_ids('transcript_name', transcript_name)

    def exon_ids_of_transcript_id(self, transcript_id):
        return self._query_exon_ids('transcript_id', transcript_id)


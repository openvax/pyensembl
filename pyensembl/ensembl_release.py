"""
Contains the EnsemblRelease class, with its millions of accessors and wrappers
around the Ensembl annotation database.
"""

from glob import glob
import logging
from os.path import join
from os import remove

from common import CACHE_SUBDIR
from compute_cache import cached_object
from database import Database
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

    def __init__(self, release, server=ENSEMBL_FTP_SERVER):
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.release = check_release_number(release)
        self.species = "homo_sapiens"
        self.server = server
        self.gtf = GTF(self.release, self.species, server)
        self.db = Database(gtf = self.gtf)
        self.reference = ReferenceTranscripts(
            self.release, self.species, server)

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
        base = self.gtf.base_filename()
        dirpath = self.local_gtf_dir()
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
    #         Locations: (contig, start, stop)
    #
    ###################################################

    def _get_locations(self, property_name, property_value, feature):
        return self.db.query(
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
        """
        Construct a Gene object for the given gene ID.
        """
        return Gene(gene_id, self.db)

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

    def gene_name_of_exon_id(self, transcript_id):
        return self._query_gene_name("exon_id", exon_id, 'exon')

    ###################################################
    #
    #             Gene IDs
    #
    ###################################################

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
        results = self.db.query(
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

    def _query_transcript_ids(self, property_name, value):
        results = self.db.query(
            select_column_names=['transcript_id'],
            filter_column=property_name,
            filter_value=value,
            feature='transcript',
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


    ###################################################
    #
    #                Start Codons
    #
    ###################################################

    def _query_start_codon_location(self, property_name, value):
        return self.db.query(
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
        return self.db.query(
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


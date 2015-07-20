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

"""
Contains the Genome class, with its millions of accessors and wrappers
around an arbitrary genomic database.
"""

from __future__ import print_function, division, absolute_import

from glob import glob
import logging
from os.path import join
from os import remove

from .common import (
    require_human_transcript_id,
    require_human_protein_id,
    memoize
)
from .compute_cache import cached_object
from .database import Database
from .exon import Exon
from .gene import Gene
from .gtf import GTF
from .sequence_data import SequenceData
from .transcript import Transcript

class Genome(object):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular genomic database source (e.g. a single Ensembl release) and
    provides a wide variety of helper methods for accessing this data.
    """
    def __init__(self,
                 reference_name,
                 genome_source,
                 name=None,
                 version=None,
                 only_human=False,
                 auto_download=False,
                 local_fasta_filename_func=None,
                 require_ensembl_ids=True):
        self.reference_name = reference_name
        self.genome_source = genome_source
        self.name = name
        self.version = version
        self.only_human = only_human
        self.auto_download = auto_download

        # GTF object wraps the source GTF file from which we get
        # genome annotations. Presents access to each feature
        # annotations as a pandas.DataFrame.
        self.gtf = GTF(
            genome_source,
            auto_download=auto_download)

        # Database object turns the GTF dataframes into sqlite3 tables
        # and wraps them with methods like `query_one`
        self.db = Database(gtf=self.gtf, auto_download=auto_download)

        # get the path for the cDNA FASTA file containing
        # this genome database's transcript sequences
        transcript_sequences = None
        if genome_source.transcript_fasta_path:
            transcript_sequences = SequenceData(
                genome_source=genome_source,
                fasta_type="transcript",
                local_filename_func=local_fasta_filename_func,
                require_ensembl_ids=require_ensembl_ids,
                auto_download=auto_download)
        self.transcript_sequences = transcript_sequences

        protein_sequences = None
        if genome_source.protein_fasta_path:
            protein_sequences = SequenceData(
                genome_source=genome_source,
                fasta_type="protein",
                local_filename_func=local_fasta_filename_func,
                require_ensembl_ids=require_ensembl_ids,
                auto_download=auto_download)
        self.protein_sequences = protein_sequences

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

    def __str__(self):
        return ("Genome(name=%s, version=%s, reference_name=%s, "
                "only_human=%s, genome_source=%s)") % (
            self.name,
            self.version,
            self.reference_name,
            self.only_human,
            self.genome_source)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is Genome and
            self.name == other.name and
            self.version == other.version and
            self.only_human == other.only_human and
            self.genome_source == other.genome_source)

    def __hash__(self):
        return hash((self.name, self.version, self.only_human,
                     self.genome_source))

    def _delete_cached_files(self):
        """
        Any files which start with the same name as our GTF file
        is assumed to be some view of the this genome database's data
        and thus safe to delete.
        """
        base = self.gtf.base_filename()
        dirpath = self.gtf.local_gtf_path()
        for path in glob(join(dirpath, base + "*")):
            logging.info("Deleting cached file %s", path)
            remove(path)

    def clear_cache(self):
        for maybe_fn in self.__dict__.values():
            # clear cache associated with all memoization decorators,
            # GTF and SequenceData objects
            if hasattr(maybe_fn, "clear_cache"):
                maybe_fn.clear_cache()

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
            Restrict results to "+" or "-" strands

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
        Explicitely download and index any data for this genome source that
        is not yet downloaded and/or indexed.
        """
        self.download(force=False)
        self.index(force=False)

    def index(self, force=True):
        """
        Create databases and indices for all data for this genome source

        Parameters
        ----------
        force : bool
            Recreate databases even if already indexed (default = True)

        Raises an error if data is not downloaded.
        """
        self.db.create(force=force)
        if self.transcript_sequences:
            self.transcript_sequences.index(force=force)
        if self.protein_sequences:
            self.protein_sequences.index(force=force)

    def transcript_sequence(self, transcript_id):
        """Return cDNA nucleotide sequence of transcript, or None if
        transcript doesn't have cDNA sequence.
        """
        assert self.transcript_sequences, (
            "This genome source does not include transcript FASTA data: %s"
            % self.genome_source)
        if self.only_human:
            require_human_transcript_id(transcript_id)
        return self.transcript_sequences.get(transcript_id)

    def protein_sequence(self, protein_id):
        """Return cDNA nucleotide sequence of transcript, or None if
        transcript doesn't have cDNA sequence.
        """
        assert self.protein_sequences, (
            "This genome source does not include protein FASTA data: %s"
            % self.genome_source)
        if self.only_human:
            require_human_protein_id(protein_id)
        return self.protein_sequences.get(protein_id)

    def download(self, force=True):
        """
        Download all data for this genome source.

        Parameters
        ----------
        force : bool
            Download data even if we already have a local copy.
        """
        self.gtf.download(force=force)
        if self.transcript_sequences:
            self.transcript_sequences.download(force=force)
        if self.protein_sequences:
            self.protein_sequences.download(force=force)

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
            column="gene_id",
            feature="gene",
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def gene_names_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
             column="gene_name",
             feature="gene",
             contig=contig,
             position=position,
             end=end,
             strand=strand)

    def exon_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column="exon_id",
            feature="exon",
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column="transcript_id",
            feature="transcript",
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def transcript_names_at_locus(
            self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column="transcript_name",
            feature="transcript",
            contig=contig,
            position=position,
            end=end,
            strand=strand)

    def protein_ids_at_locus(self, contig, position, end=None, strand=None):
        return self.db.distinct_column_values_at_locus(
            column="protein_id",
            feature="transcript",
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

    @memoize
    def locus_of_gene_id(self, gene_id):
        """
        Given a gene ID returns Locus with: chromosome, start, stop, strand
        """
        return self.db.query_locus(
            filter_column="gene_id",
            filter_value=gene_id,
            feature="gene")

    @memoize
    def loci_of_gene_names(self, gene_name):
        """
        Given a gene name returns list of Locus objects with fields:
            chromosome, start, stop, strand
        You can get multiple results since a gene might have multiple copies
        in the genome.
        """
        return self.db.query_loci("gene_name", gene_name, "gene")

    @memoize
    def locus_of_transcript_id(self, transcript_id):
        return self.db.query_locus(
            filter_column="transcript_id",
            filter_value=transcript_id,
            feature="transcript")

    @memoize
    def locus_of_exon_id(self, exon_id):
        """
        Given an exon ID returns Locus
        """
        return self.db.query_locus("exon_id", exon_id, feature="exon")

    ###################################################
    #
    #             Gene Info Objects
    #
    ###################################################

    @memoize
    def genes(self, contig=None, strand=None):
        """
        Returns all Gene objects in the database. Can be restricted to a
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

    @memoize
    def gene_by_id(self, gene_id):
        """
        Construct a Gene object for the given gene ID.
        """
        field_names = [
            "seqname",
            "start",
            "end",
            "strand",
        ]
        optional_field_names = [
            "gene_name",
            "gene_biotype",
        ]
        # Do not look for gene_name and gene_biotype if they are
        # not in the database.
        field_names.extend([name for name in optional_field_names
                            if self.db.column_exists("gene", name)])
        result = self.db.query_one(
            field_names,
            filter_column="gene_id",
            filter_value=gene_id,
            feature="gene")
        if not result:
            raise ValueError("Gene not found: %s" % (gene_id,))

        gene_name, gene_biotype = None, None
        assert len(result) >= 4 and len(result) <= 6, \
            "Result is not the expected length: %d" % len(result)
        contig, start, end, strand = result[:4]
        if len(result) == 5:
            if "gene_name" in field_names:
                gene_name = result[4]
            else:
                gene_biotype = result[4]
        elif len(result) == 6:
            gene_name, gene_biotype = result[4:]

        return Gene(
            gene_id=gene_id,
            gene_name=gene_name,
            contig=contig,
            start=start,
            end=end,
            strand=strand,
            biotype=gene_biotype,
            ensembl=self,
            require_valid_biotype=("gene_biotype" in field_names))

    @memoize
    def genes_by_name(self, gene_name):
        """
        Get all the unqiue genes with the given name (there might be multiple
        due to copies in the genome), return a list containing a Gene object
        for each distinct ID.
        """
        gene_ids = self.gene_ids_of_gene_name(gene_name)
        return [self.gene_by_id(gene_id) for gene_id in gene_ids]

    @memoize
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

    @memoize
    def _query_gene_name(self, property_name, property_value, feature_type):
        results = self.db.query(
            select_column_names=["gene_name"],
            filter_column=property_name,
            filter_value=property_value,
            feature=feature_type,
            distinct=True,
            required=True)
        return str(results[0][0])

    @memoize
    def gene_names(self, contig=None, strand=None):
        """
        Return all genes in the database,
        optionally restrict to a chromosome and/or strand.
        """
        return self.all_feature_values(
            column="gene_name",
            feature="gene",
            contig=contig,
            strand=strand)

    @memoize
    def gene_name_of_gene_id(self, gene_id):
        return self._query_gene_name("gene_id", gene_id, "gene")

    @memoize
    def gene_name_of_transcript_id(self, transcript_id):
        return self._query_gene_name(
            "transcript_id", transcript_id, "transcript")

    @memoize
    def gene_name_of_transcript_name(self, transcript_name):
        return self._query_gene_name(
            "transcript_name", transcript_name, "transcript")

    @memoize
    def gene_name_of_exon_id(self, exon_id):
        return self._query_gene_name("exon_id", exon_id, "exon")

    ###################################################
    #
    #             Gene IDs
    #
    ###################################################

    @memoize
    def _query_gene_ids(self, property_name, value, feature="gene"):
        results = self.db.query(
            select_column_names=["gene_id"],
            filter_column=property_name,
            filter_value=value,
            feature=feature,
            distinct=True,
            required=True)
        return [str(result[0]) for result in results if result[0]]

    @memoize
    def gene_ids(self, contig=None, strand=None):
        """
        What are all the gene IDs
        (optionally restrict to a given chromosome/contig and/or strand)
        """
        return self.all_feature_values(
            column="gene_id",
            feature="gene",
            contig=contig,
            strand=strand)

    @memoize
    def gene_ids_of_gene_name(self, gene_name):
        """
        What are the gene IDs associated with a given gene name?
        (due to copy events, there might be multiple genes per name)
        """
        results = self._query_gene_ids("gene_name", gene_name)
        if len(results) == 0:
            raise ValueError("Gene name not found: %s" % gene_name)
        return results

    @memoize
    def gene_id_of_protein_id(self, protein_id):
        """
        What is the gene ID associated with a given protein ID?
        """
        results = self._query_gene_ids("protein_id", protein_id,
                                       feature="CDS")
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

    @memoize
    def transcripts(self, contig=None, strand=None):
        """
        Construct Transcript object for every transcript entry in
        the database. Optionally restrict to a particular
        chromosome using the `contig` argument.
        """
        transcript_ids = self.transcript_ids(contig=contig, strand=strand)
        return [
            self.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    @memoize
    def transcript_by_id(self, transcript_id):
        """Construct Transcript object with given transcript ID"""

        optional_field_names = [
            "transcript_name",
            "transcript_biotype",
        ]
        field_names = [
            "seqname",
            "start",
            "end",
            "strand",
            "gene_id",
        ]
        # Do not look for transcript_name and transcript_biotype if
        # they are not in the database.
        field_names.extend([name for name in optional_field_names
                            if self.db.column_exists("transcript", name)])
        result = self.db.query_one(
            select_column_names=field_names,
            filter_column="transcript_id",
            filter_value=transcript_id,
            feature="transcript",
            distinct=True)
        if not result:
            raise ValueError("Transcript not found: %s" % (transcript_id,))

        transcript_name, transcript_biotype = None, None
        assert len(result) >= 5 and len(result) <= 7, \
            "Result is not the expected length: %d" % len(result)
        contig, start, end, strand, gene_id = result[:5]
        if len(result) == 6:
            if "transcript_name" in field_names:
                transcript_name = result[5]
            else:
                transcript_biotype = result[5]
        elif len(result) == 7:
            transcript_name, transcript_biotype = result[5:]

        return Transcript(
            transcript_id=transcript_id,
            transcript_name=transcript_name,
            contig=contig,
            start=start,
            end=end,
            strand=strand,
            biotype=transcript_biotype,
            gene_id=gene_id,
            ensembl=self,
            require_valid_biotype=("transcript_biotype" in field_names))

    @memoize
    def transcripts_by_name(self, transcript_name):
        transcript_ids = self.transcript_ids_of_transcript_name(transcript_name)
        return [
            self.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    @memoize
    def transcript_by_protein_id(self, protein_id):
        transcript_id = self.transcript_id_of_protein_id(protein_id)
        return self.transcript_by_id(transcript_id)

    ###################################################
    #
    #            Transcript Names
    #
    ###################################################

    @memoize
    def _query_transcript_names(self, property_name, value):
        results = self.db.query(
            select_column_names=["transcript_name"],
            filter_column=property_name,
            filter_value=value,
            feature="transcript",
            distinct=True,
            required=True)
        return [result[0] for result in results]

    @memoize
    def transcript_names(self, contig=None, strand=None):
        """
        What are all the transcript names in the database
        (optionally, restrict to a given chromosome and/or strand)
        """
        return self.all_feature_values(
            column="transcript_name",
            feature="transcript",
            contig=contig,
            strand=strand)

    @memoize
    def transcript_names_of_gene_name(self, gene_name):
        return self._query_transcript_names("gene_name", gene_name)

    @memoize
    def transcript_name_of_transcript_id(self, transcript_id):
        transcript_names = self._query_transcript_names(
            "transcript_id", transcript_id)
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

    @memoize
    def _query_transcript_ids(
            self,
            property_name,
            value,
            feature="transcript"):
        results = self.db.query(
            select_column_names=["transcript_id"],
            filter_column=property_name,
            filter_value=value,
            feature=feature,
            distinct=True,
            required=True)
        return [result[0] for result in results]

    @memoize
    def transcript_ids(self, contig=None, strand=None):
        return self.all_feature_values(
            column="transcript_id",
            feature="transcript",
            contig=contig,
            strand=strand)

    @memoize
    def transcript_ids_of_gene_id(self, gene_id):
        return self._query_transcript_ids("gene_id", gene_id)

    @memoize
    def transcript_ids_of_gene_name(self, gene_name):
        return self._query_transcript_ids("gene_name", gene_name)

    @memoize
    def transcript_ids_of_transcript_name(self, transcript_name):
        return self._query_transcript_ids("transcript_name", transcript_name)

    @memoize
    def transcript_ids_of_exon_id(self, exon_id):
        return self._query_transcript_ids("exon_id", exon_id)

    @memoize
    def transcript_id_of_protein_id(self, protein_id):
        """
        What is the transcript ID associated with a given protein ID?
        """
        results = self._query_transcript_ids("protein_id", protein_id,
                                             feature="CDS")
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

    @memoize
    def exons(self, contig=None, strand=None):
        """
        Create exon object for all exons in the database, optionally
        restrict to a particular chromosome using the `contig` argument.
        """
        # DataFrame with single column called "exon_id"
        exon_ids = self.exon_ids(contig=contig, strand=strand)
        return [
            self.exon_by_id(exon_id)
            for exon_id in exon_ids
        ]

    @memoize
    def exon_by_id(self, exon_id):
        """Construct an Exon object from its ID by looking up the exon"s
        properties in the given Database.
        """
        field_names = [
            "seqname",
            "start",
            "end",
            "strand",
            "gene_name",
            "gene_id",
        ]

        contig, start, end, strand, gene_name, gene_id = self.db.query_one(
            select_column_names=field_names,
            filter_column="exon_id",
            filter_value=exon_id,
            feature="exon",
            distinct=True)

        return Exon(
            exon_id=exon_id,
            contig=contig,
            start=start,
            end=end,
            strand=strand,
            gene_name=gene_name,
            gene_id=gene_id)

    ###################################################
    #
    #                Exon IDs
    #
    ###################################################

    @memoize
    def _query_exon_ids(self, property_name, value):
        results = self.db.query(
            select_column_names=["exon_id"],
            filter_column=property_name,
            filter_value=value,
            feature="exon",
            distinct=True,
            required=True)
        return [result[0] for result in results]

    @memoize
    def exon_ids(self, contig=None, strand=None):
        return self.all_feature_values(
            column="exon_id",
            feature="exon",
            contig=contig,
            strand=strand)

    @memoize
    def exon_ids_of_gene_id(self, gene_id):
        return self._query_exon_ids("gene_id", gene_id)

    @memoize
    def exon_ids_of_gene_name(self, gene_name):
        return self._query_exon_ids("gene_name", gene_name)

    @memoize
    def exon_ids_of_transcript_name(self, transcript_name):
        return self._query_exon_ids("transcript_name", transcript_name)

    @memoize
    def exon_ids_of_transcript_id(self, transcript_id):
        return self._query_exon_ids("transcript_id", transcript_id)

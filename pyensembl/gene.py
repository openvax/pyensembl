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

from memoized_property import memoized_property
from typechecks import require_string, require_instance

from .biotypes import is_valid_biotype
from .database import Database
from .locus import Locus
from .transcript import Transcript

class Gene(Locus):

    def __init__(self, gene_id, ensembl):
        require_string(gene_id, "gene ID")
        self.id = gene_id
        # can't check the type of ensembl since it will create a circular
        # dependency between this module and ensembl_release but note that
        # it's expected to be an EnsemblRelease object
        self.ensembl = ensembl
        self.db = ensembl.db
        require_instance(self.db, Database, "db")
        columns = [
            'gene_name',
            'seqname',
            'start',
            'end',
            'strand',
            'gene_biotype'
        ]
        gene_name, contig, start, end, strand, biotype = self.db.query_one(
            columns,
            filter_column='gene_id',
            filter_value=gene_id,
            feature='gene')
        if not gene_name:
            raise ValueError("Missing name for gene with ID = %s" % gene_id)
        self.name = gene_name

        Locus.__init__(self, contig, start, end, strand)

        if not is_valid_biotype(biotype):
            raise ValueError(
                "Invalid gene_biotype %s for gene with ID = %s" % (
                    biotype, gene_id))
        elif biotype:
            self.biotype = biotype
        else:
            raise ValueError(
                "Missing gene_biotype for gene with ID = %s" % gene_id)

    def __str__(self):
        return "Gene(id=%s, name=%s, location=%s:%d-%d)" % (
            self.id,
            self.name,
            self.contig,
            self.start,
            self.end)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            isinstance(other, Gene) and
            self.id == other.id and
            self.ensembl == other.ensembl)

    def __hash__(self):
        return hash(self.id)

    @memoized_property
    def transcripts(self):
        """
        Property which dynamically construct transcript objects for all
        transcript IDs associated with this gene.
        """
        transcript_id_results = self.db.query(
            select_column_names=['transcript_id'],
            filter_column='gene_id',
            filter_value=self.id,
            feature='transcript',
            distinct=False,
            required=False)

        # We're doing a SQL query for each transcript ID to fetch
        # its particular information, might be more efficient if we
        # just get all the columns here, but how do we keep that modular?
        return [
            Transcript(result[0], self.ensembl)
            for result in transcript_id_results
        ]

    @memoized_property
    def exons(self):
        exons_dict = {}
        for transcript in self.transcripts:
            for exon in transcript.exons:
                if exon.id not in exons_dict:
                    exons_dict[exon.id] = exon
        return list(exons_dict.values())

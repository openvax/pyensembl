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
from typechecks import require_instance

from .biotypes import is_valid_biotype
from .database import Database
from .locus import Locus

class Gene(Locus):

    def __init__(
            self,
            gene_id,
            gene_name,
            contig,
            start,
            end,
            strand,
            biotype,
            genome,
            require_valid_biotype=True):
        """
        genome is a Genome object.
        """
        self.id = gene_id
        self.genome = genome
        self.db = genome.db
        require_instance(self.db, Database, "db")

        self.name = gene_name

        Locus.__init__(self, contig, start, end, strand)

        if require_valid_biotype and not is_valid_biotype(biotype):
            raise ValueError(
                "Invalid gene_biotype %s for gene with ID = %s" % (
                    biotype, gene_id))
        self.biotype = biotype

    def __str__(self):
        return "Gene(id=%s, name=%s, biotype=%s, location=%s:%d-%d)" % (
            self.id,
            self.name,
            self.biotype,
            self.contig,
            self.start,
            self.end)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is Gene and
            self.id == other.id and
            self.genome == other.genome)

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
            self.genome.transcript_by_id(result[0])
            for result in transcript_id_results
        ]

    @memoized_property
    def exons(self):
        exon_set = set([])
        for transcript in self.transcripts:
            for exon in transcript.exons:
                exon_set.add(exon)
        return list(sorted(exon_set))

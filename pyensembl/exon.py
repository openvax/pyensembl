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

from .locus import Locus

class Exon(Locus):
    def __init__(
            self,
            exon_id,
            contig,
            start,
            end,
            strand,
            gene_name,
            gene_id):
        Locus.__init__(self, contig, start, end, strand)
        self.id = exon_id
        self.gene_name = gene_name
        self.gene_id = gene_id

    def __str__(self):
        return "Exon(exon_id=%s, gene_name=%s, contig=%s, start=%d, end=%s)" % (
            self.id, self.gene_name, self.contig, self.start, self.end)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return (
            other.__class__ is Exon and
            self.contig == other.contig and
            self.start == other.start and
            self.end == other.end and
            self.strand == other.strand and
            self.id == other.id)

# Copyright (c) 2015-2016. Mount Sinai School of Medicine
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

from .locus import Locus

class LocusWithGenome(Locus):
    """
    Common base class for Gene and Transcript to avoid copying
    their shared logic.
    """
    def __init__(self, contig, start, end, strand, genome):
        Locus.__init__(self, contig, start, end, strand)
        self.genome = genome
        self.db = self.genome.db

    def to_dict(self):
        return dict(
            contig=self.contig,
            start=self.start,
            end=self.end,
            strand=self.strand,
            genome=self.genome)

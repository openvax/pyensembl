# Copyright (c) 2016. Mount Sinai School of Medicine
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

from __future__ import absolute_import
import pickle
from nose.tools import eq_, ok_

from .common import test_ensembl_releases
from .data import TP53_gene_id

@test_ensembl_releases()
def test_gene(ensembl):
    gene = ensembl.gene_by_id(TP53_gene_id)
    gene_pickled = pickle.dumps(gene)
    gene_new = pickle.loads(gene_pickled)
    eq_(gene, gene_new)
    ok_(gene.db is not None)

@test_ensembl_releases()
def test_transcript(ensembl):
    gene = ensembl.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    transcript_pickled = pickle.dumps(transcript)
    transcript_new = pickle.loads(transcript_pickled)
    eq_(transcript, transcript_new)
    ok_(transcript.db is not None)

@test_ensembl_releases()
def test_genome(ensembl):
    gene = ensembl.gene_by_id(TP53_gene_id)
    genome = gene.genome
    genome_pickled = pickle.dumps(genome)
    genome_new = pickle.loads(genome_pickled)
    eq_(genome, genome_new)
    ok_(genome.db is not None)

    # This Genome happens to be an EnsemblRelease; test that too.
    eq_(genome.release, genome_new.release)
    eq_(genome.species, genome_new.species)

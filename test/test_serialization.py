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
from nose.tools import eq_, ok_, with_setup
from pyensembl import Genome

from .common import test_ensembl_releases
from .data import (
    TP53_gene_id,
    custom_mouse_genome_grcm38_subset,
    setup_init_custom_mouse_genome
)


@test_ensembl_releases()
def test_pickle_ensembl_gene(ensembl_genome):
    gene = ensembl_genome.gene_by_id(TP53_gene_id)
    gene_pickled = pickle.dumps(gene)
    gene_new = pickle.loads(gene_pickled)
    eq_(gene, gene_new)
    ok_(gene.db is not None)

@test_ensembl_releases()
def test_pickle_ensembl_transcript(ensembl_genome):
    gene = ensembl_genome.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    transcript_pickled = pickle.dumps(transcript)
    transcript_reconstructed = pickle.loads(transcript_pickled)
    eq_(transcript, transcript_reconstructed)
    ok_(transcript.db is not None)

@test_ensembl_releases()
def test_pickle_ensembl_genome(ensembl_genome):
    genome_pickled = pickle.dumps(ensembl_genome)
    genome_reconstructed = pickle.loads(genome_pickled)
    eq_(ensembl_genome, genome_reconstructed)

    # This Genome happens to be an EnsemblRelease; test that too.
    eq_(ensembl_genome.release, genome_reconstructed.release)
    eq_(ensembl_genome.species, genome_reconstructed.species)

@test_ensembl_releases()
def test_ensembl_genome_to_dict(ensembl_genome):
    genome_dict = ensembl_genome.to_dict()
    genome_reconstructed = ensembl_genome.__class__.from_dict(genome_dict)
    eq_(ensembl_genome, genome_reconstructed)

@test_ensembl_releases()
def test_ensembl_genome_to_json(ensembl_genome):
    genome_json = ensembl_genome.to_json()
    genome_reconstructed = ensembl_genome.__class__.from_json(genome_json)
    eq_(ensembl_genome, genome_reconstructed)


@with_setup(setup=setup_init_custom_mouse_genome)
def test_custom_genome_to_json():
    json = custom_mouse_genome_grcm38_subset.to_json()
    reconstructed = Genome.from_json(json)
    eq_(custom_mouse_genome_grcm38_subset, reconstructed)

@with_setup(setup=setup_init_custom_mouse_genome)
def test_custom_genome_to_dict():
    reconstructed = Genome.from_dict(custom_mouse_genome_grcm38_subset.to_dict())
    eq_(custom_mouse_genome_grcm38_subset, reconstructed)

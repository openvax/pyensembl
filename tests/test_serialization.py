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

import pickle

from pyensembl import Genome, Transcript, Gene, Exon
from pyensembl.species import Species, human

from .common import run_multiple_genomes, eq_
from .data import (
    TP53_gene_id,
    custom_mouse_genome_grcm38_subset,
    setup_init_custom_mouse_genome,
)


@run_multiple_genomes
def test_pickle_ensembl_gene(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    gene_new = pickle.loads(pickle.dumps(gene))
    assert gene == gene_new


@run_multiple_genomes()
def test_pickle_ensembl_transcript(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    transcript_reconstructed = pickle.loads(pickle.dumps(transcript))
    eq_(transcript, transcript_reconstructed)


@run_multiple_genomes()
def test_pickle_ensembl_exon(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    exon = transcript.exons[0]
    exon_reconstructed = pickle.loads(pickle.dumps(exon))
    eq_(exon, exon_reconstructed)


@run_multiple_genomes()
def test_json_ensembl_gene(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    gene_reconstructed = Gene.from_json(gene.to_json())
    eq_(gene, gene_reconstructed)


@run_multiple_genomes()
def test_json_ensembl_transcript(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    transcript_reconstructed = Transcript.from_json(transcript.to_json())
    eq_(transcript, transcript_reconstructed)


@run_multiple_genomes()
def test_json_ensembl_exon(genome):
    gene = genome.gene_by_id(TP53_gene_id)
    transcript = gene.transcripts[0]
    exon = transcript.exons[0]
    exon_reconstructed = Exon.from_json(exon.to_json())
    eq_(exon, exon_reconstructed)


@run_multiple_genomes()
def test_pickle_ensembl_genome(genome):
    genome_pickled = pickle.dumps(genome)
    genome_reconstructed = pickle.loads(genome_pickled)
    eq_(genome, genome_reconstructed)

    # This Genome happens to be an EnsemblRelease; test that too.
    eq_(genome.release, genome_reconstructed.release)
    eq_(genome.species, genome_reconstructed.species)


@run_multiple_genomes()
def test_ensembl_genome_to_dict(genome):
    genome_dict = genome.to_dict()
    genome_reconstructed = genome.__class__.from_dict(genome_dict)
    eq_(genome, genome_reconstructed)


@run_multiple_genomes()
def test_ensembl_genome_to_json(genome):
    genome_json = genome.to_json()
    genome_class = genome.__class__
    genome_reconstructed = genome_class.from_json(genome_json)
    eq_(genome, genome_reconstructed)


def test_custom_genome_to_json():
    setup_init_custom_mouse_genome()
    json = custom_mouse_genome_grcm38_subset.to_json()
    reconstructed = Genome.from_json(json)
    eq_(custom_mouse_genome_grcm38_subset, reconstructed)


def test_custom_genome_to_dict():
    setup_init_custom_mouse_genome()
    reconstructed = Genome.from_dict(custom_mouse_genome_grcm38_subset.to_dict())
    eq_(custom_mouse_genome_grcm38_subset, reconstructed)


def test_species_to_dict():
    eq_(human, Species.from_dict(human.to_dict()))


def test_species_to_json():
    eq_(human, Species.from_json(human.to_json()))


def test_species_to_pickle():
    eq_(human, pickle.loads(pickle.dumps(human)))


@run_multiple_genomes()
def test_unique_memory_address_of_unpickled_genomes(genome):
    unpickled = pickle.loads(pickle.dumps(genome))
    assert (
        genome is unpickled
    ), "Expected same object for %s but got two different instances" % (unpickled,)

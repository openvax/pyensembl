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

import json
import pickle

import pytest

from pyensembl import EnsemblRelease, Genome, Transcript, Gene, Exon
from pyensembl.species import Species, human, mouse

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


def test_genome_species_and_assembly_accession_roundtrip():
    genome = Genome(
        reference_name="T2T-CHM13v2.0",
        annotation_name="Ensembl",
        annotation_version="2022_07",
        species="human",
        assembly_accession="GCA_009914755.4",
    )
    assert genome.species == human
    assert genome.assembly_accession == "GCA_009914755.4"
    assert genome.taxon_id == 9606
    assert genome.species_taxon_id == 9606
    assert genome.nuclear_translation_table == 1
    assert genome.mitochondrial_translation_table == 2
    assert genome.genetic_code("MT").table_id == 2

    state = genome.to_dict()
    assert state["species"] == "homo_sapiens"
    assert state["assembly_accession"] == "GCA_009914755.4"

    for restored in (
        Genome.from_dict(state),
        Genome.from_json(genome.to_json()),
        pickle.loads(pickle.dumps(genome)),
    ):
        assert restored == genome
        assert hash(restored) == hash(genome)
        assert restored.species == human
        assert restored.assembly_accession == "GCA_009914755.4"
        assert restored.taxon_id == 9606
        assert restored.mitochondrial_translation_table == 2

    # Species object input matches string input and differs from another species/accession.
    same = Genome(
        reference_name="T2T-CHM13v2.0",
        annotation_name="Ensembl",
        annotation_version="2022_07",
        species=human,
        assembly_accession="GCA_009914755.4",
    )
    assert same == genome
    different_accession = Genome(
        reference_name="T2T-CHM13v2.0",
        annotation_name="Ensembl",
        annotation_version="2022_07",
        species=human,
        assembly_accession="GCF_009914755.1",
    )
    assert different_accession != genome
    different_species = Genome(
        reference_name="T2T-CHM13v2.0",
        annotation_name="Ensembl",
        annotation_version="2022_07",
        species=mouse,
        assembly_accession="GCA_009914755.4",
    )
    assert different_species != genome
    assert different_species.taxon_id == 10090


def test_genome_species_inference_and_legacy_deserialization():
    inferred = Genome(reference_name="GRCh38", annotation_name="custom")
    assert inferred.species == human
    assert inferred.assembly_accession is None
    assert inferred.taxon_id == 9606

    # Legacy serialized payloads without species or assembly_accession still load cleanly.
    legacy_state = inferred.to_dict()
    legacy_state.pop("species")
    legacy_state.pop("assembly_accession")
    restored = Genome.from_dict(legacy_state)
    assert restored == inferred
    assert restored.species == human
    assert restored.assembly_accession is None

    # Unknown reference_name stays None when species is omitted.
    unregistered = Genome(reference_name="CustomAsm1", annotation_name="custom")
    assert unregistered.species is None
    assert unregistered.assembly_accession is None
    assert unregistered.taxon_id is None
    legacy_payload = json.loads(unregistered.to_json())
    legacy_payload.pop("species")
    legacy_payload.pop("assembly_accession")
    restored_unregistered = Genome.from_json(json.dumps(legacy_payload))
    assert restored_unregistered == unregistered
    assert restored_unregistered.species is None
    assert restored_unregistered.assembly_accession is None


def test_genome_species_and_accession_validation_and_subclasses():
    with pytest.raises(ValueError, match="Species not found"):
        Genome(reference_name="CustomAsm1", annotation_name="custom", species="not_a_species")
    with pytest.raises(ValueError, match="disagree"):
        Genome(
            reference_name="T2T-CHM13v2.0",
            annotation_name="custom",
            species="human",
            reference_metadata={"taxon_id": 10090},
        )
    with pytest.raises(TypeError, match="assembly_accession"):
        Genome(reference_name="CustomAsm1", annotation_name="custom", assembly_accession=123)
    with pytest.raises(ValueError, match="assembly_accession"):
        Genome(reference_name="CustomAsm1", annotation_name="custom", assembly_accession="   ")

    # Legacy subclasses (e.g. orfdata.PinnedEnsemblGenome) may set species before
    # super().__init__ and assembly_accession after super().__init__.
    class LegacySubclass(Genome):
        def __init__(self, *args, species=None, assembly_accession=None, **kwargs):
            self.species = human if species == "human" else species
            super().__init__(*args, **kwargs)
            self.assembly_accession = assembly_accession

    sub = LegacySubclass(
        reference_name="T2T-CHM13v2.0",
        annotation_name="Ensembl",
        species="human",
        assembly_accession="GCA_009914755.4",
    )
    assert sub.species == human
    assert sub.assembly_accession == "GCA_009914755.4"
    assert sub.taxon_id == 9606

    # EnsemblRelease.with_mirror preserves species and assembly_accession on base Genome.
    release = EnsemblRelease(75, species="human")
    mirrored = release.with_mirror("https://mirror.example/release/")
    assert type(mirrored) is Genome
    assert mirrored.species == human
    assert mirrored.assembly_accession == release.assembly_accession
    assert Genome.from_json(mirrored.to_json()) == mirrored

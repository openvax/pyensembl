"""Genetic-code assignments from captured NCBI/Ensembl evidence, offline."""

from copy import deepcopy
import importlib
import json
from pathlib import Path
import pickle
import runpy

import pytest
import requests

from pyensembl import EnsemblRelease, Genome
from pyensembl import ensembl_metadata
from pyensembl.species import Species


codes = importlib.import_module("pyensembl.genetic_code")
DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("species,release,taxon,parent,nuclear,mitochondrial,plastid", [
    ("human", 75, 9606, 9606, 1, 2, None),
    ("mouse", 81, 10090, 10090, 1, 2, None),
    ("dog", 104, 9615, 9612, 1, 2, None),
    ("dog_basenji", 116, 9615, 9612, 1, 2, None),
    ("yeast", 116, 559292, 4932, 1, 3, None),
    ("fly", 116, 7227, 7227, 1, 5, None),
    ("caenorhabditis_elegans", 116, 6239, 6239, 1, 5, None),
    ("candida_albicans", 63, 237561, 5476, 12, 4, None),
    ("arabidopsis_thaliana", 63, 3702, 3702, 1, 1, 11),
    ("rice", 63, 39947, 4530, 1, 1, 11),
    ("tomato", 63, 4081, 4081, 1, 1, 11),
])
def test_reference_assignments(species, release, taxon, parent, nuclear, mitochondrial, plastid):
    genome = EnsemblRelease(release, species=species)
    assert (genome.taxon_id, genome.species_taxon_id) == (taxon, parent)
    assert genome.nuclear_translation_table == nuclear
    assert genome.mitochondrial_translation_table == mitochondrial
    assert genome.plastid_translation_table == plastid
    for contig in (None, "unplaced_42", "GL000220.1", "chr1_KI270706v1_random"):
        assignment = genome.genetic_code(contig)
        assert (assignment.table_id, assignment.compartment, assignment.source) == (nuclear, "nuclear", "taxonomy")
        assert assignment.taxon_id == taxon
        assert assignment.reason is None
    for compartment, expected in zip(codes.COMPARTMENTS, (nuclear, mitochondrial, plastid)):
        assignment = genome.genetic_code(compartment=compartment)
        assert assignment.table_id == expected
        assert assignment.taxon_id == taxon
        assert assignment.source == ("taxonomy" if expected else "unresolved")
        assert bool(assignment.reason) == (expected is None)
    record = genome.reference_metadata
    assert "ncbi.nlm.nih.gov" in record["taxonomy_source"]["url"]
    assert len(record["taxonomy_source"]["sha256"]) == 64
    assert "species_Ensembl" in record["taxon_source"]


def test_snapshot_rebuilds_from_captured_primary_evidence():
    # Exercises the parser against verbatim independent upstream records;
    # does not merely compare a declarative list with another hardcoded list.
    builder = runpy.run_path(str(Path(__file__).parents[1] / "scripts/build_genetic_code_catalog.py"))
    assert builder["build"](DATA / "genetic_codes") == codes.catalogue()
    for species in Species._latin_names_to_species.values():
        release = max(end for _, end in species.reference_assemblies.values())
        genome = EnsemblRelease(release, species=species)
        assert genome.taxon_id and genome.species_taxon_id
        assert genome.nuclear_translation_table is not None
        assert genome.genetic_code("scaffold_42").table_id == genome.nuclear_translation_table


@pytest.mark.parametrize("contig,compartment", [
    ("MT", "mitochondrial"), ("chrM", "mitochondrial"), ("mitochondrion_genome", "mitochondrial"),
    ("Pt", "plastid"), ("chloroplast", "plastid"), ("apicoplast", "plastid"),
    ("1", "nuclear"), ("chrX", "nuclear"), ("XVI", "nuclear"), ("2L", "nuclear"),
])
def test_documented_aliases(contig, compartment):
    assert EnsemblRelease(63, species="tomato").genetic_code(contig).compartment == compartment


def test_explicit_compartments_precede_aliases():
    genome = Genome("custom", "test", reference_metadata={
        "taxon_id": 7227, "contig_compartments": {"MT": "nuclear", "scaffold_42": "mitochondrial"},
    })
    assert genome.genetic_code("MT").table_id == 1
    assert genome.genetic_code("SCAFFOLD_42").table_id == 5
    assert genome.genetic_code("unplaced").table_id == 1
    assert genome.genetic_code("unplaced", compartment="nuclear").table_id == 1
    assert genome.genetic_code("unplaced", compartment="mitochondrial").table_id == 5
    assert genome.genetic_code(compartment="mitochondrial").table_id == 5
    with pytest.raises(ValueError, match="conflicts"):
        genome.genetic_code("MT", compartment="mitochondrial")
    with pytest.raises(ValueError, match="conflicts"):
        genome.genetic_code("chr1", compartment="mitochondrial")
    with pytest.raises(ValueError, match="Unknown cellular"):
        genome.genetic_code(compartment="cytoplasm")
    ambiguous = Genome("custom", "test", reference_metadata={"contig_compartments": {
        "Scaffold": "nuclear", "scaffold": "mitochondrial",
    }})
    with pytest.raises(ValueError, match="Conflicting reference compartments"):
        ambiguous.genetic_code("SCAFFOLD")


def test_unknowns_are_not_silently_standard_or_vertebrate():
    unknown = Genome("custom", "test", reference_metadata={"taxon_id": 999999999})
    assert unknown.taxon_id == 999999999
    assert unknown.species_taxon_id is None
    assert unknown.genetic_code("1").table_id is None
    assert unknown.genetic_code("MT").table_id is None
    assert unknown.genetic_code("scaffold_7").table_id is None
    assert unknown.genetic_code().table_id is None
    assert Genome("custom", "test").genetic_code("chr1").table_id is None
    # A display label on an unregistered assembly is not evidence of taxonomy.
    unknown = EnsemblRelease("2026_04", assembly="GCA_999999999.1", reference_name="GRCh38")
    assert unknown.taxon_id is None
    assert unknown.genetic_code("MT").table_id is None
    assert Genome("GRCh38", "legacy").genetic_code("MT").table_id == 2


def test_explicit_reference_codes_provenance_and_copies():
    options = {"taxon_id": "7227", "translation_tables": {"mitochondrial": 4},
               "source_url": "https://example.org/reference", "contig_compartments": {"scaffold": "mitochondrial"}}
    genome = Genome("custom", "test", reference_metadata=options)
    options["translation_tables"]["mitochondrial"] = 2
    record = genome.reference_metadata
    record["translation_tables"]["mitochondrial"] = 2
    genome.contig_compartments["scaffold"] = "nuclear"
    result = genome.genetic_code("scaffold")
    assert (result.table_id, result.source, result.source_url) == (4, "reference", "https://example.org/reference")
    assert genome.nuclear_translation_table == 1
    with pytest.raises(AttributeError):
        result.table_id = 2


@pytest.mark.parametrize("metadata", [
    {"taxon_id": True}, {"taxon_id": 0}, {"taxon_id": "human"}, {"taxon_id": 9606.0},
    {"species_taxon_id": False}, {"translation_tables": {"nuclear": 999}},
    {"translation_tables": {"nuclear": True}}, {"translation_tables": {"nuclear": "1"}},
    {"translation_tables": {"cytoplasm": 1}}, {"translation_tables": []},
    {"contig_compartments": {"MT": "cytoplasm"}}, {"contig_compartments": []},
    {"source_url": 7}, {"misspelled": 1}, {"schema_version": 2}, {"schema_version": 1},
])
def test_invalid_reference_metadata_rejected(metadata):
    with pytest.raises(ValueError):
        Genome("custom", "test", reference_metadata=metadata)


def test_taxon_parent_constraints_and_unrelated_conflicts():
    assert EnsemblRelease(116, species="yeast", reference_metadata={"taxon_id": 4932}).taxon_id == 4932
    assert EnsemblRelease(116, species="dog", reference_metadata={"species_taxon_id": 9612}).taxon_id == 9615
    for metadata in ({"taxon_id": 7227}, {"species_taxon_id": 7227}):
        with pytest.raises(ValueError, match="disagree"):
            EnsemblRelease(75, reference_metadata=metadata)
    with pytest.raises(ValueError, match="disagree"):
        Genome("custom", "test", reference_metadata={"taxon_id": 559292, "species_taxon_id": 9606})
    with pytest.raises(ValueError, match="verified species ancestor"):
        EnsemblRelease(75, reference_metadata={"species_taxon_id": 9605})


@pytest.mark.parametrize("contig", [None, "1", "scaffold_42"])
def test_explicit_unknown_table_remains_unknown(contig):
    genome = Genome("GRCh38", "test", reference_metadata={"translation_tables": {"nuclear": None}})
    assert genome.nuclear_translation_table is None
    assignment = genome.genetic_code(contig)
    assert assignment.table_id is None and assignment.source == "reference"
    assert assignment.reason
    assert genome.genetic_code("MT").table_id == 2


def test_missing_organelle_table_does_not_use_nuclear_default():
    genome = Genome("GRCh38", "test", reference_metadata={
        "translation_tables": {"mitochondrial": None},
        "contig_compartments": {"scaffold_42": "mitochondrial"},
    })
    assert genome.genetic_code("unplaced_1").table_id == 1
    for contig in ("MT", "scaffold_42", "Pt"):
        assignment = genome.genetic_code(contig)
        assert assignment.table_id is None
        assert assignment.reason


def test_pinned_json_pickle_mirror_and_older_objects(monkeypatch):
    genome = EnsemblRelease("2020_10", species="dog", reference_name="study")
    encoded, payload, expected = genome.to_json(), pickle.dumps(genome), genome.reference_metadata
    old = genome.to_dict()
    old.pop("reference_metadata")
    assert EnsemblRelease.from_dict(old).taxon_id == 9615
    changed = deepcopy(codes.catalogue())
    changed["taxa"]["9615"]["translation_tables"]["mitochondrial"] = 4
    changed["taxa"]["9615"]["translation_tables"]["nuclear"] = 12
    monkeypatch.setattr(codes, "catalogue", lambda: changed)
    for restored in [EnsemblRelease.from_json(encoded), pickle.loads(payload), genome.with_mirror("https://example.org/data")]:
        assert restored.reference_metadata == expected
        assert restored.genetic_code("MT").table_id == 2
        assert restored.genetic_code("unplaced_42").table_id == 1
        assert restored.genetic_code().table_id == 1
        assert type(restored).from_json(restored.to_json()).reference_metadata == expected
    assert EnsemblRelease("2020_10", species="dog").genetic_code("MT").table_id == 4
    assert EnsemblRelease("2020_10", species="dog").genetic_code("unplaced_42").table_id == 12


def test_overrides_distinguish_objects_but_preserve_storage_identity(tmp_path):
    options = dict(cache_directory_path=str(tmp_path))
    first = EnsemblRelease.cached(75, **options)
    second = EnsemblRelease.cached(75, **options, reference_metadata={"translation_tables": {"mitochondrial": 4}})
    assert first is not second and first != second
    assert len({first, second}) == 2
    assert first.gtf_url == second.gtf_url
    assert first.required_local_files() == second.required_local_files()
    assert not list(tmp_path.iterdir())


def test_api_taxa_are_retained_and_validated_offline_after_selection(monkeypatch):
    details = json.loads((DATA / "assembly_metadata/dog-details.json").read_text())
    links = json.loads((DATA / "assembly_metadata/dog-ftplinks.json").read_text())
    monkeypatch.setattr(ensembl_metadata, "_get_json", lambda path: details if path.endswith("details") else links)
    genome = EnsemblRelease.from_ensembl(details["genome_id"])
    assert (genome.taxon_id, genome.species_taxon_id) == (9615, 9612)
    assert genome.reference_metadata["taxon_source"].endswith("/details")
    encoded = genome.to_json()
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("network"))
    with monkeypatch.context() as offline:
        offline.setattr(ensembl_metadata, "_get_json", lambda *a: pytest.fail("metadata request while loading"))
        assert EnsemblRelease.from_json(encoded).genetic_code("MT").table_id == 2
    details["taxonomy_id"] = "9606"
    with pytest.raises(ValueError, match="disagree"):
        EnsemblRelease.from_ensembl(details["genome_id"])
    details["taxonomy_id"] = True
    with pytest.raises(ValueError, match="positive NCBI"):
        EnsemblRelease.from_ensembl(details["genome_id"])

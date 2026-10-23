"""Selection policy, scoped evidence, old-cache reuse and concrete snapshots."""

from copy import deepcopy
import json
from pathlib import Path
import pickle
import sys

import pytest

from pyensembl import EnsemblRelease, available_releases, genome_for_reference_name, shell
from pyensembl import release_catalogue
from pyensembl.dated_releases import _listing_location
from .test_legacy_cache_compatibility import snapshot
from . import test_release_cache_combinations as cache_matrix

sources = cache_matrix.sources


@pytest.mark.parametrize("options, release, assembly", [
    ({}, 116, "GRCh38"), ({"assembly": "GRCh37"}, 75, "GRCh37"),
    ({"assembly": "GRCh37.p13"}, 75, "GRCh37"),
    ({"assembly_accession": "GCA_000001405.14"}, 75, "GRCh37"),
    ({"species": "dog"}, 116, "ROS_Cfam_1.0"),
    ({"species": "dog_basenji"}, 116, "Basenji_breed-1.1"),
    ({"species": "mouse"}, 116, "GRCm39"),
    ({"species": "fly"}, 116, "BDGP6.54"),
    ({"species": "yeast"}, 116, "R64-1-1"),
    ({"species": "arabidopsis_thaliana"}, 63, "TAIR10"),
    ({"species": "tomato"}, 63, "SL4.0"),
    ({"species": "human", "provider": "ensembl"}, "2026_04", "GRCh38"),
])
def test_latest_is_scoped_offline_and_uses_existing_concrete_paths(tmp_path, monkeypatch, options, release, assembly):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "absent"))
    monkeypatch.setattr(release_catalogue, "available_annotation_dates",
                        lambda *a, **kw: pytest.fail("offline selection used network"))
    selected = EnsemblRelease("latest", **options)
    concrete = EnsemblRelease(release, **options)
    assert selected.release == release
    assert selected.assembly == assembly
    assert selected.gtf_url == concrete.gtf_url
    assert selected.download_cache.cache_directory_path == concrete.download_cache.cache_directory_path
    assert selected.to_dict()["_resolved"]["selection"]["policy"] == "latest"
    assert not (tmp_path / "absent").exists()
    assert sum(row["latest"] for row in available_releases(**options)) == 1


def test_snapshot_does_not_resolve_latest_again(monkeypatch):
    genome = EnsemblRelease("latest", species="dog")
    saved = genome.to_json()
    monkeypatch.setattr(release_catalogue, "resolve_latest", lambda **kw: pytest.fail("resolved again"))
    for restored in (EnsemblRelease.from_json(saved), pickle.loads(pickle.dumps(genome))):
        assert restored.release == 116
        assert restored.assembly_accession == genome.assembly_accession
        assert restored.to_dict()["_resolved"] == genome.to_dict()["_resolved"]


def test_discovered_dates_are_visible_but_cannot_displace_frozen_release(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    url, path = _listing_location("GCA_000001405.29", "ensembl")
    path = Path(path)
    path.parent.mkdir(parents=True)
    path.write_text('<a href="2099_01/">2099_01</a>')
    before = snapshot(tmp_path)
    rows = available_releases()
    future, = [row for row in rows if row["release"] == "2099_01"]
    assert not future["latest"] and not future["comparisons"]
    assert EnsemblRelease("latest").release == 116
    assert snapshot(tmp_path) == before
    refreshed = []
    monkeypatch.setattr(release_catalogue, "available_annotation_dates",
                        lambda *a, **kw: refreshed.append((a, kw)) or ["2099_02"])
    assert "2099_02" in {r["release"] for r in available_releases(refresh=True)}
    assert refreshed and all(kw == {"refresh": True} for a, kw in refreshed)


def test_future_reviewed_preference_is_concrete_and_serializable(monkeypatch):
    evidence = deepcopy(release_catalogue.comparison_evidence())
    evidence["preferred"] = [dict(release="2026_04", assembly_accession="GCA_000001405.29",
                                  provider="ensembl", coverage="including-alt", reason="Reviewed newer annotation")]
    monkeypatch.setattr(release_catalogue, "comparison_evidence", lambda: evidence)
    genome = EnsemblRelease("latest")
    assert genome.release == "2026_04"
    assert EnsemblRelease.from_json(genome.to_json()).release == "2026_04"
    assert genome.gtf_url.endswith("genes-including_alt.gtf.gz")


def test_evidence_distinguishes_identifiers_content_and_coverage():
    human = {r["release"]: r for r in available_releases()}
    assert human["2026_04"]["comparisons"][0]["relationship"] == "same-gtf-records"
    subset = human["2024_11"]["comparisons"][0]
    assert subset["relationship"] == "transcript-subset"
    assert subset["transcripts"]["only_dated"] == 0
    assert subset["transcripts"]["only_numbered"] > 0
    dog = [r for r in available_releases("dog") if r["kind"] == "dated"]
    assert dog[0]["comparisons"][0]["relationship"] == "same-identifiers"
    primary = available_releases(assembly_accession="GCA_000001405.29")
    assert all(not r["comparisons"] for r in primary if r["kind"] == "dated")
    with pytest.warns(UserWarning, match="transcript subset"):
        assert EnsemblRelease("2024_11").release == "2024_11"


@pytest.mark.parametrize("selection, reference", [
    ({}, "GRCh38"),
    ({"assembly_accession": "GCA_000001405.29", "include_alt": True}, "GRCh38"),
    ({"assembly": "GRCh37"}, "GRCh37"),
    ({"species": "dog"}, "ROS_Cfam_1.0"),
])
def test_reference_lookup_reuses_installed_dated_data_read_only(sources, monkeypatch, selection, reference):
    date = "2013_09" if reference == "GRCh37" else "2020_10" if reference == "ROS_Cfam_1.0" else "2026_04"
    genome = EnsemblRelease(date, **selection)
    expected = sources.register(genome, "INSTALLED", alt=genome.include_alt)
    genome.download()
    genome.index()
    before = snapshot(sources.root)
    sources.offline = True
    selected = genome_for_reference_name(reference)
    assert selected.release == date
    assert selected.gene_by_id(expected["gene"]).name == "INSTALLED"
    assert selected.download_cache.cache_directory_path == genome.download_cache.cache_directory_path
    assert any(g.release == date for g in shell.collect_all_installed_ensembl_releases())
    assert snapshot(sources.root) == before
    assert isinstance(genome_for_reference_name(reference, False).release, int)


def test_numbered_installation_keeps_precedence(sources):
    for release in (115, "2026_04"):
        genome = EnsemblRelease(release)
        sources.register(genome, str(release), alt=True)
        genome.download()
        genome.index()
    sources.offline = True
    assert genome_for_reference_name("GRCh38").release == 115


def test_cli_available_json_and_latest_reference_selection(monkeypatch, capsys):
    monkeypatch.setattr(shell, "configure_logging", lambda **kw: None)
    monkeypatch.setattr(sys, "argv", ["pyensembl", "available", "--species", "dog", "--json"])
    shell.run()
    rows = json.loads(capsys.readouterr().out)
    assert any(r["release"] == "2020_10" for r in rows)
    assert [r["release"] for r in rows if r["latest"]] == [116]
    args = shell.parser.parse_args(["install", "--reference-name", "GRCh37", "--release", "latest"])
    genome, = shell.collect_selected_genomes(args)
    assert genome.release == 75
    monkeypatch.setattr(sys, "argv", ["pyensembl", "delete-all-files", "--release", "latest"])
    with pytest.raises(SystemExit) as error:
        shell.run()
    assert error.value.code == 2
    assert "concrete" in capsys.readouterr().err


def test_custom_dated_scope_uses_cached_dates_without_assuming_human(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    accession = "GCF_999999999.1"
    url, path = _listing_location(accession, "community")
    path = Path(path)
    path.parent.mkdir(parents=True)
    path.write_text('<a href="2025_01/">old</a><a href="2026_02/">new</a>')
    for selection in ({"assembly_accession": accession}, {"assembly": accession}):
        genome = EnsemblRelease("latest", provider="community", **selection)
        assert genome.release == "2026_02"
        assert genome.species is None and genome.assembly_accession == accession
        assert genome.provider == "community"
        assert EnsemblRelease.from_json(genome.to_json()).release == "2026_02"


def test_api_subset_label_is_preserved_without_replacing_release(monkeypatch):
    from pyensembl import ensembl_metadata
    fixture = Path(__file__).parent / "data/assembly_metadata"
    details = json.loads((fixture / "human-details.json").read_text())
    links = json.loads((fixture / "human-ftplinks.json").read_text())
    details["annotation_version"] = "51_subset"
    monkeypatch.setattr(ensembl_metadata, "_get_json",
                        lambda path: details if path.endswith("details") else links)
    with pytest.warns(UserWarning, match="transcript subset"):
        genome = EnsemblRelease.from_ensembl(details["genome_id"])
    assert genome.release == genome.annotation_version == "2026_04"
    assert genome.annotation_label == "51_subset" and genome.annotation_is_subset
    monkeypatch.setattr(ensembl_metadata, "_get_json", lambda *a: pytest.fail("network"))
    with pytest.warns(UserWarning, match="transcript subset"):
        restored = EnsemblRelease.from_json(genome.to_json())
    assert restored.annotation_label == "51_subset" and restored.annotation_is_subset

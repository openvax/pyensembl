"""Captured responses from the new API, not the frozen legacy REST service."""

from copy import deepcopy
import json
from pathlib import Path
import pickle

import pytest

from pyensembl import EnsemblRelease
from pyensembl import ensembl_metadata as api


DATA = Path(__file__).parent / "data" / "assembly_metadata"


@pytest.fixture
def responses(monkeypatch, tmp_path):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    documents = {}
    for species in ["human", "dog"]:
        explain = json.loads((DATA / ("%s-explain.json" % species)).read_text())
        identifier = explain["genome_id"]
        accession = explain["assembly"]["accession_id"]
        documents["/genome/%s/explain" % accession] = explain
        for kind in ["details", "ftplinks"]:
            documents["/genome/%s/%s" % (identifier, kind)] = json.loads(
                (DATA / ("%s-%s.json" % (species, kind))).read_text()
            )
    monkeypatch.setattr(api, "_get_json", lambda path: deepcopy(documents[path]))
    return documents


@pytest.mark.parametrize("identifier, name, accession, date", [
    ("GCA_000001405.29", "GRCh38.p14", "GCA_000001405.29", "2026_04"),
    ("GCA_014441545.1", "ROS_Cfam_1.0", "GCA_014441545.1", "2020_10"),
    ("a6edaf0c-daed-4f0f-8de3-f732502973b2", "ROS_Cfam_1.0", "GCA_014441545.1", "2020_10"),
])
def test_api_discovery_pins_dataset_and_serializes_offline(
    responses, monkeypatch, identifier, name, accession, date, tmp_path,
):
    genome = EnsemblRelease.from_ensembl(identifier)
    assert genome.assembly_name == name
    assert genome.assembly_accession == accession
    assert genome.release == date
    assert genome.website_release == "2026-07"
    assert genome.genome_id
    assert not (tmp_path / "cache").exists()
    encoded, pickled = genome.to_json(), pickle.dumps(genome)
    monkeypatch.setattr(api, "_get_json", lambda path: pytest.fail("unexpected network"))
    for restored in [EnsemblRelease.from_json(encoded), pickle.loads(pickled)]:
        assert restored.assembly_name == name
        assert restored.release == date
        assert restored.genome_id == genome.genome_id
        assert restored.website_release == genome.website_release
        assert restored.download_cache.cache_directory_path == genome.download_cache.cache_directory_path


@pytest.mark.parametrize("options", [
    {"species": "mouse"}, {"assembly": "GRCh37"}, {"provider": "community"},
])
def test_api_conflicting_constraints_fail(responses, options):
    with pytest.raises(ValueError):
        EnsemblRelease.from_ensembl("GCA_000001405.29", **options)


@pytest.mark.parametrize("change", [
    "accession", "suppressed", "wrong_uuid", "no_name", "date", "ftp_host",
    "ftp_accession", "ftp_missing", "ftp_multiple", "dna_conflict",
    "malformed_assembly", "malformed_links", "malformed_date",
])
def test_inconsistent_api_metadata_is_rejected(responses, change):
    identifier = responses["/genome/GCA_000001405.29/explain"]["genome_id"]
    details = responses["/genome/%s/details" % identifier]
    links = responses["/genome/%s/ftplinks" % identifier]
    if change == "accession":
        details["assembly"]["accession_id"] = "GCA_000001405.14"
    elif change == "suppressed":
        details["is_suppressed"] = True
    elif change == "wrong_uuid":
        details["genome_id"] = "a6edaf0c-daed-4f0f-8de3-f732502973b2"
    elif change == "no_name":
        details["assembly"].pop("name")
    elif change == "date":
        details["annotation_date"] = "2025-12"
    elif change == "malformed_date":
        details["annotation_date"] = 202604
    elif change == "malformed_assembly":
        details["assembly"] = None
    elif change == "malformed_links":
        links.append(None)
    elif change == "ftp_missing":
        links[:] = [r for r in links if r["dataset"] != "genebuild"]
    elif change == "ftp_multiple":
        links.append(next(r for r in links if r["dataset"] == "genebuild"))
    else:
        target = next(r for r in links if r["dataset"] ==
                      ("assembly" if change == "dna_conflict" else "genebuild"))
        if change == "ftp_host":
            target["url"] = target["url"].replace("ftp.ebi.ac.uk", "untrusted.example")
        elif change == "ftp_accession":
            target["url"] = target["url"].replace("/405/29/", "/405/14/")
        else:
            target["url"] = target["url"].replace("2026_04", "2025_12")
    with pytest.raises(ValueError):
        EnsemblRelease.from_ensembl("GCA_000001405.29")


def test_api_outage_is_not_an_assembly_guess(monkeypatch):
    def unavailable(path):
        raise OSError("metadata unavailable")
    monkeypatch.setattr(api, "_get_json", unavailable)
    with pytest.raises(OSError, match="metadata unavailable"):
        EnsemblRelease.from_ensembl("GCA_014441545.1")
    assert EnsemblRelease("2020_10", species="dog").assembly_name == "ROS_Cfam_1.0"


def test_actual_http_client_uses_new_api(monkeypatch):
    calls = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def read(self, size):
            return b'{"ok":true}'

    def open_request(request, timeout):
        calls.append((request.full_url, timeout))
        return Response()

    monkeypatch.setattr(api, "urlopen", open_request)
    assert api._get_json("/genome/test/details") == {"ok": True}
    assert calls == [("https://www.ensembl.org/api/metadata/genome/test/details", 30)]

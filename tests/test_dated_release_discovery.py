"""Discovering dated releases on the new Ensembl platform (#447 phase 2)."""

from datetime import datetime, timedelta, timezone
import json
import logging
import os
from pathlib import Path
import re

import pytest
import requests

from pyensembl import (
    EnsemblAnnotation,
    EnsemblRelease,
    UnpublishedDateError,
    available_annotation_dates,
    available_dated_releases,
    shell,
)
from pyensembl import dated_releases
from pyensembl.dated_releases import (
    explain_unpublished_date,
    fetch_all_dated_releases,
    parse_dated_release_listing,
)
from pyensembl.ensembl_url_templates import make_dated_releases_directory
from pyensembl.genome import Genome
from pyensembl.genome_fasta import GenomeFasta
from pyensembl.species import Species, find_species_by_name

PLATFORM = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"
HUMAN_LISTING = PLATFORM + "/GCA/000/001/405/29/ensembl/"


def listing(*dates):
    links = "".join('<a href="%s/">%s/</a>' % (date, date) for date in dates)
    return '<a href="../">Parent</a><a href="README">README</a>%s' % links


def http_error(status, headers=None):
    response = requests.Response()
    response.status_code = status
    response.headers.update(headers or {})
    return requests.HTTPError("%d" % status, response=response)


class Server:
    """Stands in for requests.get: serves listings by URL, records requests.

    A listing may also be an exception, or a list of responses served in turn.
    """

    def __init__(self, **listings):
        self.listings = listings
        self.requests = []
        self.offline = False

    def get(self, url, timeout=None):
        self.requests.append(url)
        if self.offline:
            raise requests.ConnectionError("offline")
        served = self.listings.get(url, http_error(404))
        if isinstance(served, list):
            served = served.pop(0) if len(served) > 1 else served[0]
        if isinstance(served, Exception):
            raise served
        response = requests.Response()
        response.status_code = 200
        response._content = served.encode()
        return response


@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    fake = Server(**{HUMAN_LISTING: listing("2023_03", "2026_04")})
    monkeypatch.setattr(requests, "get", fake.get)
    monkeypatch.setattr(dated_releases.time, "sleep", lambda seconds: None)
    return fake


def dated_species():
    return [
        s for s in map(find_species_by_name, Species.all_registered_latin_names())
        if s.dated_releases is not None
    ]


def missing_file(*args, **kwargs):
    raise http_error(404)


# Listings


def test_listing_keeps_only_annotation_dates():
    html = listing("2026_04", "2023_03", "2026_13", "latest") + '<a href="2023_03/">again</a>'
    assert parse_dated_release_listing(html) == ["2023_03", "2026_04"]


def test_dates_are_fetched_once_then_read_offline(server):
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    server.offline = True
    assert available_dated_releases("homo_sapiens") == ["2023_03", "2026_04"]
    assert server.requests == [HUMAN_LISTING]
    with pytest.raises(OSError):
        available_dated_releases("human", refresh=True)


def test_refresh_finds_new_dates(server):
    available_dated_releases("human")
    server.listings[HUMAN_LISTING] = listing("2023_03", "2026_04", "2026_10")
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    assert available_dated_releases("human", refresh=True)[-1] == "2026_10"
    assert available_dated_releases("human")[-1] == "2026_10"


def test_dates_by_accession_and_provider(server):
    grch37 = PLATFORM + "/GCA/000/001/405/14/ensembl/"
    server.listings[grch37] = listing("2013_09")
    assert available_annotation_dates("GCA_000001405.14") == ["2013_09"]
    assert available_annotation_dates("GCA_000001405.29", "ensembl") == ["2023_03", "2026_04"]


def test_unusable_cache_is_fetched_again(server, tmp_path):
    path = tmp_path / "cache/pyensembl/dated_releases/GCA_000001405.29_ensembl.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"url": "https://mirror/", "dates": ["1999_01"]}))
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    assert json.loads(path.read_text())["url"] == HUMAN_LISTING


def test_species_without_dated_releases():
    with pytest.raises(ValueError, match="No dated Ensembl releases"):
        available_dated_releases("toxoplasma_gondii")


def test_page_without_dates_is_not_trusted(server):
    server.listings[HUMAN_LISTING] = "<html>Proxy error</html>"
    with pytest.raises(OSError, match="No annotation dates"):
        available_dated_releases("human")
    server.listings[HUMAN_LISTING] = listing("2026_04")
    assert available_dated_releases("human") == ["2026_04"]


@pytest.mark.parametrize("failure", [
    requests.ConnectionError("refused"),
    http_error(503, {"Retry-After": "1"}),
    http_error(429),
])
def test_transient_listing_failures_are_retried(server, failure):
    server.listings[HUMAN_LISTING] = [failure, failure, listing("2026_04")]
    assert available_dated_releases("human") == ["2026_04"]
    assert len(server.requests) == 3


@pytest.mark.parametrize("failure", [requests.Timeout("slow"), http_error(404)])
def test_timeouts_and_missing_listings_are_not_retried(server, failure):
    server.listings[HUMAN_LISTING] = failure
    with pytest.raises(OSError):
        available_dated_releases("human")
    assert len(server.requests) == 1


def test_unwritable_cache_returns_dates_and_warns_once(server, monkeypatch, caplog):
    def read_only(path, value):
        raise PermissionError(path)

    monkeypatch.setattr(dated_releases, "_write_json", read_only)
    monkeypatch.setattr(dated_releases, "_warned_unwritable_cache", False)
    with caplog.at_level(logging.DEBUG, logger="pyensembl"):
        assert available_dated_releases("human") == ["2023_03", "2026_04"]
        assert available_dated_releases("human") == ["2023_03", "2026_04"]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1 and "Couldn't cache annotation dates" in warnings[0].message


# Explaining downloads of unpublished dates


def test_missing_download_of_an_unpublished_date_lists_the_dates(server, monkeypatch):
    monkeypatch.setattr(Genome, "_get_gtf_path", missing_file)
    message = (
        "Ensembl publishes no homo_sapiens GRCh38 annotation (GCA_000001405.29, "
        "provider ensembl) dated 2013_09; available dates: 2023_03, 2026_04"
    )
    with pytest.raises(UnpublishedDateError, match=re.escape(message)) as raised:
        EnsemblRelease("2013_09").download()
    assert raised.value.__cause__.response.status_code == 404
    assert isinstance(raised.value, ValueError)


def test_missing_dna_of_an_unpublished_date_lists_the_dates(server, monkeypatch):
    monkeypatch.setattr(GenomeFasta, "prepare", missing_file)
    with pytest.raises(UnpublishedDateError, match="dated 2013_09"):
        EnsemblRelease("2013_09", genome_fasta=True).download_genome_fasta()


def test_accession_datasets_are_explained_too(server, monkeypatch):
    monkeypatch.setattr(Genome, "_get_gtf_path", missing_file)
    message = (
        "Ensembl publishes no annotation of GCA_000001405.29 (provider ensembl) "
        "dated 2013_09; available dates: 2023_03, 2026_04"
    )
    with pytest.raises(UnpublishedDateError, match=re.escape(message)):
        EnsemblAnnotation("GCA_000001405.29", "2013_09").download()


def test_missing_file_of_a_published_date_is_not_explained(server):
    with pytest.raises(requests.HTTPError):
        with explain_unpublished_date("GCA_000001405.29", "ensembl", "2026_04", PLATFORM):
            raise http_error(404)


def test_newly_published_date_is_found_by_refreshing(server):
    available_dated_releases("human")
    server.listings[HUMAN_LISTING] = listing("2023_03", "2026_04", "2026_10")
    with pytest.raises(requests.HTTPError):
        with explain_unpublished_date("GCA_000001405.29", "ensembl", "2026_10", PLATFORM):
            raise http_error(404)
    assert len(server.requests) == 2


@pytest.mark.parametrize("server_url", [PLATFORM + "/", "https://mirror.example"])
def test_only_the_official_server_is_explained(server, server_url):
    expected = UnpublishedDateError if server_url.rstrip("/") == PLATFORM else requests.HTTPError
    with pytest.raises(expected):
        with explain_unpublished_date("GCA_000001405.29", "ensembl", "2013_09", server_url):
            raise http_error(404)


def test_other_failures_and_unreachable_listings_keep_the_original_error(server):
    with pytest.raises(requests.HTTPError, match="503"):
        with explain_unpublished_date("GCA_000001405.29", "ensembl", "2013_09", PLATFORM):
            raise http_error(503)
    assert server.requests == []
    server.offline = True
    with pytest.raises(requests.HTTPError, match="404"):
        with explain_unpublished_date("GCA_000001405.29", "ensembl", "2013_09", PLATFORM):
            raise http_error(404)


def test_install_reports_an_unpublished_date(server, monkeypatch, capsys):
    # run() would attach a handler to this test's captured stderr for good.
    monkeypatch.setattr(shell, "configure_logging", lambda **kwargs: None)
    monkeypatch.setattr(Genome, "_get_gtf_path", missing_file)
    monkeypatch.setattr("sys.argv", ["pyensembl", "install", "--release", "2013_09"])
    with pytest.raises(SystemExit):
        shell.run()
    assert "dated 2013_09; available dates: 2023_03, 2026_04" in capsys.readouterr().err


# `pyensembl available`


def cache_listing(species, dates, age):
    url = make_dated_releases_directory(*species.dated_releases) + "/"
    path = Path(dated_releases._listing_location(*species.dated_releases)[1])
    path.parent.mkdir(parents=True, exist_ok=True)
    fetched = datetime.now(timezone.utc) - timedelta(seconds=age)
    path.write_text(json.dumps({"url": url, "fetched": fetched.isoformat(), "dates": dates}))


def test_available_uses_fresh_cached_dates_offline(server):
    for species in dated_species():
        cache_listing(species, ["2020_01"], age=60)
    server.offline = True
    dates = fetch_all_dated_releases()
    assert server.requests == []
    assert set(dates) == {s.latin_name for s in dated_species()}
    assert set(map(tuple, dates.values())) == {("2020_01",)}


def test_available_refreshes_stale_dates(server):
    for species in dated_species():
        cache_listing(species, ["2020_01"], age=2 * 24 * 60 * 60)
    dates = fetch_all_dated_releases()
    assert server.requests[0] == HUMAN_LISTING  # Human is the probe.
    assert dates["homo_sapiens"] == ["2023_03", "2026_04"]
    assert dates["mus_musculus"] == ["2020_01"]  # Not served here: cached dates.
    assert len(server.requests) == len(dated_species())


def test_available_offline_only_probes_human(server):
    cache_listing(find_species_by_name("human"), ["2020_01"], age=2 * 24 * 60 * 60)
    server.offline = True
    dates = fetch_all_dated_releases()
    # Connection errors are retried, but no other species is requested.
    assert set(server.requests) == {HUMAN_LISTING}
    assert dates["homo_sapiens"] == ["2020_01"]
    assert dates["mus_musculus"] is None


def test_available_table_shows_current_assembly_dates():
    table = shell.format_available_species(
        use_color=False,
        dated_releases={"homo_sapiens": ["2023_03", "2026_04"], "mus_musculus": None},
    )
    rows = {line.split()[0]: line for line in table.splitlines() if line.strip()}
    assert "Dated releases" in table.splitlines()[0]
    assert "2023_03, 2026_04" in rows["human"]  # GRCh38 is human's first row.
    grcm39 = next(line for line in table.splitlines() if "GRCm39" in line)
    assert " ? " in grcm39
    assert table.endswith("? = dates could not be fetched and none are cached")


def test_available_table_explains_species_without_dates():
    table = shell.format_available_species(use_color=False, dated_releases={})
    assert table.endswith("? = dates could not be fetched and none are cached")
    assert "?" not in shell.format_available_species(use_color=False)


# `pyensembl list`


def test_installed_dated_releases_are_listed(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    for release in ("2026_04", 116, "2023_03"):
        genome = EnsemblRelease(release)
        gtf = Path(genome._local_source_path(genome._gtf_path_or_url))
        gtf.parent.mkdir(parents=True, exist_ok=True)
        gtf.write_bytes(b"gtf")
    installed = [genome.release for genome in shell.collect_all_installed_ensembl_releases()]
    assert installed == [116, "2023_03", "2026_04"]
    table = shell.format_installed_genomes(use_color=False)
    releases = [line.split()[2] for line in table.splitlines()[1:]]
    assert releases == ["116", "2023_03", "2026_04"]


def test_other_files_in_a_release_directory_are_listed_as_custom(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    directory = Path(EnsemblRelease("2026_04").download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    (directory / "custom.gtf").write_text("")
    rows = shell.format_installed_genomes(use_color=False).splitlines()[1:]
    assert [row.split()[:3] for row in rows] == [["custom", "GRCh38", "ensembl2026_04"]]
    assert shell.collect_all_installed_ensembl_releases() == []


@pytest.mark.skipif(
    not os.environ.get("PYENSEMBL_NETWORK_TESTS"),
    reason="set PYENSEMBL_NETWORK_TESTS=1 to check the new Ensembl platform",
)
def test_dates_on_the_platform(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    assert {"2023_03", "2025_12", "2026_04"} <= set(available_dated_releases("human"))
    assert available_annotation_dates("GCA_000001405.14") == ["2013_09"]

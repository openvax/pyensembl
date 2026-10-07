"""Discovering dated releases on the new Ensembl platform (#447 phase 2)."""

import http.client
import io
import json
import os
from pathlib import Path
import re
import urllib.error

import pytest

from pyensembl import EnsemblRelease, available_dated_releases, shell
from pyensembl import dated_releases
from pyensembl.dated_releases import (
    fetch_all_dated_releases,
    parse_dated_release_listing,
    require_published_date,
)

HUMAN_LISTING = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/"


def listing(*dates):
    links = "".join('<a href="%s/">%s/</a>' % (date, date) for date in dates)
    return '<a href="../">Parent</a><a href="README">README</a>%s' % links


class Server:
    """Stands in for urllib: serves listings by URL and records requests."""

    def __init__(self, **listings):
        self.listings = listings
        self.requests = []
        self.offline = False

    def urlopen(self, url, timeout=None):
        self.requests.append(url)
        if self.offline:
            raise OSError("offline")
        if url not in self.listings:
            raise OSError("404 %s" % url)
        response = self.listings[url]
        if isinstance(response, Exception):
            raise response
        return io.BytesIO(response.encode())


@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    fake = Server(**{HUMAN_LISTING: listing("2023_03", "2026_04")})
    monkeypatch.setattr(dated_releases.urllib.request, "urlopen", fake.urlopen)
    return fake


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


def test_unusable_cache_is_fetched_again(server, tmp_path):
    path = tmp_path / "cache/pyensembl/dated_releases/GCA_000001405.29_ensembl.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"url": "https://mirror/", "dates": ["1999_01"]}))
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    assert json.loads(path.read_text())["url"] == HUMAN_LISTING


def test_species_without_dated_releases():
    with pytest.raises(ValueError, match="No dated Ensembl releases"):
        available_dated_releases("toxoplasma_gondii")


def test_unknown_date_fails_before_downloading(server):
    genome = EnsemblRelease("2013_09")
    message = (
        "Ensembl publishes no homo_sapiens annotation of GRCh38 "
        "(GCA_000001405.29, provider ensembl) dated 2013_09; "
        "available dates: 2023_03, 2026_04"
    )
    with pytest.raises(ValueError, match=re.escape(message)):
        genome.download()
    with pytest.raises(ValueError, match="dated 2013_09"):
        EnsemblRelease("2013_09", genome_fasta=True).download_genome_fasta()
    assert not os.path.exists(genome.download_cache.cache_directory_path)
    # The cached dates lacked 2013_09, so they were refreshed once per check.
    assert server.requests == [HUMAN_LISTING] * 3


def test_published_dates_are_checked_from_the_cache(server):
    available_dated_releases("human")
    server.offline = True
    require_published_date(EnsemblRelease("2026_04"))
    assert len(server.requests) == 1


def test_new_dates_are_found_by_refreshing(server):
    available_dated_releases("human")
    server.listings[HUMAN_LISTING] = listing("2023_03", "2026_04", "2026_10")
    require_published_date(EnsemblRelease("2026_10"))


def test_unreachable_listing_does_not_block_downloads(server):
    server.offline = True
    require_published_date(EnsemblRelease("2013_09"))


def test_mirrors_are_not_checked(server):
    require_published_date(EnsemblRelease("2013_09", server="https://mirror.example"))
    assert server.requests == []


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


def test_page_without_dates_is_not_trusted(server):
    server.listings[HUMAN_LISTING] = "<html>Proxy error</html>"
    with pytest.raises(OSError, match="No annotation dates"):
        available_dated_releases("human")
    require_published_date(EnsemblRelease("2026_04"))  # Fails open.
    server.listings[HUMAN_LISTING] = listing("2026_04")
    assert available_dated_releases("human") == ["2026_04"]


def test_truncated_response_counts_as_unreachable(server):
    server.listings[HUMAN_LISTING] = http.client.IncompleteRead(b"")
    with pytest.raises(OSError):
        available_dated_releases("human")
    require_published_date(EnsemblRelease("2013_09"))


def test_refused_connections_are_retried(server, monkeypatch):
    monkeypatch.setattr(dated_releases.time, "sleep", lambda seconds: None)
    refusals = [urllib.error.URLError(ConnectionRefusedError()), ConnectionResetError()]
    original = server.urlopen

    def flaky(url, timeout=None):
        if refusals:
            server.requests.append(url)
            raise refusals.pop(0)
        return original(url, timeout)

    monkeypatch.setattr(dated_releases.urllib.request, "urlopen", flaky)
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    assert len(server.requests) == 3


def test_timeouts_are_not_retried(server, monkeypatch):
    def timeout(url, timeout=None):
        server.requests.append(url)
        raise urllib.error.URLError(TimeoutError())

    monkeypatch.setattr(dated_releases.urllib.request, "urlopen", timeout)
    with pytest.raises(OSError):
        available_dated_releases("human")
    assert len(server.requests) == 1


def test_dates_are_returned_when_the_cache_is_not_writable(server, monkeypatch):
    def read_only(path, value):
        raise PermissionError(path)

    monkeypatch.setattr(dated_releases, "_write_json", read_only)
    assert available_dated_releases("human") == ["2023_03", "2026_04"]
    with pytest.raises(ValueError, match="dated 2013_09"):
        require_published_date(EnsemblRelease("2013_09"))


def test_local_dna_is_not_checked(server, tmp_path):
    fasta = tmp_path / "local.fa"
    fasta.write_text(">1\nACGT\n")
    EnsemblRelease("2013_09", genome_fasta=str(fasta)).download_genome_fasta()
    assert server.requests == []


def test_all_dates_fall_back_to_the_cache_after_one_failed_request(server):
    available_dated_releases("human")
    server.offline = True
    dates = fetch_all_dated_releases()
    assert len(server.requests) == 2  # The cached fetch, then one failed request.
    assert dates["homo_sapiens"] == ["2023_03", "2026_04"]
    assert dates["mus_musculus"] is None


def test_all_dates_are_fetched_with_per_species_fallback(server):
    dates = fetch_all_dated_releases()
    assert dates["homo_sapiens"] == ["2023_03", "2026_04"]
    assert dates["mus_musculus"] is None  # Not served here, and never cached.
    assert len(dates) == 42


def test_install_rejects_an_unpublished_date_before_installing(server, monkeypatch, capsys):
    # run() would attach a handler to this test's captured stderr for good.
    monkeypatch.setattr(shell, "configure_logging", lambda **kwargs: None)
    monkeypatch.setattr("sys.argv", ["pyensembl", "install", "--release", "2026_04", "2013_09"])
    with pytest.raises(SystemExit):
        shell.run()
    assert "dated 2013_09; available dates: 2023_03, 2026_04" in capsys.readouterr().err
    assert not os.path.exists(EnsemblRelease("2026_04").download_cache.cache_directory_path)


def test_custom_genome_named_like_a_dated_release_stays_listed(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    directory = Path(EnsemblRelease("2026_04").download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    (directory / "custom.gtf").write_text("")
    rows = shell.format_installed_genomes(use_color=False).splitlines()[1:]
    assert [row.split()[:3] for row in rows] == [["human", "GRCh38", "2026_04"]]
    assert shell.collect_all_installed_ensembl_releases() == []


@pytest.mark.skipif(
    not os.environ.get("PYENSEMBL_NETWORK_TESTS"),
    reason="set PYENSEMBL_NETWORK_TESTS=1 to check the new Ensembl platform",
)
def test_human_dates_on_the_platform(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    assert {"2023_03", "2025_12", "2026_04"} <= set(available_dated_releases("human"))

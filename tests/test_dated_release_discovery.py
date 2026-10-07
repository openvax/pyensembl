"""Discovering dated releases on the new Ensembl platform (#447 phase 2)."""

import io
import json
import os
from pathlib import Path
import re

import pytest

from pyensembl import EnsemblRelease, available_dated_releases, shell
from pyensembl import dated_releases
from pyensembl.dated_releases import parse_dated_release_listing, require_published_date

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
        return io.BytesIO(self.listings[url].encode())


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


@pytest.mark.skipif(
    not os.environ.get("PYENSEMBL_NETWORK_TESTS"),
    reason="set PYENSEMBL_NETWORK_TESTS=1 to check the new Ensembl platform",
)
def test_human_dates_on_the_platform(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    assert {"2023_03", "2025_12", "2026_04"} <= set(available_dated_releases("human"))

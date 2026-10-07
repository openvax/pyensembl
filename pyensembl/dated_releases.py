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

"""
Which annotation dates the new Ensembl platform publishes for a species.

Dates come from the directory listing of the species' current assembly and
provider, e.g. https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/,
which is exactly what dated releases download from. Ensembl's species.json
catalogue is 37 MB and omits some installable datasets, so it is not used.
"""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import http.client
import logging
import os
import re
import time
import urllib.error
import urllib.request

from .download_cache import cache_root
from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER, make_dated_releases_directory
from .ensembl_versions import is_dated_release
from .genome_fasta import _read_json, _write_json
from .species import check_species_object, find_species_by_name, Species

logger = logging.getLogger(__name__)

LISTING_TIMEOUT_SECONDS = 30


def available_dated_releases(species="human", refresh=False):
    """
    Annotation dates published for the species' current assembly, oldest first.

    Parameters
    ----------
    species : str or Species, optional
        Species name, e.g. "human" or "mus_musculus".
    refresh : bool, optional
        Fetch the dates again. Otherwise a cached listing is used without
        network access; dates are fetched and cached on first use.

    Returns
    -------
    list of str
        Dates in YYYY_MM form, each usable as ``EnsemblRelease(date)``.

    Raises
    ------
    ValueError
        If the species has no dated releases.
    OSError
        If the dates must be fetched and the server can't be reached or
        returns a page that lists no dates.
    """
    if not refresh:
        cached = cached_dated_releases(species)
        if cached is not None:
            return cached
    return _fetch_dated_releases(species, LISTING_TIMEOUT_SECONDS)


def cached_dated_releases(species="human"):
    """Dates from the last listing fetched for this species, or None; offline."""
    url, path = _listing_location(species)
    record = _read_json(path)
    if not isinstance(record, dict) or record.get("url") != url:
        return None  # From another source or unreadable; fetch again.
    dates = record.get("dates")
    if not isinstance(dates, list) or not dates or not all(map(is_dated_release, dates)):
        return None
    return dates


def _listing_location(species):
    species = check_species_object(species)
    if species.dated_releases is None:
        raise ValueError("No dated Ensembl releases for %s" % (species.latin_name,))
    accession, provider = species.dated_releases
    url = make_dated_releases_directory(accession, provider) + "/"
    path = os.path.join(cache_root(), "dated_releases", "%s_%s.json" % (accession, provider))
    return url, path


def _fetch_dated_releases(species, timeout):
    url, path = _listing_location(species)
    html = _read_listing_page(url, timeout)
    dates = parse_dated_release_listing(html)
    if not dates:
        # Every species with dated releases has at least one; this page is
        # something else (a proxy error, a new listing format), so don't
        # cache it as the truth.
        raise OSError("No annotation dates found at %s" % (url,))
    record = {
        "url": url,
        "fetched": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dates": dates,
    }
    try:
        _write_json(path, record)
    except OSError as error:  # e.g. a read-only shared cache
        logger.warning("Couldn't cache dated releases in %s: %s", path, error)
    return dates


def fetch_all_dated_releases(timeout=10):
    """
    Dates of every species with dated releases, by latin name, fetched again
    and falling back to cached dates, or None when none are cached. Human is
    fetched first: if its listing can't be read, the platform is treated as
    unreachable and every species is read from the cache, so going offline
    costs one timeout rather than one per species.
    """
    species = [
        find_species_by_name(name) for name in sorted(Species.all_registered_latin_names())
    ]
    species = [s for s in species if s.dated_releases is not None]
    species.sort(key=lambda s: s.latin_name != "homo_sapiens")  # Stable: human first.

    def fetch(s):
        try:
            return _fetch_dated_releases(s, timeout)
        except OSError as error:
            logger.debug("Dated releases of %s: %s", s.latin_name, error)
            return cached_dated_releases(s)

    try:
        dates = {species[0].latin_name: _fetch_dated_releases(species[0], timeout)}
    except OSError as error:
        logger.info("Using cached dated releases: %s", error)
        return {s.latin_name: cached_dated_releases(s) for s in species}
    with ThreadPoolExecutor(4) as pool:
        dates.update(zip((s.latin_name for s in species[1:]), pool.map(fetch, species[1:])))
    return dates


def _read_listing_page(url, timeout):
    """The listing's HTML. EBI sometimes refuses bursts of connections, so a
    refused or reset connection is retried; other failures raise OSError."""
    for delay in (0.5, 1.0, None):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read().decode("utf-8", "replace")
        except http.client.HTTPException as error:  # e.g. a truncated response
            raise OSError("Couldn't read %s: %s" % (url, error)) from error
        except OSError as error:
            reason = error.reason if isinstance(error, urllib.error.URLError) else error
            if delay is None or not isinstance(
                reason, (ConnectionRefusedError, ConnectionResetError)
            ):
                raise
            time.sleep(delay)


def parse_dated_release_listing(html):
    """Annotation dates linked from an HTML directory listing, oldest first."""
    links = re.findall(r'href="([^"/?]+)/"', html)
    return sorted({link for link in links if is_dated_release(link)})


def require_published_date(genome):
    """
    Raise ValueError before a download if Ensembl doesn't publish the dated
    release's date, instead of failing on a missing file.

    Uses cached dates, refreshing once in case the date is new. Mirrors are
    not checked, and a listing that can't be fetched is not an error here.
    """
    if genome.server != ENSEMBL_PLATFORM_FTP_SERVER:
        return
    try:
        dates = available_dated_releases(genome.species)
        if genome.release not in dates:
            dates = available_dated_releases(genome.species, refresh=True)
    except OSError:
        return  # Offline: the download reports its own error.
    if genome.release not in dates:
        accession, provider = genome.species.dated_releases
        raise ValueError(
            "Ensembl publishes no %s annotation of %s (%s, provider %s) dated %s; "
            "available dates: %s"
            % (
                genome.species.latin_name, genome.reference_name, accession, provider,
                genome.release, ", ".join(dates),
            )
        )

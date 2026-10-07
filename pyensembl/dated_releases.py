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
Which annotation dates the new Ensembl platform publishes.

Dates come from the directory listing of an assembly accession and provider,
e.g. https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/,
which is exactly where dated annotations download from. Ensembl's species.json
catalogue is 37 MB and omits some installable datasets, so it is not used.
"""

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
import logging
import os
import re
import time

from .download_cache import cache_root
from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER, make_dated_releases_directory
from .ensembl_versions import is_dated_release
from .genome_fasta import _read_json, _write_json
from .species import Species, check_species_object, find_species_by_name

logger = logging.getLogger(__name__)

LISTING_TIMEOUT_SECONDS = 30
_LISTING_RETRIES = 2
_LISTING_RETRY_MAX_DELAY = 10.0
# `pyensembl available` checks listings again once they are this old.
AVAILABLE_MAX_AGE_SECONDS = 24 * 60 * 60

_warned_unwritable_cache = False


class UnpublishedDateError(ValueError):
    """Ensembl doesn't publish an annotation with this date; the message lists
    the dates it does publish."""


def available_dated_releases(species="human", refresh=False):
    """
    Annotation dates published for the species' current assembly, oldest first.

    Parameters
    ----------
    species : str or Species, optional
        Species name, e.g. "human" or "mus_musculus".
    refresh : bool, optional
        Fetch the dates again. Otherwise cached dates are used without network
        access; dates are fetched and cached on first use.

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
    accession, provider = _dated_release_dataset(species)
    return available_annotation_dates(accession, provider, refresh=refresh)


def available_annotation_dates(assembly_accession, provider="ensembl", refresh=False):
    """
    Annotation dates published for one assembly accession and provider on the
    new Ensembl platform, oldest first; see ``EnsemblAnnotation``.

    Parameters
    ----------
    assembly_accession : str
        Versioned GCA or GCF accession, e.g. "GCA_000001405.29".
    provider : str, optional
        Provider directory, e.g. "ensembl" or "community".
    refresh : bool, optional
        Fetch the dates again rather than using cached dates.

    Raises
    ------
    OSError
        If the dates must be fetched and the server can't be reached or
        returns a page that lists no dates.
    """
    if not refresh:
        record = _cached_listing(assembly_accession, provider)
        if record is not None:
            return record["dates"]
    return _fetch_listing(assembly_accession, provider, LISTING_TIMEOUT_SECONDS)


def fetch_all_dated_releases(max_age=AVAILABLE_MAX_AGE_SECONDS, timeout=10):
    """
    Dates of every species with dated releases, by latin name, or None when
    they can't be fetched and none are cached.

    Cached dates younger than ``max_age`` seconds are used as they are; older
    ones are fetched again. Human is fetched first: if its listing can't be
    read, the platform is treated as unreachable and cached dates are used,
    so going offline costs one timeout rather than one per species.
    """
    species = [
        find_species_by_name(name) for name in sorted(Species.all_registered_latin_names())
    ]
    species = [s for s in species if s.dated_releases is not None]
    species.sort(key=lambda s: s.latin_name != "homo_sapiens")  # Stable: human first.
    records = {s.latin_name: _cached_listing(*s.dated_releases) for s in species}
    dates = {
        name: None if record is None else record["dates"] for name, record in records.items()
    }
    stale = [s for s in species if _age(records[s.latin_name]) > max_age]
    if not stale:
        return dates

    probe = species[0]
    try:
        dates[probe.latin_name] = _fetch_listing(*probe.dated_releases, timeout=timeout)
    except OSError as error:
        logger.info("Using cached dated releases: %s", error)
        return dates

    def fetch(s):
        try:
            return _fetch_listing(*s.dated_releases, timeout=timeout)
        except OSError as error:
            logger.debug("Dated releases of %s: %s", s.latin_name, error)
            return dates[s.latin_name]

    rest = [s for s in stale if s is not probe]
    with ThreadPoolExecutor(4) as pool:
        dates.update(zip((s.latin_name for s in rest), pool.map(fetch, rest)))
    return dates


@contextmanager
def explain_unpublished_date(assembly_accession, provider, date, server, description=None):
    """
    Turn a download failure into UnpublishedDateError when Ensembl doesn't
    publish the annotation date, listing the dates it does publish.

    Only a missing file (HTTP 404) on the new platform's server is explained.
    The original error propagates if the date is published (another file is
    missing) or the listing can't be read. ``description`` names the dataset,
    e.g. "homo_sapiens GRCh38".
    """
    try:
        yield
    except Exception as error:
        status = getattr(getattr(error, "response", None), "status_code", None)
        if status != 404 or server.rstrip("/") != ENSEMBL_PLATFORM_FTP_SERVER:
            raise
        try:
            dates = available_annotation_dates(assembly_accession, provider)
            if date not in dates:  # Possibly new since the dates were cached.
                dates = available_annotation_dates(assembly_accession, provider, refresh=True)
        except OSError as listing_error:
            logger.debug("Can't check annotation dates: %s", listing_error)
            dates = None
        if dates is None or date in dates:
            raise
        if description:
            dataset = "%s annotation (%s, provider %s)" % (
                description, assembly_accession, provider
            )
        else:
            dataset = "annotation of %s (provider %s)" % (assembly_accession, provider)
        raise UnpublishedDateError(
            "Ensembl publishes no %s dated %s; available dates: %s"
            % (dataset, date, ", ".join(dates))
        ) from error


def parse_dated_release_listing(html):
    """Annotation dates linked from an HTML directory listing, oldest first."""
    links = re.findall(r'href="([^"/?]+)/"', html)
    return sorted({link for link in links if is_dated_release(link)})


def _dated_release_dataset(species):
    species = check_species_object(species)
    if species.dated_releases is None:
        raise ValueError("No dated Ensembl releases for %s" % (species.latin_name,))
    return species.dated_releases


def _listing_location(assembly_accession, provider):
    url = make_dated_releases_directory(assembly_accession, provider) + "/"
    path = os.path.join(
        cache_root(), "dated_releases", "%s_%s.json" % (assembly_accession, provider)
    )
    return url, path


def _cached_listing(assembly_accession, provider):
    """The cached {"url", "fetched", "dates"} record, or None; offline."""
    url, path = _listing_location(assembly_accession, provider)
    record = _read_json(path)
    if not isinstance(record, dict) or record.get("url") != url:
        return None  # From another source or unreadable; fetch again.
    dates = record.get("dates")
    if not isinstance(dates, list) or not dates or not all(map(is_dated_release, dates)):
        return None
    return record


def _age(record):
    """Seconds since a cached listing was fetched; infinite if unknown."""
    try:
        fetched = datetime.fromisoformat(record["fetched"])
    except (TypeError, KeyError, ValueError):
        return float("inf")
    return (datetime.now(timezone.utc) - fetched).total_seconds()


def _fetch_listing(assembly_accession, provider, timeout):
    global _warned_unwritable_cache
    url, path = _listing_location(assembly_accession, provider)
    dates = parse_dated_release_listing(_get(url, timeout))
    if not dates:
        # Every dataset directory has at least one date; this page is
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
        level = logging.DEBUG if _warned_unwritable_cache else logging.WARNING
        logger.log(level, "Couldn't cache annotation dates in %s: %s", path, error)
        _warned_unwritable_cache = True
    return dates


def _get(url, timeout):
    """
    A listing page, fetched with requests (so TLS trust matches downloads)
    and retried like datacache downloads: connection errors, 408, 429 and
    5xx, honoring Retry-After. A timeout is not retried; a listing is a few
    kilobytes, so it means the server is unreachable. Errors are OSErrors.
    """
    import requests
    from datacache.retries import is_retryable_http_error, retry_delay

    backoff = 1.0
    for attempt in range(_LISTING_RETRIES + 1):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            return response.text
        except requests.Timeout:
            raise
        except requests.RequestException as error:
            delay = None
            if attempt < _LISTING_RETRIES and is_retryable_http_error(error):
                delay = retry_delay(error, backoff, _LISTING_RETRY_MAX_DELAY)
            if delay is None:
                raise
            time.sleep(delay)
            backoff *= 2

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
from datetime import timedelta
import logging
import os
import re

from .download_cache import cache_root
from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER, make_dated_releases_directory
from .ensembl_versions import is_dated_release
from .species import Species, check_species_object, find_species_by_name

logger = logging.getLogger(__name__)

LISTING_TIMEOUT_SECONDS = 30
# `pyensembl available` checks listings again once they are this old.
AVAILABLE_EXPIRE_AFTER = timedelta(days=1)

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
    return _listing_dates(assembly_accession, provider, refresh=refresh)


def fetch_all_dated_releases(expire_after=AVAILABLE_EXPIRE_AFTER, timeout=10):
    """
    Dates of every species with dated releases, by latin name, or None when
    they can't be fetched and none are cached.

    Cached dates newer than ``expire_after`` are used as they are; older ones
    are fetched again, keeping the cached dates if that fails. Human is
    fetched first: if its listing can't be read, the platform is treated as
    unreachable and every species uses its cached dates, so going offline
    costs one failed request rather than one per species.
    """
    species = [
        find_species_by_name(name) for name in sorted(Species.all_registered_latin_names())
    ]
    species = [s for s in species if s.dated_releases is not None]
    species.sort(key=lambda s: s.latin_name != "homo_sapiens")  # Stable: human first.
    probe = species[0]
    try:
        dates = {probe.latin_name: _listing_dates(
            *probe.dated_releases, timeout=timeout, expire_after=expire_after)}
    except OSError as error:
        logger.info("Using cached dated releases: %s", error)
        return {s.latin_name: _cached_dates(*s.dated_releases) for s in species}

    def fetch(s):
        try:
            return _listing_dates(
                *s.dated_releases, timeout=timeout, expire_after=expire_after,
                return_stale_on_error=True)
        except OSError as error:
            logger.debug("Dated releases of %s: %s", s.latin_name, error)
            return None

    rest = species[1:]
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
        cache_root(), "dated_releases", "%s_%s.html" % (assembly_accession, provider)
    )
    return url, path


def _dates_in(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return parse_dated_release_listing(handle.read())


def _require_dates(path):
    """datacache validator: every dataset directory lists at least one date,
    so a page without one (a proxy error, a new listing format) is rejected
    rather than cached."""
    if not _dates_in(path):
        raise ValueError("no annotation dates listed")


def _cached_dates(assembly_accession, provider):
    """Dates from the cached listing, or None; never uses the network."""
    _, path = _listing_location(assembly_accession, provider)
    try:
        return _dates_in(path) or None
    except OSError:
        return None


def _listing_dates(
    assembly_accession, provider, refresh=False, *,
    timeout=LISTING_TIMEOUT_SECONDS, expire_after=None, return_stale_on_error=False,
):
    """
    Dates from the cached listing, fetched when it's missing, when refresh is
    true, or when it's older than expire_after. Raises OSError if the listing
    can't be fetched or lists no dates.
    """
    global _warned_unwritable_cache
    import datacache

    url, path = _listing_location(assembly_accession, provider)
    try:
        path = datacache.fetch_file(
            url, destination=path, raw=True, force=refresh, timeout=timeout,
            record_provenance=True, expire_after=expire_after,
            return_stale_on_error=return_stale_on_error, validator=_require_dates,
        )
    except datacache.FileValidationError as error:
        if refresh or not os.path.exists(path):
            raise OSError("No annotation dates found at %s" % (url,)) from error
        # A cached listing without dates, e.g. from another tool: fetch it again.
        return _listing_dates(assembly_accession, provider, refresh=True, timeout=timeout)
    except PermissionError as error:
        # A read-only shared cache: read the listing without caching it.
        level = logging.DEBUG if _warned_unwritable_cache else logging.WARNING
        logger.log(level, "Couldn't cache annotation dates in %s: %s", path, error)
        _warned_unwritable_cache = True
        listing = datacache.fetch_bytes(url, timeout=timeout)
        dates = parse_dated_release_listing(listing.decode("utf-8", "replace"))
        if not dates:
            raise OSError("No annotation dates found at %s" % (url,)) from error
        return dates
    return _dates_in(path)

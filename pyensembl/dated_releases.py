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

from datetime import datetime, timezone
import json
import os
import re
import urllib.request

from .common import _atomic_output
from .download_cache import cache_root
from .ensembl_url_templates import ENSEMBL_PLATFORM_FTP_SERVER, make_dated_releases_directory
from .ensembl_versions import is_dated_release
from .species import check_species_object

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
        If the dates must be fetched and the server can't be reached.
    """
    species = check_species_object(species)
    if species.dated_releases is None:
        raise ValueError("No dated Ensembl releases for %s" % (species.latin_name,))
    accession, provider = species.dated_releases
    url = make_dated_releases_directory(accession, provider) + "/"
    path = os.path.join(cache_root(), "dated_releases", "%s_%s.json" % (accession, provider))
    if not refresh:
        cached = _read_listing(path, url)
        if cached is not None:
            return cached
    with urllib.request.urlopen(url, timeout=LISTING_TIMEOUT_SECONDS) as response:
        dates = parse_dated_release_listing(response.read().decode("utf-8", "replace"))
    record = {
        "url": url,
        "fetched": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dates": dates,
    }
    with _atomic_output(path, "w") as handle:
        json.dump(record, handle, indent=1)
    return dates


def parse_dated_release_listing(html):
    """Annotation dates linked from an HTML directory listing, oldest first."""
    links = re.findall(r'href="([^"/?]+)/"', html)
    return sorted({link for link in links if is_dated_release(link)})


def _read_listing(path, url):
    try:
        with open(path) as handle:
            record = json.load(handle)
    except (OSError, ValueError):
        return None
    if not isinstance(record, dict) or record.get("url") != url:
        return None  # From another source or unreadable; fetch again.
    dates = record.get("dates")
    if not isinstance(dates, list) or not all(is_dated_release(date) for date in dates):
        return None
    return dates


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
                genome.release, ", ".join(dates) or "none",
            )
        )

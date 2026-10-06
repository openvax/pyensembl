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

import re

MIN_ENSEMBL_RELEASE = 40
MAX_ENSEMBL_RELEASE = 116
# Ensembl Genomes (plants, fungi, metazoa, protists, bacteria) has its own
# release numbering that runs separately from the main Ensembl release.
MAX_ENSEMBL_GENOMES_RELEASE = 63


def is_dated_release(release):
    """True for an annotation date on the new Ensembl platform, e.g. "2026_04"."""
    return (
        isinstance(release, str)
        and re.fullmatch(r"\d{4}_(0[1-9]|1[0-2])", release) is not None
    )


def normalize_release(release):
    """
    A numbered Ensembl release (e.g. 116) as an int, or the annotation date
    of a dated release on the new Ensembl platform (e.g. "2026_04").
    """
    if isinstance(release, str):
        release = release.strip()
        if re.fullmatch(r"\d{4}-\d{2}(-\d{2})?", release):
            raise ValueError(
                "%r looks like an Ensembl website release label, not an "
                "annotation date; dated releases use YYYY_MM, e.g. '2026_04'"
                % (release,)
            )
        if is_dated_release(release):
            return release
        if re.fullmatch(r"\d{4}_\d{2}", release):
            raise ValueError("Invalid annotation date: %s" % (release,))
    return check_release_number(release)


def check_release_number(release):
    """
    Check to make sure a release is in the valid range of
    Ensembl releases.
    """
    # int() accepts digit separators, which would read "2026_04" as 202604.
    if isinstance(release, str) and not release.strip().isdigit():
        raise ValueError("Invalid Ensembl release: %s" % (release,))
    try:
        release = int(release)
    except (ValueError, TypeError):
        raise ValueError("Invalid Ensembl release: %s" % (release,))

    if release < MIN_ENSEMBL_RELEASE:
        raise ValueError(
            "Invalid Ensembl releases %d, must be greater than %d"
            % (release, MIN_ENSEMBL_RELEASE)
        )
    return release

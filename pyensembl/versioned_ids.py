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
Matching caller-supplied Ensembl IDs against the form an annotation stores.

A version suffix is optional. ``ENST00000269305`` means whatever version the
annotation has; ``ENST00000269305.8`` means exactly version 8, and is an error
if the annotation records another version or none at all. Ensembl GTFs store
bare IDs with separate ``*_version`` attributes, while GENCODE GTFs and newer
Ensembl FASTA headers embed the version in the ID.
"""

import re

# An Ensembl ID with a version: a plain decimal suffix without leading zeros
_VERSIONED_ID = re.compile(r"(ENS.*)\.(0|[1-9][0-9]*)")


def _split_ens_version(identifier):
    """
    Split an ENS-prefix identifier into ``(bare_id, version_int)``.

    Returns ``(identifier, None)`` for IDs that don't carry a parseable
    ENS version. Non-ENS IDs (e.g. TAIR ``AT1G01010.1``) are always
    returned as-is with version ``None`` — the ``.N`` in those is an
    isoform suffix, not a version.
    """
    match = _VERSIONED_ID.fullmatch(identifier) if isinstance(identifier, str) else None
    if match is None:
        return identifier, None
    return match.group(1), int(match.group(2))


def match_version(identifier, installed):
    """
    Choose the stored form of ``identifier`` from ``installed``.

    ``installed`` maps each stored form of one stable ID (bare or versioned)
    to the version recorded for it, or ``None`` when none is recorded.
    Returns the stored form to use, or ``None`` if the ID isn't installed.
    Raises ``ValueError`` if a supplied version differs from the recorded one
    or can't be checked, or if a bare ID matches several versions.
    """
    if identifier in installed:
        return identifier
    if not installed:
        return None
    bare, version = _split_ens_version(identifier)
    if version is None:
        if len(installed) == 1:
            return next(iter(installed))
        raise ValueError(
            "%s matches several versions: %s; pass one of them"
            % (identifier, ", ".join(sorted(installed)))
        )
    for stored, recorded in installed.items():
        if recorded == version:
            return stored
    recorded = sorted(v for v in installed.values() if v is not None)
    if not recorded:
        raise ValueError(
            "%s: no version is recorded for %s, so version %d can't be "
            "checked; pass %s without a version" % (identifier, bare, version, bare)
        )
    raise ValueError(
        "%s is not in this annotation, which has %s"
        % (identifier, ", ".join("%s.%d" % (bare, v) for v in recorded))
    )

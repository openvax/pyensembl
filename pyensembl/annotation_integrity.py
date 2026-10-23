"""Provider checksums for newly acquired dated annotation sources.

Legacy caches need no receipts. Verification is explicit and read-only; normal
reuse never fetches manifests, hashes existing files, or writes metadata.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.request import urlopen

from filelock import FileLock

from .common import _atomic_output


class AnnotationIntegrityError(ValueError):
    """A provider checksum is missing, ambiguous, or differs from a pinned one."""


def file_md5(path):
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_checksum(url):
    """Read one unambiguous checksum from the source directory's manifest."""
    directory, name = url.rsplit("/", 1)
    manifest_url = directory + "/md5sum.txt"
    with urlopen(manifest_url, timeout=60) as response:
        data = response.read(2 * 1024 * 1024 + 1)
    if len(data) > 2 * 1024 * 1024:
        raise AnnotationIntegrityError("Checksum manifest is unexpectedly large: " + manifest_url)
    matches = []
    for line in data.decode("utf-8").splitlines():
        match = re.fullmatch(r"([a-fA-F0-9]{32})\s+\*?(?:\./)?([^/]+)", line.strip())
        if match and match[2] == name:
            matches.append(match[1].lower())
    if len(matches) != 1:
        raise AnnotationIntegrityError("Expected one checksum for %s in %s" % (name, manifest_url))
    return dict(
        source=url, md5=matches[0], manifest_url=manifest_url,
        manifest_sha256=hashlib.sha256(data).hexdigest(),
        retrieved_at=datetime.now(timezone.utc).isoformat(),
    )


def receipt_path(path):
    return Path(str(path) + ".ensembl-md5.json")


def read_receipt(path, source):
    try:
        value = json.loads(receipt_path(path).read_text())
    except (OSError, ValueError):
        return None
    if (not isinstance(value, dict) or value.get("source") != source
            or not re.fullmatch(r"[a-f0-9]{32}", str(value.get("md5", "")))):
        return None
    return value


def acquire(genome, field_name, source, overwrite, show_progress):
    cache = genome.download_cache
    destination = Path(cache._cache_path(source))
    destination.parent.mkdir(parents=True, exist_ok=True)
    # The source and its receipt are published under the same lock. Neither
    # this lock nor the directory is touched by the existing-cache fast path.
    with FileLock(str(destination) + ".ensembl-download.lock", mode=-1):
        if not overwrite and cache.inspect(source).status != "missing":
            return cache.local_path_or_install_error(field_name, source)
        identity = source_checksum(source)
        pinned = genome._resolved.get("source_checksums", {}).get(source)
        if pinned and not overwrite and pinned["md5"] != identity["md5"]:
            raise AnnotationIntegrityError("Upstream checksum changed for pinned source: " + source)

        def validate(path):
            if file_md5(path) != identity["md5"]:
                raise AnnotationIntegrityError("MD5 mismatch for " + source)

        path = cache.local_path_or_install_error(
            field_name, source, download_if_missing=True, overwrite=overwrite,
            show_progress=show_progress, validator=validate,
        )
        with _atomic_output(receipt_path(path), "w") as handle:
            json.dump(identity, handle, sort_keys=True)
            handle.write("\n")
        # Keep the acquisition pin in this object as well as its on-disk receipt.
        genome._resolved = {**genome._resolved, "source_checksums": {
            **genome._resolved.get("source_checksums", {}), source: identity,
        }}
        return path


def source_checksums(genome):
    pins = dict(genome._resolved.get("source_checksums", {}))
    for source in genome._annotation_sources().values():
        identity = read_receipt(genome.download_cache.local_path(source), source)
        # A serialized pin is authoritative even if the shared file was later
        # explicitly overwritten by a different object.
        if identity is not None:
            pins.setdefault(source, identity)
    return pins


def verify_sources(genome, check_remote=False):
    """Report local and optional upstream integrity without changing any file."""
    pins = source_checksums(genome)
    reports = []
    for role, source in genome._annotation_sources().items():
        path = genome.download_cache.local_path(source)
        inspection = genome.download_cache.inspect(source)
        identity = pins.get(source)
        result = dict(role=role, source=source, path=path, status=inspection.status,
                      expected_md5=identity["md5"] if identity else None)
        if inspection.status == "available":
            if identity is None:
                result["status"] = "unrecorded"
            else:
                result["actual_md5"] = file_md5(path)
                result["status"] = "verified" if result["actual_md5"] == identity["md5"] else "changed"
        if check_remote:
            if genome._uses_annotation_checksums():
                current = source_checksum(source)
                result["upstream_md5"] = current["md5"]
                result["upstream_status"] = (
                    "unrecorded" if identity is None else
                    "unchanged" if identity["md5"] == current["md5"] else "changed"
                )
            else:
                result["upstream_status"] = "unsupported"
        reports.append(result)
    return reports

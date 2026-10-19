"""Reject package versions already present on main, a tag, or PyPI."""

import argparse
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path


VERSION_FILE = "pyensembl/version.py"
VERSION_ASSIGNMENT = re.compile(r'^__version__ = ["\'](\d+\.\d+\.\d+)["\']$', re.MULTILINE)
VERSION_NUMBER = re.compile(r"^\d+\.\d+\.\d+$")
PYPI_URL = "https://pypi.org/pypi/pyensembl/json"


def version_tuple(value):
    """Parse a stable three-part PyEnsembl version."""
    if not VERSION_NUMBER.fullmatch(value):
        raise ValueError(f"Invalid stable version: {value!r}")
    return tuple(int(part) for part in value.split("."))


def version_from_text(text):
    """Read the package version without importing its dependencies."""
    match = VERSION_ASSIGNMENT.search(text)
    if not match:
        raise ValueError(f"Missing stable __version__ in {VERSION_FILE}")
    return match.group(1)


def git(*args):
    """Return a Git command's standard output."""
    return subprocess.check_output(["git", *args], text=True).strip()


def published_version():
    """Read the most recently published version from PyPI."""
    with urllib.request.urlopen(PYPI_URL, timeout=15) as response:
        return json.load(response)["info"]["version"]


def check_version(proposed, *, base=None, tags=(), published):
    """Raise if the proposed version has already been used."""
    candidate = version_tuple(proposed)
    existing = {"PyPI": published}
    if base is not None:
        existing["current main"] = base
    existing.update({f"tag v{tag}": tag for tag in tags})
    for source, value in existing.items():
        if candidate <= version_tuple(value):
            raise ValueError(
                f"Version {proposed} must be newer than {value} ({source})"
            )


def main(argv=None):
    """Check a PR against main or a clean main checkout before deployment."""
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--base-ref", help="Freshly fetched main ref for a PR")
    group.add_argument("--deploy", action="store_true", help="Require clean, current main")
    args = parser.parse_args(argv)

    try:
        if args.deploy:
            if git("branch", "--show-current") != "main":
                raise ValueError("Deployment requires the main branch")
            if git("status", "--porcelain"):
                raise ValueError("Deployment requires a clean working tree")
            if git("rev-parse", "HEAD") != git("rev-parse", "origin/main"):
                raise ValueError("Deployment requires current origin/main")

        proposed = version_from_text(Path(VERSION_FILE).read_text())
        base = None
        if args.base_ref:
            base = version_from_text(git("show", f"{args.base_ref}:{VERSION_FILE}"))
        tags = [tag[1:] for tag in git("tag", "--list", "v*").splitlines()
                if VERSION_NUMBER.fullmatch(tag[1:])]
        check_version(proposed, base=base, tags=tags, published=published_version())
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(f"Release version check failed: {error}", file=sys.stderr)
        return 1
    print(f"Release version {proposed} is new relative to main, tags, and PyPI")
    return 0


if __name__ == "__main__":
    sys.exit(main())

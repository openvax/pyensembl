"""Explicit, offline gene-symbol aliases mapped to annotation gene IDs."""

from collections.abc import Mapping
import csv
import gzip
from hashlib import sha256
from io import TextIOWrapper
import os
import re

from .download_cache import DownloadCache, cache_root
from .versioned_ids import _split_ens_version


HGNC_COMPLETE_SET_URL = (
    "https://storage.googleapis.com/public-download-files/"
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
)


def normalize_aliases(aliases):
    """Copy an alias mapping, treating one string ID as one value."""
    if not isinstance(aliases, Mapping):
        raise TypeError("Gene aliases must be a mapping from names to gene IDs")
    result = {}
    for name, ids in aliases.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Gene alias names must be nonempty strings")
        if isinstance(ids, str):
            ids = [ids]
        try:
            ids = list(ids)
        except TypeError as error:
            raise ValueError("Gene alias values must contain gene ID strings") from error
        if any(not isinstance(gene_id, str) or not gene_id.strip() for gene_id in ids):
            raise ValueError("Gene alias values must contain nonempty gene ID strings")
        result[name] = tuple(sorted(set(ids)))
    return result


def stable_ensembl_gene_id(gene_id):
    """Strip an Ensembl version suffix, preserving other identifier formats."""
    return _split_ens_version(gene_id)[0]


class GeneNameAliases(Mapping):
    """Read-only name-to-ID mapping with optional species/source metadata.

    Ordinary mappings also work with ``Genome.genes_by_name(aliases=...)``.
    Values may name several genes: aliases are not necessarily unique.
    """

    def __init__(self, aliases, *, species=None, source=None):
        self._aliases = normalize_aliases(aliases)
        self.species = species
        self.source = source

    def __getitem__(self, name):
        return self._aliases[name]

    def __iter__(self):
        return iter(self._aliases)

    def __len__(self):
        return len(self._aliases)

    @classmethod
    def download_hgnc(cls, *, cache_directory_path=None,
                      source_url=HGNC_COMPLETE_SET_URL, overwrite=False):
        """Download and load human HGNC aliases, reusing a cached snapshot.

        Parameters
        ----------
        cache_directory_path : str or Path, optional
            Storage base directory. Defaults to aliases/homo_sapiens/hgnc
            under PyEnsembl's configured cache root. Each source URL has
            a separate subdirectory.
        source_url : str, optional
            Complete-set TSV or gzip TSV URL. Defaults to HGNC's current
            complete set; an archived snapshot or mirror can be selected.
        overwrite : bool, optional
            Download again even when this source is already cached.

        Returns
        -------
        GeneNameAliases
            Human alias mapping whose source is the cached file path.

        Notes
        -----
        Existing cached files are read without network access. Refresh is
        explicit; HGNC nomenclature is independent of the Ensembl release.
        Use from_hgnc(path) to load an already selected local snapshot.
        """
        if not isinstance(source_url, str) or "://" not in source_url:
            raise ValueError("HGNC source_url must be a URL; use from_hgnc(path) for a local file")
        if cache_directory_path is None:
            cache_directory_path = os.path.join(cache_root(), "aliases", "homo_sapiens", "hgnc")
        source_key = sha256(source_url.encode("utf-8")).hexdigest()[:16]
        cache = DownloadCache(
            reference_name="homo_sapiens", annotation_name="hgnc",
            cache_directory_path=os.path.join(os.fspath(cache_directory_path), source_key),
        )
        path = cache.download_or_copy_if_necessary(
            source_url, download_if_missing=True, overwrite=overwrite,
        )
        return cls.from_hgnc(path)

    @classmethod
    def from_hgnc(cls, path):
        """Load an HGNC complete-set TSV snapshot (plain or gzip).

        Include symbols, aliases and previous symbols of approved human
        records with an Ensembl gene ID. Do not download or modify the file.
        HGNC nomenclature may be newer than the selected Ensembl annotation.
        """
        path = os.fspath(path)
        aliases = {}
        required = {"symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id", "status"}
        with open(path, "rb") as raw:
            binary = gzip.GzipFile(fileobj=raw) if raw.peek(2)[:2] == b"\x1f\x8b" else raw
            with TextIOWrapper(binary, encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t", strict=True)
                if not required.issubset(reader.fieldnames or ()):
                    raise ValueError("HGNC TSV is missing required columns: %s" %
                                     ", ".join(sorted(required - set(reader.fieldnames or ()))))
                for row in reader:
                    if None in row or any(row[column] is None for column in required):
                        raise ValueError("Malformed HGNC TSV row %d" % reader.line_num)
                    if row["status"] != "Approved" or not row["ensembl_gene_id"].strip():
                        continue
                    ids = [value.strip() for value in row["ensembl_gene_id"].split("|")]
                    if any(not re.fullmatch(r"ENSG\d+(?:\.\d+)?", value) for value in ids):
                        raise ValueError("Invalid Ensembl gene ID in HGNC row %d" % reader.line_num)
                    names = [row["symbol"], *row["alias_symbol"].split("|"),
                             *row["prev_symbol"].split("|")]
                    for name in names:
                        name = name.strip()
                        if name:
                            aliases.setdefault(name, set()).update(ids)
        return cls(aliases, species="homo_sapiens", source=os.path.abspath(path))

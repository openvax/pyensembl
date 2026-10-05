"""Explicit, offline gene-symbol aliases mapped to annotation gene IDs."""

from collections.abc import Mapping
import csv
import gzip
import os
import re


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
    match = re.fullmatch(r"(ENS[A-Z]*G\d+)\.\d+", gene_id)
    return match.group(1) if match else gene_id


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
    def from_hgnc(cls, path):
        """Load an HGNC complete-set TSV snapshot (plain or gzip).

        Include symbols, aliases and previous symbols of approved human
        records with an Ensembl gene ID. Do not download or modify the file.
        HGNC nomenclature may be newer than the selected Ensembl annotation.
        """
        path = os.fspath(path)
        opener = gzip.open if path.endswith(".gz") else open
        aliases = {}
        required = {"symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id", "status"}
        with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
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

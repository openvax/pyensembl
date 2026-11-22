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

"""Pinned rapid-release and non-numbered Ensembl assemblies (e.g. T2T-CHM13v2.0).

Some Ensembl gene annotations, such as human ``T2T-CHM13v2.0``
(``GCA_009914755.4`` / ``hs1``, Ensembl rapid-release ``2022_07``), are not
published under numbered ``release-*`` directories or the newer Ensembl
Platform dated-release tree. This module catalogues those assemblies with
SHA-256 and byte-size pins, assembly/contig aliases, explicit rejection of
incompatible assembly versions (such as ``T2T-CHM13v1.0``), and
mitochondrion disambiguation (CHM13's native ``CP068254.1`` vs the rCRS
``NC_012920.1`` carried by ``chm13v2.0_maskedY_rCRS``).
"""

from copy import deepcopy
from dataclasses import dataclass, field
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import re
from typing import Mapping, Optional
from weakref import WeakValueDictionary

from filelock import FileLock

from .annotation_integrity import AnnotationIntegrityError
from .common import _atomic_output, file_digest
from .download_cache import _default_cache_directory, inspect_data_file
from .genome import Genome
from .species import check_species_object

# The revised Cambridge Reference Sequence (GRCh38's MT).
RCRS_ACCESSION = "NC_012920.1"

_ASSET_ROLES = ("gtf", "transcripts", "proteins")


@dataclass(frozen=True)
class ContigAlias:
    """One contig in a pinned assembly across annotation, UCSC, GenBank, and RefSeq names."""

    name: str
    ucsc: Optional[str]
    genbank: Optional[str]
    refseq: Optional[str]
    length: int


@dataclass(frozen=True)
class ReferenceAlias:
    """One reference name or alias associated with a pinned assembly (or an unsupported version)."""

    name: str
    assembly: Optional[str]
    reference_name: Optional[str]
    assembly_accession: Optional[str]
    mitochondrion: Optional[str]
    without: tuple[str, ...] = ()
    note: str = ""


@dataclass(frozen=True)
class PinnedAssembly:
    """A pinned non-numbered / rapid-release Ensembl assembly with checksum-verified assets."""

    name: str
    species: str
    reference_name: str
    annotation_name: str
    annotation_version: str
    title: str
    assets: Mapping = field(hash=False)
    provider: str
    terms: str
    terms_url: str
    assembly_accession: Optional[str] = None
    names: Mapping = field(default_factory=dict, hash=False)
    other_versions: Mapping = field(default_factory=dict, hash=False)
    contigs: tuple[ContigAlias, ...] = field(default=(), hash=False)
    contigs_source: Mapping = field(default_factory=dict, hash=False)
    mitochondria: Mapping = field(default_factory=dict, hash=False)


def _normalize_alias(name):
    """Lowercase, convert dashes/spaces to underscores, and join version suffixes.

    ``T2T-CHM13 v1.0`` and ``t2t-chm13.v1.0`` normalize to ``t2t_chm13v1.0``,
    preventing accidental matches against ``t2t-chm13``.
    """
    value = name.strip().lower().replace("-", "_").replace(" ", "_")
    return re.sub(r"[_.]+(?=v\d+\.\d)", "", value)


@lru_cache(maxsize=1)
def _raw_catalog():
    return json.loads(files("pyensembl").joinpath("pinned_assemblies.json").read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _all_pinned_assemblies():
    result = []
    for name, entry in sorted(_raw_catalog().items()):
        contigs = tuple(
            ContigAlias(
                name=c["name"],
                ucsc=c.get("ucsc"),
                genbank=c.get("genbank"),
                refseq=c.get("refseq"),
                length=int(c["length"]),
            )
            for c in entry.get("contigs", ())
        )
        kwargs = {k: v for k, v in entry.items() if k != "contigs"}
        result.append(PinnedAssembly(name=name, contigs=contigs, **kwargs))
    return tuple(result)


def pinned_assemblies(*, species=None):
    """Return all catalogued pinned assemblies, optionally filtered by ``species``."""
    items = _all_pinned_assemblies()
    if species is None:
        return items
    target = check_species_object(species).latin_name
    return tuple(item for item in items if item.species == target)


@lru_cache(maxsize=1)
def _all_reference_aliases():
    rows = []
    for assembly in _all_pinned_assemblies():
        default_mito = next((c.genbank for c in assembly.contigs if c.name == "MT"), None)
        own = {
            assembly.name: "pyensembl pinned assembly name",
            assembly.reference_name: "NCBI assembly name",
        }
        if assembly.assembly_accession:
            own[assembly.assembly_accession] = "Assembly accession"
        for label, note in own.items():
            if label not in assembly.names:
                rows.append(
                    ReferenceAlias(
                        name=label,
                        assembly=assembly.name,
                        reference_name=assembly.reference_name,
                        assembly_accession=assembly.assembly_accession,
                        mitochondrion=default_mito,
                        without=(),
                        note=note,
                    )
                )
        for label, entry in assembly.names.items():
            rows.append(
                ReferenceAlias(
                    name=label,
                    assembly=assembly.name,
                    reference_name=assembly.reference_name,
                    assembly_accession=assembly.assembly_accession,
                    mitochondrion=entry.get("mitochondrion"),
                    without=tuple(entry.get("without", ())),
                    note=entry.get("note", ""),
                )
            )
        for label, reason in assembly.other_versions.items():
            rows.append(
                ReferenceAlias(
                    name=label,
                    assembly=None,
                    reference_name=None,
                    assembly_accession=None,
                    mitochondrion=None,
                    without=(),
                    note=reason,
                )
            )
    return tuple(rows)


def pinned_reference_aliases(assembly=None):
    """Return known reference aliases (and rejected draft/older versions) for pinned assemblies.

    Parameters
    ----------
    assembly : str, optional
        If provided, return aliases belonging to the specified pinned assembly
        (plus its rejected ``other_versions``).
    """
    rows = _all_reference_aliases()
    if assembly is None:
        return rows
    resolved = resolve_pinned_assembly(assembly)
    other = set(resolved.other_versions)
    return tuple(
        row for row in rows
        if row.assembly == resolved.name or (row.assembly is None and row.name in other)
    )


@lru_cache(maxsize=1)
def _alias_index():
    return {_normalize_alias(row.name): row for row in _all_reference_aliases()}


@lru_cache(maxsize=1)
def _alias_patterns():
    rows = []
    for row in _all_reference_aliases():
        key = _normalize_alias(row.name)
        pattern = re.compile(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(key))
        rows.append((pattern, row))
    return tuple(sorted(rows, key=lambda item: -len(item[1].name)))


def lookup_pinned_alias(name):
    """Return the :class:`ReferenceAlias` matching ``name`` (case/punctuation-insensitive), or None."""
    if not isinstance(name, str) or not name.strip():
        return None
    return _alias_index().get(_normalize_alias(name))


def find_pinned_reference_name(text):
    """Return the longest whole-word pinned reference name occurring in ``text``, or None.

    Matches ignore case, dashes, and spaces, and only match whole alphanumeric
    tokens (so ``hs1`` matches ``hs1.fa`` and ``/refs/hs1/``, never ``hs10``).
    The longest match wins, so ``chm13v2.0_maskedY_rCRS.fa`` returns
    ``"chm13v2.0_maskedY_rCRS"`` rather than ``"chm13v2.0"``.
    """
    if not isinstance(text, str) or not text.strip():
        return None
    normalized = _normalize_alias(text)
    for pattern, row in _alias_patterns():
        if pattern.search(normalized):
            return row.name
    return None


def resolve_pinned_assembly(name):
    """Resolve any pinned assembly name or alias to its :class:`PinnedAssembly`.

    Raises
    ------
    TypeError
        If ``name`` is not a string.
    ValueError
        If ``name`` identifies an unsupported version of a pinned assembly
        (such as ``T2T-CHM13v1.0``, with a detailed explanation of why its
        coordinates differ) or is not a known pinned assembly.
    """
    if not isinstance(name, str):
        raise TypeError("Pinned assembly name must be a string, not %s" % type(name).__name__)
    row = lookup_pinned_alias(name)
    if row is None:
        known = ", ".join(item.reference_name for item in _all_pinned_assemblies()) or "none"
        raise ValueError("Unknown pinned assembly %r (known: %s)" % (name, known))
    if row.assembly is None:
        raise ValueError("%s is not a supported reference: %s" % (name, row.note))
    for item in _all_pinned_assemblies():
        if item.name == row.assembly:
            return item
    raise ValueError("Unknown pinned assembly %r" % (name,))


def is_rcrs_reference(name):
    """Return True if files named ``name`` carry the rCRS (``NC_012920.1``) as their mitochondrion."""
    row = lookup_pinned_alias(name)
    return row is not None and row.mitochondrion == RCRS_ACCESSION


def contig_aliases(name, *, chrm_is_rcrs=None):
    """Return the :class:`ContigAlias` table for a pinned assembly or alias.

    Parameters
    ----------
    name : str
        Any supported pinned reference name or alias, e.g. ``"T2T-CHM13v2.0"``,
        ``"hs1"``, ``"chm13v2.0_noY"``, or ``"chm13v2.0_maskedY_rCRS"``.
    chrm_is_rcrs : bool, optional
        Override whether ``chrM`` refers to the rCRS (``NC_012920.1``) instead
        of the assembly's native mitochondrion (``CP068254.1``). When omitted,
        inferred from ``name`` (True for ``_rCRS`` analysis sets).
    """
    assembly = resolve_pinned_assembly(name)
    row = lookup_pinned_alias(name)
    without = set(row.without) if row is not None else set()
    use_rcrs = (row is not None and row.mitochondrion == RCRS_ACCESSION) if chrm_is_rcrs is None else bool(chrm_is_rcrs)
    if use_rcrs:
        without.add("MT")
    elif chrm_is_rcrs is False:
        without.discard("MT")
    contigs = [c for c in assembly.contigs if c.name not in without]
    if use_rcrs and RCRS_ACCESSION in assembly.mitochondria:
        rcrs = assembly.mitochondria[RCRS_ACCESSION]["contig"]
        contigs.append(
            ContigAlias(
                name=rcrs["name"],
                ucsc=rcrs.get("ucsc"),
                genbank=rcrs.get("genbank"),
                refseq=rcrs.get("refseq"),
                length=int(rcrs["length"]),
            )
        )
    return tuple(contigs)


@lru_cache(maxsize=None)
def _contig_maps(reference_name, chrm_is_rcrs):
    rows = contig_aliases(reference_name, chrm_is_rcrs=bool(chrm_is_rcrs))
    names = {}
    lengths = {}
    for row in rows:
        for alias in (row.name, row.ucsc, row.genbank, row.refseq):
            if alias:
                names[alias] = row.name
        if row.ucsc == "chrM":
            names["M"] = row.name
        lengths[row.name] = row.length
    return names, lengths


def normalize_pinned_contig(contig, reference_name, *, chrm_is_rcrs=False):
    """Map any UCSC/GenBank/RefSeq/annotation contig spelling to the pinned annotation's name.

    Returns ``None`` if ``reference_name`` is not a supported pinned assembly or
    if ``contig`` is not a known contig in that assembly.
    """
    row = lookup_pinned_alias(reference_name)
    if row is None or row.assembly is None:
        return None
    names = _contig_maps(row.reference_name, bool(chrm_is_rcrs))[0]
    contig_str = str(contig)
    if contig_str in names:
        return names[contig_str]
    match = re.fullmatch(r"chr([0-9]+|[xym]|mt)", contig_str, flags=re.IGNORECASE)
    if match:
        return names.get("chr" + match[1].upper())
    return names.get(contig_str.upper())


def pinned_contig_spellings(reference_name, *, chrm_is_rcrs=False, include_mitochondrion=True):
    """Return ``{annotation_contig: set_of_spellings}`` for a pinned assembly."""
    assembly = resolve_pinned_assembly(reference_name)
    names, _ = _contig_maps(assembly.reference_name, bool(chrm_is_rcrs))
    excluded = {names.get("chrM")} if not include_mitochondrion else set()
    spellings = {}
    for alias, canonical in names.items():
        if canonical not in excluded:
            spellings.setdefault(canonical, set()).add(alias)
    return spellings


def pinned_contig_lengths(reference_name, *, chrm_is_rcrs=False):
    """Return ``{annotation_contig: length}`` for a pinned assembly."""
    assembly = resolve_pinned_assembly(reference_name)
    return dict(_contig_maps(assembly.reference_name, bool(chrm_is_rcrs))[1])


def mitochondrion_note(contig, reference_name):
    """Explain why a pinned assembly's mitochondrion has no annotation, or return None."""
    row = lookup_pinned_alias(reference_name)
    if row is None or row.assembly is None:
        return None
    contig_str = str(contig)
    if contig_str not in ("MT", "chrM", "M", RCRS_ACCESSION):
        return None
    canonical = row.reference_name
    if contig_str == RCRS_ACCESSION:
        return (
            "this %s annotation has no genes on the rCRS; PyEnsembl's pinned %s "
            "annotates %s with Ensembl GRCh38's MT genes"
            % (canonical, row.assembly, RCRS_ACCESSION)
        )
    return (
        "%s's own mitochondrion (CP068254.1) has no annotation: Ensembl's 2022_07 "
        "annotation of %s has no mitochondrial genes, and its position 1 corresponds "
        "to rCRS position 577, so it cannot be annotated as the rCRS. Only the "
        "chm13v2.0_maskedY_rCRS analysis set carries the rCRS as chrM, annotated on %s"
        % (canonical, canonical, RCRS_ACCESSION)
    )


def _packaged_bytes(resource, expected_sha256, expected_size=None):
    raw = files("pyensembl").joinpath(resource).read_bytes()
    if expected_size is not None and len(raw) != int(expected_size):
        raise AnnotationIntegrityError(
            "Packaged resource %s has size %d, expected %d" % (resource, len(raw), int(expected_size))
        )
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected_sha256:
        raise AnnotationIntegrityError(
            "Packaged resource %s has SHA-256 %s, expected %s" % (resource, actual, expected_sha256)
        )
    return raw


def _validate_pinned_file(path, expected_sha256, expected_size, label):
    actual_size = os.path.getsize(path)
    if expected_size is not None and actual_size != int(expected_size):
        raise AnnotationIntegrityError(
            "Size mismatch for %s: got %d bytes, expected %d bytes"
            % (label, actual_size, int(expected_size))
        )
    actual_sha256 = file_digest(path, "sha256")
    if actual_sha256 != expected_sha256:
        raise AnnotationIntegrityError(
            "SHA-256 mismatch for %s: got %s, expected %s"
            % (label, actual_sha256, expected_sha256)
        )


class PinnedGenome(Genome):
    """A :class:`~pyensembl.Genome` backed by a checksum-pinned rapid-release assembly.

    Supports offline alias resolution (e.g. ``T2T-CHM13v2.0``, ``hs1``,
    ``GCA_009914755.4``, ``t2t-chm13``, ``chm13v2.0_maskedY_rCRS``), SHA-256
    and byte-size verification on download and via :meth:`verify_sources`,
    and bundled rCRS (``NC_012920.1``) mitochondrial annotation alongside
    Ensembl's rapid-release nuclear annotation.
    """

    _genome_cache = WeakValueDictionary()

    @classmethod
    def cached(
        cls,
        name="T2T-CHM13v2.0",
        *,
        genome_fasta_path_or_url=None,
        cache_directory_path=None,
        reference_metadata=None,
    ):
        """Construct a :class:`PinnedGenome` or return an existing cached instance."""
        genome = cls(
            name=name,
            genome_fasta_path_or_url=genome_fasta_path_or_url,
            cache_directory_path=cache_directory_path,
            reference_metadata=reference_metadata,
        )
        cache_key = (
            cls,
            genome.pinned_name,
            genome.download_cache.cache_directory_path,
            genome._genome_fasta_path_or_url,
            genome._reference_metadata_key,
        )
        if cache_key in cls._genome_cache:
            return cls._genome_cache[cache_key]
        cls._genome_cache[cache_key] = genome
        return genome

    def __init__(
        self,
        name="T2T-CHM13v2.0",
        *,
        genome_fasta_path_or_url=None,
        cache_directory_path=None,
        reference_metadata=None,
    ):
        assembly = resolve_pinned_assembly(name)
        self.pinned_assembly = assembly
        self.pinned_name = assembly.name
        self.assembly = assembly.reference_name
        self.release = assembly.annotation_version
        self._explicit_reference_metadata = deepcopy(reference_metadata) if reference_metadata is not None else None

        cache_dir = (
            os.path.abspath(os.fspath(cache_directory_path))
            if cache_directory_path is not None
            else _default_cache_directory(
                assembly.reference_name,
                assembly.annotation_name,
                assembly.annotation_version,
            )
        )
        self._raw_asset_urls = {role: assembly.assets[role]["url"] for role in _ASSET_ROLES}

        gtf_filename = self._raw_asset_urls["gtf"].rsplit("/", 1)[-1]
        if assembly.mitochondria:
            stem = gtf_filename[:-7] if gtf_filename.endswith(".gtf.gz") else gtf_filename
            gtf_source = os.path.join(cache_dir, "%s+%s.gtf.gz" % (stem, "+".join(assembly.mitochondria)))
            extra_transcripts = [
                os.path.join(cache_dir, m["assets"]["transcripts"]["resource"].rsplit("/", 1)[-1])
                for m in assembly.mitochondria.values()
            ]
            extra_proteins = [
                os.path.join(cache_dir, m["assets"]["proteins"]["resource"].rsplit("/", 1)[-1])
                for m in assembly.mitochondria.values()
            ]
        else:
            gtf_source = self._raw_asset_urls["gtf"]
            extra_transcripts = []
            extra_proteins = []

        merged_metadata = deepcopy(reference_metadata) if reference_metadata is not None else {}
        if assembly.mitochondria:
            compartments = dict(merged_metadata.get("contig_compartments") or {})
            for m in assembly.mitochondria.values():
                compartments.setdefault(m["contig"]["name"], "mitochondrial")
            merged_metadata["contig_compartments"] = compartments

        super().__init__(
            reference_name=assembly.reference_name,
            annotation_name=assembly.annotation_name,
            annotation_version=assembly.annotation_version,
            gtf_path_or_url=gtf_source,
            transcript_fasta_paths_or_urls=[self._raw_asset_urls["transcripts"], *extra_transcripts],
            protein_fasta_paths_or_urls=[self._raw_asset_urls["proteins"], *extra_proteins],
            cache_directory_path=cache_directory_path,
            genome_fasta_path_or_url=genome_fasta_path_or_url,
            reference_metadata=merged_metadata or None,
            species=assembly.species,
            assembly_accession=assembly.assembly_accession,
        )

    def install_string(self):
        return "pyensembl install --reference-name %s" % self.reference_name

    def _genome_fasta_setup_hint(self):
        return "Pass genome_fasta_path_or_url to PinnedGenome or --genome-fasta-path to pyensembl install."

    def __str__(self):
        return "PinnedGenome(name=%r, reference_name=%r, annotation_version=%r)" % (
            self.pinned_name,
            self.reference_name,
            self.annotation_version,
        )

    def __repr__(self):
        return str(self)

    def to_dict(self):
        return dict(
            name=self.pinned_name,
            genome_fasta_path_or_url=self._genome_fasta_path_or_url,
            cache_directory_path=self.cache_directory_path,
            reference_metadata=self._explicit_reference_metadata,
        )

    def _acquire_pinned_asset(self, role, *, overwrite=False, show_progress=False, source_path=None):
        spec = self.pinned_assembly.assets[role]
        url = spec["url"]
        expected_sha256 = spec["sha256"]
        expected_size = spec["size"]
        field_name = {"gtf": "gtf", "transcripts": "transcript-fasta", "proteins": "protein-fasta"}[role]

        destination = Path(self.download_cache._cache_path(url))
        destination.parent.mkdir(parents=True, exist_ok=True)

        def validate(path):
            _validate_pinned_file(path, expected_sha256, expected_size, url)

        with FileLock(str(destination) + ".pinned-download.lock", mode=-1):
            if source_path is not None:
                src = Path(source_path)
                if not src.is_file():
                    raise FileNotFoundError("Pinned %s source path does not exist: %s" % (role, src))
                validate(str(src))
                if overwrite or not destination.is_file() or file_digest(str(destination), "sha256") != expected_sha256:
                    with open(src, "rb") as reader, _atomic_output(destination, "wb") as writer:
                        while chunk := reader.read(1 << 20):
                            writer.write(chunk)
                return str(destination)

            if not overwrite and self.download_cache.inspect(url).status == "available":
                validate(str(destination))
                return str(destination)

            return self.download_cache.local_path_or_install_error(
                field_name=field_name,
                path_or_url=url,
                download_if_missing=True,
                overwrite=overwrite,
                show_progress=show_progress,
                validator=validate,
            )

    def _materialize_mitochondria(self, raw_gtf_path, *, overwrite=False):
        if not self.pinned_assembly.mitochondria:
            return
        cache_dir = Path(self.download_cache.cache_directory_path)
        cache_dir.mkdir(parents=True, exist_ok=True)
        gtf_layers = []
        for mito in self.pinned_assembly.mitochondria.values():
            for role, asset in mito["assets"].items():
                raw = _packaged_bytes(asset["resource"], asset["sha256"], asset.get("size"))
                if role == "gtf":
                    gtf_layers.append(raw)
                    continue
                target = cache_dir / asset["resource"].rsplit("/", 1)[-1]
                if overwrite or not target.is_file() or target.stat().st_size != len(raw):
                    with _atomic_output(target, "wb") as handle:
                        handle.write(raw)

        merged = Path(self._gtf_path_or_url)
        expected_merged_size = os.path.getsize(raw_gtf_path) + sum(len(layer) for layer in gtf_layers)
        if overwrite or not merged.is_file() or merged.stat().st_size != expected_merged_size:
            with open(raw_gtf_path, "rb") as source, _atomic_output(merged, "wb") as handle:
                while chunk := source.read(1 << 20):
                    handle.write(chunk)
                for layer in gtf_layers:
                    handle.write(layer)

    def download(self, overwrite=False, show_progress=False, *, source_paths=None):
        """Download (or import from ``source_paths``) and SHA-256-verify all pinned assets."""
        if source_paths is not None:
            unknown = set(source_paths) - set(_ASSET_ROLES)
            if unknown:
                raise ValueError(
                    "Unknown pinned source_paths keys %s; expected a subset of %s"
                    % (sorted(unknown), list(_ASSET_ROLES))
                )
        raw_paths = {}
        for role in _ASSET_ROLES:
            raw_paths[role] = self._acquire_pinned_asset(
                role,
                overwrite=overwrite,
                show_progress=show_progress,
                source_path=source_paths.get(role) if source_paths else None,
            )
        self._materialize_mitochondria(raw_paths["gtf"], overwrite=overwrite)
        self._set_local_paths(download_if_missing=False, overwrite=False, show_progress=show_progress)
        if self.requires_genome_fasta:
            self.download_genome_fasta(overwrite=overwrite, show_progress=show_progress)

    def verify_sources(self, check_remote=False):
        """Verify SHA-256 and byte size of every pinned annotation file offline."""
        reports = []
        for role in _ASSET_ROLES:
            spec = self.pinned_assembly.assets[role]
            url = spec["url"]
            path = self.download_cache.local_path(url)
            inspection = inspect_data_file(path)
            result = dict(
                role=role,
                source=url,
                path=path,
                status=inspection.status,
                expected_sha256=spec["sha256"],
                expected_size=spec["size"],
            )
            if inspection.status == "available":
                actual_size = os.path.getsize(path)
                actual_sha256 = file_digest(path, "sha256")
                result["actual_size"] = actual_size
                result["actual_sha256"] = actual_sha256
                result["status"] = (
                    "verified"
                    if actual_size == int(spec["size"]) and actual_sha256 == spec["sha256"]
                    else "changed"
                )
            if check_remote:
                result["upstream_status"] = "pinned"
            reports.append(result)

        cache_dir = Path(self.download_cache.cache_directory_path)
        for accession, mito in self.pinned_assembly.mitochondria.items():
            for role in ("transcripts", "proteins"):
                asset = mito["assets"][role]
                path = str(cache_dir / asset["resource"].rsplit("/", 1)[-1])
                inspection = inspect_data_file(path)
                result = dict(
                    role="%s_%s" % (accession, role),
                    source="package:%s" % asset["resource"],
                    path=path,
                    status=inspection.status,
                    expected_sha256=asset["sha256"],
                    expected_size=asset["size"],
                )
                if inspection.status == "available":
                    actual_size = os.path.getsize(path)
                    actual_sha256 = file_digest(path, "sha256")
                    result["actual_size"] = actual_size
                    result["actual_sha256"] = actual_sha256
                    result["status"] = (
                        "verified"
                        if actual_size == int(asset["size"]) and actual_sha256 == asset["sha256"]
                        else "changed"
                    )
                if check_remote:
                    result["upstream_status"] = "pinned"
                reports.append(result)

        if self.pinned_assembly.mitochondria:
            merged_path = self._gtf_path_or_url
            raw_gtf_path = self.download_cache.local_path(self.pinned_assembly.assets["gtf"]["url"])
            inspection = inspect_data_file(merged_path)
            gtf_layers = [
                _packaged_bytes(m["assets"]["gtf"]["resource"], m["assets"]["gtf"]["sha256"], m["assets"]["gtf"].get("size"))
                for m in self.pinned_assembly.mitochondria.values()
            ]
            expected_merged_size = int(self.pinned_assembly.assets["gtf"]["size"]) + sum(
                len(layer) for layer in gtf_layers
            )
            result = dict(
                role="merged_gtf",
                source=merged_path,
                path=merged_path,
                status=inspection.status,
                expected_size=expected_merged_size,
            )
            if inspection.status == "available":
                actual_size = os.path.getsize(merged_path)
                result["actual_size"] = actual_size
                actual_sha256 = file_digest(merged_path, "sha256")
                result["actual_sha256"] = actual_sha256
                if os.path.isfile(raw_gtf_path) and os.path.getsize(raw_gtf_path) == int(self.pinned_assembly.assets["gtf"]["size"]):
                    expected_hasher = hashlib.sha256()
                    with open(raw_gtf_path, "rb") as reader:
                        while chunk := reader.read(1 << 20):
                            expected_hasher.update(chunk)
                    for layer in gtf_layers:
                        expected_hasher.update(layer)
                    expected_merged_sha256 = expected_hasher.hexdigest()
                    result["expected_sha256"] = expected_merged_sha256
                    result["status"] = (
                        "verified"
                        if actual_size == expected_merged_size and actual_sha256 == expected_merged_sha256
                        else "changed"
                    )
                else:
                    result["status"] = "changed"
            if check_remote:
                result["upstream_status"] = "pinned"
            reports.append(result)
        return reports


def cached_pinned_genomes():
    """Discover pinned assemblies already present in the default cache root (offline)."""
    found = []
    for assembly in _all_pinned_assemblies():
        genome = PinnedGenome(assembly.name)
        directory = genome.download_cache.cache_directory_path
        if os.path.isdir(directory):
            found.append(genome)
    return found

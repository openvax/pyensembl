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
Contains the EnsemblRelease class, which extends the Genome class
to be specific to a particular release of Ensembl.
"""
from weakref import WeakValueDictionary
import json
import os
import shlex
import warnings

from .dated_releases import explain_unpublished_date
from .assembly import (
    dated_assembly, normalize_species, numbered_assembly, validate_assembly,
)
from .download_cache import _default_cache_directory
from .genome import Genome
from .genome_fasta_cache import canonical_source_options, release_genome_fasta
from .genome_fasta import GenomeFasta
from .ensembl_versions import MAX_ENSEMBL_RELEASE, is_dated_release, normalize_release
from .species import check_species_object, human

from .ensembl_url_templates import (
    ENSEMBL_FTP_SERVER,
    ENSEMBL_GENOMES_FTP_SERVER,
    ENSEMBL_PLATFORM_FTP_SERVER,
    dated_genome_fasta_mask,
    make_dated_release_urls,
    make_gtf_url,
    make_fasta_url,
    make_genome_fasta_url,
)


def _default_server_for_species(species):
    """Return the Ensembl FTP server appropriate for ``species``."""
    return ENSEMBL_GENOMES_FTP_SERVER if species.ensembl_genomes else ENSEMBL_FTP_SERVER


def _genome_fasta_option(genome_fasta, download_genome_fasta=None, genome_fasta_path=None):
    """Normalize reference DNA to None, True (Ensembl DNA), or an absolute path.

    Also accepts the 2.11.0 keywords download_genome_fasta/genome_fasta_path.
    """
    if download_genome_fasta is not None or genome_fasta_path is not None:
        warnings.warn(
            "download_genome_fasta= and genome_fasta_path= are deprecated; use "
            "genome_fasta=True for Ensembl DNA or genome_fasta=<path> for a local FASTA",
            DeprecationWarning,
            stacklevel=3,
        )
        if genome_fasta is not None:
            raise ValueError("Pass reference DNA only as genome_fasta")
        # 2.11.0 let a local path win and accepted any truthy download flag.
        genome_fasta = genome_fasta_path or bool(download_genome_fasta)
    if genome_fasta is None or genome_fasta is False:
        return None
    if genome_fasta is True:
        return True
    try:
        path = os.fspath(genome_fasta)
    except TypeError:
        path = None
    if not isinstance(path, str):
        raise TypeError("genome_fasta must be True, a local FASTA path, or None")
    if not path:
        raise ValueError("genome_fasta path must not be empty")
    if "://" in path:
        raise ValueError("genome_fasta must be True or a local path; use Genome for custom URLs")
    return os.path.abspath(path)


class EnsemblRelease(Genome):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database.

    A release is either numbered or an annotation date, e.g.
    ``EnsemblRelease("2026_04", species="human")``. Assembly identity is inferred
    from verified metadata; ``assembly=`` can select or check a particular
    build, full name, or versioned accession. Existing reference labels and
    cache layouts remain independent of the new assembly metadata.
    """

    @classmethod
    def normalize_init_values(cls, release, species, server):
        """
        Normalizes the arguments which uniquely specify an EnsemblRelease
        genome.
        """
        release = normalize_release(release)
        species = check_species_object(species)
        if is_dated_release(release):
            # Neither legacy server hosts dated releases.
            if server in (None, ENSEMBL_FTP_SERVER, ENSEMBL_GENOMES_FTP_SERVER):
                server = ENSEMBL_PLATFORM_FTP_SERVER
        elif server is None or server == ENSEMBL_FTP_SERVER:
            # Promote to the Ensembl Genomes server when the species lives
            # there; otherwise keep the main Ensembl server.
            server = _default_server_for_species(species)
        return (release, species, server)

    # Using a WeakValueDictionary instead of an ordinary dict to prevent a
    # memory leak in cases where we test many different releases in sequence.
    # When all the references to a particular EnsemblRelease die then that
    # genome should also be removed from this cache.
    _genome_cache = WeakValueDictionary()

    @classmethod
    def cached(
        cls, release=MAX_ENSEMBL_RELEASE, species=None, server=ENSEMBL_FTP_SERVER,
        *, genome_fasta=None, genome_fasta_type="toplevel", genome_fasta_mask="none",
        download_genome_fasta=None, genome_fasta_path=None, **selection,
    ):
        """
        Construct EnsemblRelease if it's never been made before, otherwise
        return an old instance.
        """
        genome = cls(
            release, species, server, genome_fasta=genome_fasta,
            genome_fasta_type=genome_fasta_type,
            genome_fasta_mask=genome_fasta_mask,
            download_genome_fasta=download_genome_fasta,
            genome_fasta_path=genome_fasta_path, **selection,
        )
        # Construction is offline. Key the fully resolved sources and cache,
        # including DNA attachment, rather than only release/species.
        cache_key = (
            cls, genome._fields(), genome.species, genome.server,
            genome.download_cache.cache_directory_path, genome._genome_fasta_option,
            genome.genome_fasta_type, genome.genome_fasta_mask,
            json.dumps(genome._resolved, sort_keys=True),
        )
        if cache_key in cls._genome_cache:
            return cls._genome_cache[cache_key]
        cls._genome_cache[cache_key] = genome
        return genome

    def __init__(
        self, release=MAX_ENSEMBL_RELEASE, species=None, server=ENSEMBL_FTP_SERVER,
        *, genome_fasta=None, genome_fasta_type="toplevel", genome_fasta_mask="none",
        download_genome_fasta=None, genome_fasta_path=None,
        assembly=None, assembly_accession=None, provider=None, include_alt=None,
        reference_name=None, cache_directory_path=None,
        _resolved=None, _preserve_legacy=False,
    ):
        """
        release : int or str
            Numbered Ensembl release, e.g. 116, or an annotation date on the
            new Ensembl platform in YYYY_MM form, e.g. "2026_04".

        genome_fasta : True or path, optional
            Reference DNA for ``sequence()``: True for Ensembl's DNA for this
            release (see genome_fasta_type and genome_fasta_mask), or a local
            plain or gzip FASTA. Nothing is downloaded until
            ``download_genome_fasta()``, ``download()``, or ``pyensembl install``.
        """
        release = normalize_release(release)
        self.release = release
        self.species = normalize_species(species, assembly, assembly_accession)
        if _resolved is not None and species is None:
            self.species = None
        dated = is_dated_release(release)
        if dated:
            if server in (None, ENSEMBL_FTP_SERVER, ENSEMBL_GENOMES_FTP_SERVER):
                server = ENSEMBL_PLATFORM_FTP_SERVER
        elif server in (None, ENSEMBL_FTP_SERVER):
            if self.species is None:
                raise ValueError("A numbered release requires a known species")
            server = _default_server_for_species(self.species)
        if not dated and self.species is None:
            raise ValueError("A numbered release requires a known species")
        # Numbered releases historically preserved the exact mirror URL,
        # including trailing slashes. Changing it changes source/cache identity.
        self.server = server
        if include_alt is not None and not isinstance(include_alt, bool):
            raise TypeError("include_alt must be bool")
        if _resolved is not None:
            # A serialized selection is a snapshot, not another catalogue query.
            self._resolved = dict(_resolved)
            if self._resolved.get("release") != release:
                raise ValueError("Serialized assembly selection has a different release")
        else:
            if dated:
                infer_historical = not _preserve_legacy
                if self.species is not None and not (assembly or assembly_accession):
                    # Existing installations and mirrors keep the old selection.
                    existing = _default_cache_directory(
                        self.species.current_reference_name, "ensembl", release
                    )
                    infer_historical = infer_historical and not (
                        os.path.isdir(existing) or self.server.rstrip("/") != ENSEMBL_PLATFORM_FTP_SERVER
                    )
                metadata, standard = dated_assembly(
                    release, self.species,
                    assembly, assembly_accession, provider,
                    infer_historical=infer_historical,
                )
                if standard and include_alt is not None:
                    standard = include_alt == (self.species.current_reference_name == "GRCh38")
                selected_reference = (
                    self.species.current_reference_name if standard
                    else metadata["assembly_accession"]
                )
                selected_coverage = (
                    selected_reference == "GRCh38" if standard else False
                ) if include_alt is None else include_alt
                coverage = "including-alt" if selected_coverage else "primary"
                annotation_name = "ensembl" if standard else "ensembl-%s-%s-%s-" % (
                    metadata["provider"], metadata["assembly_accession"], coverage,
                )
            else:
                if provider is not None or include_alt is not None:
                    raise ValueError("provider/include_alt apply only to dated releases")
                metadata = numbered_assembly(self.species, release)
                validate_assembly(metadata, assembly, assembly_accession)
                selected_reference = self.species.which_reference(release)
                selected_coverage = None
                standard, annotation_name = True, "ensembl"
            self._resolved = dict(
                release=release, metadata=metadata, standard=standard,
                reference_name=reference_name or selected_reference,
                annotation_name=annotation_name, include_alt=selected_coverage,
            )
        self.include_alt = self._resolved["include_alt"]
        self.annotation_date = release if dated else None
        self.reference_name = self._resolved["reference_name"]
        genome_fasta = _genome_fasta_option(
            genome_fasta, download_genome_fasta, genome_fasta_path
        )
        self._genome_fasta_option = genome_fasta
        self.genome_fasta_type = genome_fasta_type
        self.genome_fasta_mask = genome_fasta_mask
        if dated:
            urls = self._dated_release_urls()
            self.download_url = urls.directory
            genome_fasta_url = urls.genome_fasta
            self.gtf_url = urls.gtf
            self.transcript_fasta_urls = [urls.cdna]
            self.protein_fasta_urls = [urls.pep]
        else:
            genome_fasta_url = make_genome_fasta_url(
                self.release, self.species, fasta_type=genome_fasta_type,
                mask=genome_fasta_mask, server=self.server,
            )
            self.gtf_url = make_gtf_url(
                ensembl_release=self.release, species=self.species, server=self.server
            )
            self.transcript_fasta_urls = [
                make_fasta_url(
                    ensembl_release=self.release,
                    species=self.species,
                    sequence_type=sequence_type,
                    server=self.server,
                )
                for sequence_type in ("cdna", "ncrna")
            ]
            self.protein_fasta_urls = [
                make_fasta_url(
                    ensembl_release=self.release,
                    species=self.species,
                    sequence_type="pep",
                    server=self.server,
                )
            ]
        self.genome_fasta_urls = [genome_fasta_url] if genome_fasta is True else []
        genome_fasta_source = genome_fasta_url if genome_fasta is True else genome_fasta

        Genome.__init__(
            self,
            reference_name=self.reference_name,
            annotation_name=self._resolved["annotation_name"],
            annotation_version=self.release,
            gtf_path_or_url=self.gtf_url,
            transcript_fasta_paths_or_urls=self.transcript_fasta_urls,
            protein_fasta_paths_or_urls=self.protein_fasta_urls,
            genome_fasta_path_or_url=genome_fasta_source,
            cache_directory_path=cache_directory_path,
        )

    @property
    def assembly(self):
        """Canonical build/coordinate-system identifier, or None if unknown."""
        return self._resolved["metadata"].get("assembly")

    @property
    def assembly_name(self):
        """Complete official assembly name, including any patch suffix."""
        return self._resolved["metadata"].get("assembly_name")

    @property
    def assembly_accession(self):
        """Versioned GCA/GCF accession, or None when upstream does not supply it."""
        return self._resolved["metadata"].get("assembly_accession")

    @property
    def provider(self):
        """Dated FTP annotation provider (not its human-readable API label)."""
        return self._resolved["metadata"].get("provider", "ensembl")

    @property
    def genome_id(self):
        """Pinned genome UUID for an explicit new-API selection, otherwise None."""
        return self._resolved["metadata"].get("genome_id")

    @property
    def website_release(self):
        """Website release label from API discovery, distinct from release."""
        return self._resolved["metadata"].get("website_release")

    def with_mirror(self, mirror, *, genome_fasta_path_or_url=None):
        """Return a Genome using this release's filenames at a flat mirror.

        Use the same mirror when installing and reopening data. New mirror
        downloads have their own source identities and do not populate the
        canonical EnsemblRelease cache. This is also the constructor used by
        the CLI's ``--custom-mirror`` option.

        Parameters
        ----------
        mirror : str or Path
            URL or directory containing the Ensembl filenames for this release.
        genome_fasta_path_or_url : str or Path, optional
            Explicit reference DNA source. Otherwise mirror requested Ensembl
            DNA, or preserve a local DNA file attached to this release.

        Returns
        -------
        Genome
            A source-explicit genome, usable across processes and in Varcode.
        """
        mirror = os.fspath(mirror)
        if not mirror:
            raise ValueError("mirror must be a nonempty URL or directory")

        def mirrored(source):
            return os.path.join(mirror, os.path.basename(source))

        dna = genome_fasta_path_or_url
        if dna is None:
            if self.genome_fasta_urls:
                dna = mirrored(self.genome_fasta_urls[0])
            elif isinstance(self._genome_fasta_option, str):
                dna = self._genome_fasta_option
        return Genome(
            reference_name=self.reference_name,
            annotation_name=self.annotation_name,
            annotation_version=self.release,
            gtf_path_or_url=mirrored(self.gtf_url),
            transcript_fasta_paths_or_urls=[mirrored(u) for u in self.transcript_fasta_urls],
            protein_fasta_paths_or_urls=[mirrored(u) for u in self.protein_fasta_urls],
            genome_fasta_path_or_url=dna,
        )

    def _dated_release_urls(self):
        """Files of the resolved dated assembly/provider selection."""
        if self._genome_fasta_option is True and self.genome_fasta_type != "toplevel":
            raise ValueError(
                "Dated releases publish only toplevel reference DNA, not %r"
                % (self.genome_fasta_type,)
            )
        return make_dated_release_urls(
            self.assembly_accession, self.provider, self.release,
            include_alt=self.include_alt,
            genome_fasta_mask=self.genome_fasta_mask, server=self.server,
        )

    def _make_genome_fasta(self, source):
        return release_genome_fasta(
            source,
            self.download_cache.cache_directory_path,
            install_string_function=self.genome_fasta_install_string,
        )

    def _explain_download_failure(self):
        if not is_dated_release(self.release):
            return Genome._explain_download_failure(self)
        return explain_unpublished_date(
            self.assembly_accession, self.provider, self.release, self.server,
            description=("%s %s" % (
                self.species.latin_name if self.species else "assembly", self.reference_name
            )) if self._resolved["standard"] else None,
        )

    def install_string(self):
        return self._install_command(only_genome_fasta=False)

    def genome_fasta_install_string(self):
        return self._install_command(only_genome_fasta=True)

    def _install_command(self, only_genome_fasta):
        command = "pyensembl install --release %s" % self.release
        if self.species is not None:
            command += " --species " + self.species.latin_name
        if not self._resolved["standard"]:
            command += " --assembly-accession %s --provider %s" % (
                self.assembly_accession, self.provider,
            )
            if self.include_alt:
                command += " --include-alt"
        expected_reference = (
            self.species.which_reference(self.release) if self._resolved["standard"]
            else self.assembly_accession
        )
        if self.reference_name != expected_reference:
            command += " --cache-reference-name " + shlex.quote(self.reference_name)
        default_server = ENSEMBL_PLATFORM_FTP_SERVER if is_dated_release(self.release) else (
            _default_server_for_species(self.species)
        )
        if self.server != default_server:
            command += " --server " + shlex.quote(self.server)
        if self.cache_directory_path is not None:
            command += " --cache-directory " + shlex.quote(os.fspath(self.cache_directory_path))
        if only_genome_fasta:
            command += " --only-genome-fasta"
        if isinstance(self._genome_fasta_option, str):
            command += " --genome-fasta-path " + shlex.quote(self._genome_fasta_option)
        elif self.requires_genome_fasta:
            if not only_genome_fasta:
                command += " --with-genome-fasta"
            if self.genome_fasta_type != "toplevel":
                command += " --genome-fasta-type " + self.genome_fasta_type
            if self.genome_fasta_mask != "none":
                command += " --masked " + self.genome_fasta_mask
        return command

    def __str__(self):
        selected = (
            "species=%r" % self.species.latin_name if self._resolved["standard"] and self.species
            else "assembly_accession=%r, provider=%r" % (self.assembly_accession, self.provider)
        )
        return "EnsemblRelease(release=%r, %s)" % (self.release, selected)

    def __eq__(self, other):
        # Attached reference DNA does not change annotation identity.
        return type(self) is type(other) and self._annotation_identity() == other._annotation_identity()

    def __hash__(self):
        return hash(self._annotation_identity())

    def _annotation_identity(self):
        if self._resolved["standard"]:
            return self.release, self.species
        return self._fields()

    def _genome_fasta_setup_hint(self):
        """Explain how to enable reference DNA recorded for this release."""

        def call(genome_fasta, **options):
            arguments = [repr(self.release)]
            if not self._resolved["standard"]:
                arguments.append("assembly_accession=%r" % self.assembly_accession)
                arguments.append("provider=%r" % self.provider)
                if self.include_alt:
                    arguments.append("include_alt=True")
            expected_reference = (
                self.species.which_reference(self.release) if self._resolved["standard"]
                else self.assembly_accession
            )
            if self.reference_name != expected_reference:
                arguments.append("reference_name=%r" % self.reference_name)
            default_server = ENSEMBL_PLATFORM_FTP_SERVER if is_dated_release(self.release) else (
                _default_server_for_species(self.species)
            )
            if self.server != default_server:
                arguments.append("server=%r" % self.server)
            if self.species is not None and self.species != human:
                arguments.append("species=%r" % self.species.latin_name)
            if self.cache_directory_path is not None:
                arguments.append("cache_directory_path=%r" % self.cache_directory_path)
            arguments.append("genome_fasta=%s" % genome_fasta)
            arguments += ["%s=%r" % option for option in options.items()]
            return "EnsemblRelease(%s)" % ", ".join(arguments)

        try:
            installed = GenomeFasta.installed_source(self.download_cache.cache_directory_path)
            if installed is not None and installed.installed_path is None:
                installed = None  # Recorded, but since moved or deleted.
        except (OSError, ValueError):
            installed = None  # The hint must not replace the original error.
        if installed is not None and not installed.remote:
            return "Reference DNA is attached to this release; use %s." % call(
                repr(installed.source)
            )
        if installed is None:
            options = None
        elif is_dated_release(self.release):
            mask = dated_genome_fasta_mask(installed.source)
            options = None if mask is None else ("toplevel", mask)
        else:
            options = canonical_source_options(installed.source)
        if options is not None:
            fasta_type, mask = options
            extra = {}
            if fasta_type != "toplevel":
                extra["genome_fasta_type"] = fasta_type
            if mask != "none":
                extra["genome_fasta_mask"] = mask
            return "Reference DNA is installed for this release; use %s." % call(
                "True", **extra
            )
        return "Use %s for Ensembl DNA, or genome_fasta=<path> for a local FASTA." % (
            call("True")
        )

    def to_dict(self):
        return {
            "release": self.release, "species": self.species, "server": self.server,
            "genome_fasta": self._genome_fasta_option,
            "genome_fasta_type": self.genome_fasta_type,
            "genome_fasta_mask": self.genome_fasta_mask,
            "cache_directory_path": self.cache_directory_path,
            "_resolved": self._resolved,
        }

    @classmethod
    def from_dict(cls, state_dict):
        """
        Deserialize EnsemblRelease without creating duplicate instances.
        """
        state_dict = dict(state_dict)
        # Pickles and JSON from 2.11.0 used two keys for reference DNA.
        genome_fasta_path = state_dict.pop("genome_fasta_path", None)
        download_genome_fasta = state_dict.pop("download_genome_fasta", False)
        if "genome_fasta" not in state_dict:
            state_dict["genome_fasta"] = genome_fasta_path or bool(download_genome_fasta)
        if "_resolved" not in state_dict:
            state_dict["_preserve_legacy"] = True
        return cls.cached(**state_dict)

    @classmethod
    def from_ensembl(cls, genome_id_or_accession, **options):
        """Resolve a genome using the new Ensembl metadata API (network access).

        The returned object pins its annotation date and assembly metadata.
        Ordinary construction and loading serialized objects remain offline.
        """
        from .ensembl_metadata import resolve_genome

        metadata, date = resolve_genome(genome_id_or_accession)
        accession = metadata["assembly_accession"]
        for name in ("release", "assembly_accession", "_resolved"):
            if name in options:
                raise TypeError("%s is selected by the metadata API" % name)
        validate_assembly(metadata, options.get("assembly"), accession)
        if options.get("provider", metadata["provider"]) != metadata["provider"]:
            raise ValueError("provider conflicts with the Ensembl API dataset")
        from .assembly import species_matches, species_for_record

        selected_species = options.get("species")
        if selected_species is not None:
            selected_species = check_species_object(selected_species)
            if not species_matches(selected_species, metadata):
                raise ValueError("species conflicts with the Ensembl API dataset")
        else:
            options["species"] = species_for_record(metadata)
        include_alt = options.get("include_alt", False)
        if not isinstance(include_alt, bool):
            raise TypeError("include_alt must be bool")
        coverage = "including-alt" if include_alt else "primary"
        state = dict(
            release=date, metadata=metadata, standard=False, include_alt=include_alt,
            reference_name=options.get("reference_name") or accession,
            annotation_name="ensembl-%s-%s-%s-" % (metadata["provider"], accession, coverage),
        )
        return cls(date, assembly_accession=accession, _resolved=state, **options)


def cached_release(release, species="human"):
    """
    Create an EnsemblRelease instance only if it's hasn't already been made,
    otherwise returns the old instance.
    Keeping this function for backwards compatibility but this functionality
    has been moving into the cached method of EnsemblRelease.
    """
    return EnsemblRelease.cached(release=release, species=species)

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
import os
import shlex
import warnings

from .genome import Genome
from .genome_fasta_cache import canonical_source_options, release_genome_fasta
from .genome_fasta import GenomeFasta
from .ensembl_versions import check_release_number, MAX_ENSEMBL_RELEASE
from .species import check_species_object, human

from .ensembl_url_templates import (
    ENSEMBL_FTP_SERVER,
    ENSEMBL_GENOMES_FTP_SERVER,
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
    """

    @classmethod
    def normalize_init_values(cls, release, species, server):
        """
        Normalizes the arguments which uniquely specify an EnsemblRelease
        genome.
        """
        release = check_release_number(release)
        species = check_species_object(species)
        if server is None or server == ENSEMBL_FTP_SERVER:
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
        cls, release=MAX_ENSEMBL_RELEASE, species=human, server=ENSEMBL_FTP_SERVER,
        *, genome_fasta=None, genome_fasta_type="toplevel", genome_fasta_mask="none",
        download_genome_fasta=None, genome_fasta_path=None,
    ):
        """
        Construct EnsemblRelease if it's never been made before, otherwise
        return an old instance.
        """
        init_args_tuple = cls.normalize_init_values(release, species, server)
        options = dict(
            genome_fasta=_genome_fasta_option(
                genome_fasta, download_genome_fasta, genome_fasta_path
            ),
            genome_fasta_type=genome_fasta_type,
            genome_fasta_mask=genome_fasta_mask,
        )
        cache_key = init_args_tuple + tuple(options.values())
        if cache_key in cls._genome_cache:
            genome = cls._genome_cache[cache_key]
        else:
            genome = cls._genome_cache[cache_key] = cls(*init_args_tuple, **options)
        return genome

    def __init__(
        self, release=MAX_ENSEMBL_RELEASE, species=human, server=ENSEMBL_FTP_SERVER,
        *, genome_fasta=None, genome_fasta_type="toplevel", genome_fasta_mask="none",
        download_genome_fasta=None, genome_fasta_path=None,
    ):
        """
        genome_fasta : True or path, optional
            Reference DNA for ``sequence()``: True for Ensembl's DNA for this
            release (see genome_fasta_type and genome_fasta_mask), or a local
            plain or gzip FASTA. Nothing is downloaded until
            ``download_genome_fasta()``, ``download()``, or ``pyensembl install``.
        """
        self.release, self.species, self.server = self.normalize_init_values(
            release=release, species=species, server=server
        )
        genome_fasta = _genome_fasta_option(
            genome_fasta, download_genome_fasta, genome_fasta_path
        )
        self._genome_fasta_option = genome_fasta
        self.genome_fasta_type = genome_fasta_type
        self.genome_fasta_mask = genome_fasta_mask
        genome_fasta_url = make_genome_fasta_url(
            self.release, self.species, fasta_type=genome_fasta_type,
            mask=genome_fasta_mask, server=self.server,
        )
        self.genome_fasta_urls = [genome_fasta_url] if genome_fasta is True else []
        genome_fasta_source = genome_fasta_url if genome_fasta is True else genome_fasta

        self.gtf_url = make_gtf_url(
            ensembl_release=self.release, species=self.species, server=self.server
        )

        self.transcript_fasta_urls = [
            make_fasta_url(
                ensembl_release=self.release,
                species=self.species,
                sequence_type="cdna",
                server=self.server,
            ),
            make_fasta_url(
                ensembl_release=self.release,
                species=self.species,
                sequence_type="ncrna",
                server=self.server,
            ),
        ]

        self.protein_fasta_urls = [
            make_fasta_url(
                ensembl_release=self.release,
                species=self.species,
                sequence_type="pep",
                server=self.server,
            )
        ]

        self.reference_name = self.species.which_reference(self.release)

        Genome.__init__(
            self,
            reference_name=self.reference_name,
            annotation_name="ensembl",
            annotation_version=self.release,
            gtf_path_or_url=self.gtf_url,
            transcript_fasta_paths_or_urls=self.transcript_fasta_urls,
            protein_fasta_paths_or_urls=self.protein_fasta_urls,
            genome_fasta_path_or_url=genome_fasta_source,
        )

    def _make_genome_fasta(self, source):
        return release_genome_fasta(
            source,
            self.download_cache.cache_directory_path,
            install_string_function=self.genome_fasta_install_string,
        )

    def install_string(self):
        return self._install_command(only_genome_fasta=False)

    def genome_fasta_install_string(self):
        return self._install_command(only_genome_fasta=True)

    def _install_command(self, only_genome_fasta):
        command = "pyensembl install --release %d --species %s" % (
            self.release,
            self.species.latin_name,
        )
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
        return "EnsemblRelease(release=%d, species='%s')" % (
            self.release,
            self.species.latin_name,
        )

    def __eq__(self, other):
        # Attached reference DNA does not change annotation identity.
        return (
            type(self) is type(other)
            and self.release == other.release
            and self.species == other.species
        )

    def __hash__(self):
        return hash((self.release, self.species))

    def _genome_fasta_setup_hint(self):
        """Explain how to enable reference DNA recorded for this release."""

        def call(genome_fasta, **options):
            arguments = [str(self.release)]
            if self.species != human:
                arguments.append("species=%r" % self.species.latin_name)
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
        options = canonical_source_options(installed.source) if installed else None
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
        return cls.cached(**state_dict)


def cached_release(release, species="human"):
    """
    Create an EnsemblRelease instance only if it's hasn't already been made,
    otherwise returns the old instance.
    Keeping this function for backwards compatibility but this functionality
    has been moving into the cached method of EnsemblRelease.
    """
    return EnsemblRelease.cached(release=release, species=species)

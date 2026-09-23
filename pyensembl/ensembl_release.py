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

from .genome import Genome
from .genome_fasta_cache import SharedGenomeFasta, is_canonical_source
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
        *, download_genome_fasta=False, genome_fasta_path=None,
        genome_fasta_type="toplevel", genome_fasta_mask="none",
    ):
        """
        Construct EnsemblRelease if it's never been made before, otherwise
        return an old instance.
        """
        init_args_tuple = cls.normalize_init_values(release, species, server)
        if genome_fasta_path is not None:
            if "://" in os.fspath(genome_fasta_path):
                raise ValueError("genome_fasta_path must be a local path; use Genome for custom URLs")
            genome_fasta_path = os.path.abspath(os.fspath(genome_fasta_path))
        options = dict(
            download_genome_fasta=download_genome_fasta,
            genome_fasta_path=genome_fasta_path,
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
        *, download_genome_fasta=False, genome_fasta_path=None,
        genome_fasta_type="toplevel", genome_fasta_mask="none",
    ):
        self.release, self.species, self.server = self.normalize_init_values(
            release=release, species=species, server=server
        )
        self._download_genome_fasta = download_genome_fasta
        if genome_fasta_path is not None:
            genome_fasta_path = os.fspath(genome_fasta_path)
            if "://" in genome_fasta_path:
                raise ValueError("genome_fasta_path must be a local path; use Genome for custom URLs")
            genome_fasta_path = os.path.abspath(genome_fasta_path)
        self._local_genome_fasta_path = genome_fasta_path
        self.genome_fasta_type = genome_fasta_type
        self.genome_fasta_mask = genome_fasta_mask
        genome_fasta_url = make_genome_fasta_url(
            self.release, self.species, fasta_type=genome_fasta_type,
            mask=genome_fasta_mask, server=self.server,
        )
        self.genome_fasta_urls = [genome_fasta_url] if download_genome_fasta and genome_fasta_path is None else []
        genome_fasta_source = genome_fasta_path or (genome_fasta_url if download_genome_fasta else None)

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
        if self.genome_fasta_urls and is_canonical_source(genome_fasta_source):
            self._genome_fasta = SharedGenomeFasta(
                genome_fasta_source, self.download_cache.cache_directory_path
            )

    def install_string(self):
        command = "pyensembl install --release %d --species %s" % (
            self.release,
            self.species.latin_name,
        )
        if self._local_genome_fasta_path:
            command += " --genome-fasta-path " + shlex.quote(self._local_genome_fasta_path)
        elif self.requires_genome_fasta:
            command += " --with-genome-fasta --genome-fasta-type %s --masked %s" % (
                self.genome_fasta_type, self.genome_fasta_mask,
            )
        return command

    def __str__(self):
        return "EnsemblRelease(release=%d, species='%s')" % (
            self.release,
            self.species.latin_name,
        )

    def __eq__(self, other):
        return (
            type(self) is type(other)
            and self.release == other.release
            and self.species == other.species
            and self._genome_fasta_path_or_url == other._genome_fasta_path_or_url
        )

    def __hash__(self):
        return hash((self.release, self.species, self._genome_fasta_path_or_url))

    def to_dict(self):
        return {
            "release": self.release, "species": self.species, "server": self.server,
            "download_genome_fasta": self._download_genome_fasta,
            "genome_fasta_path": self._local_genome_fasta_path,
            "genome_fasta_type": self.genome_fasta_type,
            "genome_fasta_mask": self.genome_fasta_mask,
        }

    @classmethod
    def from_dict(cls, state_dict):
        """
        Deserialize EnsemblRelease without creating duplicate instances.
        """
        return cls.cached(**state_dict)


def cached_release(release, species="human"):
    """
    Create an EnsemblRelease instance only if it's hasn't already been made,
    otherwise returns the old instance.
    Keeping this function for backwards compatibility but this functionality
    has been moving into the cached method of EnsemblRelease.
    """
    return EnsemblRelease.cached(release=release, species=species)

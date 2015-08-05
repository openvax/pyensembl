# Copyright (c) 2015. Mount Sinai School of Medicine
#
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
to be specific to (a particular release of) Ensembl.
"""

from .genome import Genome
from .gtf import GTF
from .ensembl_release_source import EnsemblReleaseSource
from .release_info import check_release_number, MAX_ENSEMBL_RELEASE,
from .species import species_from_name, human, Species
from .sequence_data import SequenceData
from .url_templates import ENSEMBL_FTP_SERVER, make_gtf_url, make_fasta_url

class EnsemblRelease(Genome):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database.
    """
    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 species=human,
                 server=ENSEMBL_FTP_SERVER,
                 auto_download=False):
        self.release = check_release_number(release)
        if isinstance(species, Species):
            self.species = species
        elif isinstance(species, str):
            self.species = species_from_name(species)
        else:
            raise ValueError("Unexpected type for species: %s : %s" % (
                species, type(species)))
        self.server = server

        self.gtf_url = make_gtf_url(
            ensembl_release=self.release,
            species=species,
            server=server)
        self.transcript_fasta_url = make_fasta_url(
            ensembl_release=self.release,
            species=self.species.latin_name,
            sequence_type="cdna",
            server=server)
        self.protein_fasta_url = make_fasta_url(
            ensembl_release=self.release,
            species=self.species.latin_name,
            sequence_type="pep",
            server=server)

        self.reference_name = self.species.which_reference(self.release)

        Genome.__init__(self,
                        reference_name=self.reference_name,
                        gtf_path_or_url=self.gtf_url,
                        transcript_fasta_path_or_url=self.transcript_fasta_url,
                        protein_fasta_path_or_url=self.protein_fasta_url,
                        annotation_name="Ensembl",
                        annotation_version=self.release,
                        species_name=self.species.latin_name,
                        require_ensembl_ids=True,
                        auto_download=auto_download)

    def install_string_console(self):
        return "pyensembl install --release %d" % self.release

    def install_string_python(self):
        return "EnsemblRelease(%d).install()" % self.release

    def build_gtf(self):
        return GTF(gtf_source=EnsemblReleaseSource(url=self.gtf_url,
                                                   release=self.release,
                                                   file_type="gtf",
                                                   reference_name=self.reference_name),
                   auto_download=self.auto_download)

    def build_transcript_sequences(self):
        transcript_fasta_source = EnsemblReleaseSource(
            url=self.transcript_fasta_url,
            release=self.release,
            file_type="fa",
            reference_name=self.reference_name)
        return SequenceData(
            fasta_source=transcript_fasta_source,
            require_ensembl_ids=self.require_ensembl_ids,
            auto_download=self.auto_download)

    def build_protein_sequences(self):
        protein_fasta_source = EnsemblReleaseSource(
            url=self.protein_fasta_url,
            release=self.release,
            file_type="fa",
            reference_name=self.reference_name)
        return SequenceData(
            fasta_source=protein_fasta_source,
            require_ensembl_ids=self.require_ensembl_ids,
            auto_download=self.auto_download)

    def __str__(self):
        return "EnsemblRelease(release=%d, species=%s)" % (
            self.release,
            self.species)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is EnsemblRelease and
            self.release == other.release and
            self.species == other.species)

    def __hash__(self):
        return hash((self.release, self.species))

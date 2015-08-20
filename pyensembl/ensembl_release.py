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
from .ensembl_release_versions import check_release_number, MAX_ENSEMBL_RELEASE
from .species import check_species_object, human

from .ensembl_url_templates import (
    ENSEMBL_FTP_SERVER,
    make_gtf_url,
    make_fasta_url
)

class EnsemblRelease(Genome):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database.
    """
    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 species=human,
                 server=ENSEMBL_FTP_SERVER):
        self.release = check_release_number(release)
        self.species = check_species_object(species)
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

        Genome.__init__(
            self,
            reference_name=self.reference_name,
            annotation_name="ensembl",
            annotation_version=self.release,
            gtf_path_or_url=self.gtf_url,
            transcript_fasta_path_or_url=self.transcript_fasta_url,
            protein_fasta_path_or_url=self.protein_fasta_url,
            require_ensembl_ids=True)

    def install_string(self):
        return "pyensembl install --release %d --species %s" % (
            self.release,
            self.species.latin_name)

    def __str__(self):
        return "EnsemblRelease(release=%d, species='%s')" % (
            self.release,
            self.species.latin_name)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is EnsemblRelease and
            self.release == other.release and
            self.species == other.species)

    def __hash__(self):
        return hash((self.release, self.species))

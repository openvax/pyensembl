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
from .release_info import (
    check_release_number,
    MAX_ENSEMBL_RELEASE,
    which_human_reference_name
)
from .sequence_data import SequenceData
from .url_templates import ENSEMBL_FTP_SERVER, make_gtf_url, make_fasta_url


class EnsemblRelease(Genome):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database.
    """
    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 species="homo_sapiens",
                 server=ENSEMBL_FTP_SERVER,
                 reference_name=None,
                 auto_download=False):
        self.release = check_release_number(release)
        self.species = species
        self.server = server

        self.gtf_url = make_gtf_url(
            ensembl_release=release,
            species=species,
            server=server)
        self.transcript_fasta_url = make_fasta_url(
            ensembl_release=release,
            species=species,
            sequence_type="cdna",
            server=server)
        self.protein_fasta_url = make_fasta_url(
            ensembl_release=release,
            species=species,
            sequence_type="pep",
            server=server)
        only_human = species == "homo_sapiens"

        if not reference_name:
            if only_human:
                reference_name = which_human_reference_name(self.release)
            else:
                raise ValueError("Must provide a reference_name for "
                                 "non-human Ensembl releases.")
        self.reference_name = reference_name

        Genome.__init__(self,
                        reference_name=self.reference_name,
                        gtf_path_or_url=self.gtf_url,
                        transcript_fasta_path_or_url=self.transcript_fasta_url,
                        protein_fasta_path_or_url=self.protein_fasta_url,
                        name="Ensembl",
                        version=self.release,
                        only_human=only_human,
                        auto_download=auto_download,
                        require_ensembl_ids=True)

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

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

from os.path import split

from .genome_source import GenomeSource
from .release_info import MAX_ENSEMBL_RELEASE
from .url_templates import ENSEMBL_FTP_SERVER, make_gtf_url, make_fasta_url


class EnsemblReleaseSource(GenomeSource):
    """
    Represents the source (URLs) of an Ensembl release. The superclass
    represents the source of any genomic database.
    """
    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 species="homo_sapiens",
                 server=ENSEMBL_FTP_SERVER):
        self.release = release
        self.species = species
        self.server = server

        gtf_url = make_gtf_url(
            ensembl_release=release,
            species=species,
            server=server)
        remote_gtf_filename = split(gtf_url)[1]
        assert remote_gtf_filename.endswith(".gtf.gz"), \
            "Expected remote GTF file %s to end with '.gtf.gz'" % (
                remote_gtf_filename)

        transcript_fasta_url = make_fasta_url(
            ensembl_release=release,
            species=species,
            sequence_type="cdna",
            server=server)
        protein_fasta_url = make_fasta_url(
            ensembl_release=release,
            species=species,
            sequence_type="pep",
            server=server)
        GenomeSource.__init__(self,
                              gtf_path=gtf_url,
                              transcript_fasta_path=transcript_fasta_url,
                              protein_fasta_path=protein_fasta_url)

    def install_string_console(self):
        return "pyensembl install --release %d" % self.release

    def install_string_python(self):
        return "EnsemblRelease(%d).install()" % self.release

    def __str__(self):
        return "EnsemblReleaseSource(%s)" % self._arg_list_str()

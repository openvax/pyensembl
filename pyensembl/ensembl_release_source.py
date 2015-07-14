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


class EnsemblReleaseSource(GenomeSource):
    """
    Represents the source (URLs) of an Ensembl release. The superclass
    represents the source of any genomic database.
    """
    def __init__(self,
                 release=MAX_ENSEMBL_RELEASE,
                 species="homo_sapiens",
                 server=ENSEMBL_FTP_SERVER):
        gtf_url = gtf_url(
            ensembl_release=release,
            species=species,
            server=server)
        transcript_fasta_url = fasta_url(
            ensembl_release=self.release,
            species=self.species,
            sequence_type="cdna",
            server=self.server)
        protein_fasta_url = fasta_url(
            ensembl_release=self.release,
            species=self.species,
            sequence_type="pep",
            server=self.server)
        GenomeSource.__init__(gtf_url=gtf_url,
                              transcript_fasta_url=transcript_fasta_url,
                              protein_fasta_url=protein_fasta_url)

    def install_string_console(self):
        return "pyensembl install --release %d" % self.release

    def install_string_python(self):
        return "EnsemblRelease(%d).install()" % self.release

    def __str__(self):
        return "EnsemblReleaseSource(%s)" % self._arg_list_str()

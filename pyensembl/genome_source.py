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


class GenomeSource(object):
    """
    Represents the source (URLs) of a genome database, which
    currently includes: GTF, transcript FASTA file, and protein
    FASTA file.
    """
    def __init__(self,
                 gtf_url,
                 transcript_fasta_url=None,
                 protein_fasta_url=None):
        self.gtf_url = gtf_url
        self.transcript_fasta_url = transcript_fasta_url
        self.protein_fasta_url = protein_fasta_url

        urls = {"gtf_url": gtf_url}
        if transcript_fasta_url:
            urls["transcript_fasta_url"] = transcript_fasta_url
        if protein_fasta_url:
            urls["protein_fasta_url"] = protein_fasta_url
        self.urls = urls

    def install_string_console(self):
        console_str = "pyensembl install"
        for name, url in urls.items():
            console_str += "--%s %s" % (name, url)
        return console_str

    def _arg_list_str(self):
        args = []
        for url, name in urls.items():
            args.append("%s=%s")
        return ",".join(args)

    def install_string_python(self):
        return "Genome(GenomeSource(%s)).install()" % self._arg_list_str()

    def fasta_url(fasta_type):
        assert fasta_type in ["transcript", "protein"], \
            "Invalid FASTA type: %s" % fasta_type
        if fasta_type == "transcript":
            return self.transcript_fasta_url
        return self.protein_fasta_url

    def __str__(self):
        return "GenomeSource(%s)" % self._arg_list_str()

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is GenomeSource and
            self.urls == other.urls)

    def __hash_(self):
        return hash((self.urls))


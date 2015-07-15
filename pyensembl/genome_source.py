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
    Represents the source (URLs or local file paths) of a genome
    database, which currently includes: GTF, transcript FASTA file, 
    and protein FASTA file.
    """
    def __init__(self,
                 gtf_path,
                 transcript_fasta_path=None,
                 protein_fasta_path=None):
        self.gtf_path = gtf_path
        self.transcript_fasta_path = transcript_fasta_path
        self.protein_fasta_path = protein_fasta_path

        paths = {"gtf_path": gtf_path}
        if transcript_fasta_path:
            paths["transcript_fasta_path"] = transcript_fasta_path
        if protein_fasta_path:
            paths["protein_fasta_path"] = protein_fasta_path
        self.paths = paths

    def install_string_console(self):
        console_str = "pyensembl install"
        for name, path in self.paths.items():
            console_str += "--%s %s" % (name, path)
        return console_str

    def _arg_list_str(self):
        args = []
        for name, path in self.paths.items():
            args.append("%s=%s" % (name, path))
        return ",".join(args)

    def install_string_python(self):
        return "Genome(GenomeSource(%s)).install()" % self._arg_list_str()

    def fasta_path(self, fasta_type):
        assert fasta_type in ["transcript", "protein"], \
            "Invalid FASTA type: %s" % fasta_type
        if fasta_type == "transcript":
            return self.transcript_fasta_path
        return self.protein_fasta_path

    def __str__(self):
        return "GenomeSource(%s)" % self._arg_list_str()

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is GenomeSource and
            self.paths == other.paths)

    def __hash_(self):
        return hash((self.paths))

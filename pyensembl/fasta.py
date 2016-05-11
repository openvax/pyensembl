# Copyright (c) 2015-2016. Mount Sinai School of Medicine
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
The worse sin in bioinformatics is to write your own FASTA parser.
Unfortunately, small errors creep in to different FASTA files on the
Ensembl FTP server that no proper FASTA parser lets you skip over.
"""

from __future__ import print_function, division, absolute_import
from gzip import GzipFile
import logging

from six import binary_type

def _parse_header_id(line):
    """
    Pull the transcript or protein identifier from the header line
    which starts with '>'
    """
    if type(line) is not binary_type:
        raise TypeError("Expected header line to be of type %s but got %s" % (
            binary_type, type(line)))

    if len(line) <= 1:
        raise ValueError("No identifier on FASTA line")

    # split line at first space to get the unique identifier for
    # this sequence
    space_index = line.find(b" ")
    if space_index >= 0:
        identifier = line[1:space_index]
    else:
        identifier = line[1:]

    # annoyingly Ensembl83 reformatted the transcript IDs of its
    # cDNA FASTA to include sequence version numbers
    # .e.g.
    # "ENST00000448914.1" instead of "ENST00000448914"
    # So now we have to parse out the identifier
    dot_index = identifier.find(b".")
    if dot_index >= 0:
        identifier = identifier[:dot_index]

    return identifier.decode("ascii")

class FastaParser(object):
    """
    FastaParser object consumes lines of a FASTA file incrementally
    while building up a dictionary mapping sequence identifiers to sequences.
    """
    def __init__(self, sequence_type=binary_type):
        self.sequence_type = sequence_type
        self.current_id = None
        self.current_lines = []
        self.fasta_dictionary = {}

    def read_file(self, fasta_path):
        if fasta_path.endswith("gz") or fasta_path.endswith("gzip"):
            f = GzipFile(fasta_path, 'rb')
        else:
            f = open(fasta_path, 'rb')

        for line in f:
            line = line.rstrip()

            if len(line) == 0:
                continue

            # have to slice into a bytes object or else I get a single integer
            first_char = line[0:1]

            if first_char == b">":
                self._read_header(line)
            elif first_char == b";":
                # semicolon are comment characters
                continue
            else:
                self.current_lines.append(line)
        self._end_of_file()
        return self.fasta_dictionary

    def _end_of_file(self):
        self._add_current_sequence_to_dictionary()

    def _add_current_sequence_to_dictionary(self):
        # when we hit a new entry, if this isn't the first
        # entry of the file then put the last one in the dictionary
        if self.current_id:
            if len(self.current_lines) == 0:
                logging.warn("No sequence data for '%s'" % self.current_id)
            else:
                self.fasta_dictionary[
                    self.current_id] = self.sequence_type(
                        b"".join(self.current_lines))

    def _read_header(self, line):
        self._add_current_sequence_to_dictionary()

        self.current_id = _parse_header_id(line)

        if len(self.current_id) == 0:
            logging.warn("Unable to parse ID from header line: %s" % line)

        self.current_lines = []

def parse_fasta_dictionary(fasta_path, sequence_type=binary_type):
    """
    Given a path to a FASTA (or compressed FASTA) file, returns a dictionary
    mapping its sequence identifiers to sequences.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    sequence_type : type
        Default is str.

    Returns dictionary from string identifiers to sequences of type bytes.
    """
    parser = FastaParser(sequence_type=sequence_type)
    return parser.read_file(fasta_path)

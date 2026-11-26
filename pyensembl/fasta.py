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

from contextlib import nullcontext
from gzip import GzipFile
import io
import logging
import os
from os.path import basename

from datacache.fasta import (
    first_fasta_id as datacache_first_fasta_id,
    iter_fasta_records,
    parse_fasta_header_id,
)

from .common import read_progress

logger = logging.getLogger(__name__)


def _parse_header_id(line):
    """
    Pull the transcript or protein identifier from the header line
    which starts with '>'.

    The full versioned form (e.g. ``ENSP00000123456.3``) is returned when
    the header carries a version, preserving it as the authoritative
    identity of the sequence. Lookups match bare and versioned IDs against
    it (see :meth:`pyensembl.SequenceData.stored_id`).

    Non-ENS IDs (e.g. TAIR ``AT1G01010.1``, where ``.1`` is an isoform
    suffix rather than a version) are returned verbatim — this function
    does no surgery on them.
    """
    if type(line) is not bytes:
        raise TypeError(
            "Expected header line to be of type %s but got %s" % (bytes, type(line))
        )
    return parse_fasta_header_id(line[1:].decode("ascii"), split_pipe=True)


class FastaParser(object):
    """
    FastaParser object consumes lines of a FASTA file incrementally
    while building up a dictionary mapping sequence identifiers to sequences.
    """

    def __init__(self):
        self.current_id = None
        self.current_lines = []

    def read_file(self, fasta_path, show_progress=False):
        """
        Read the contents of a FASTA file into a dictionary
        """
        fasta_dictionary = {}
        for identifier, sequence in self.iterate_over_file(fasta_path, show_progress):
            fasta_dictionary[identifier] = sequence
        return fasta_dictionary

    def iterate_over_file(self, fasta_path, show_progress=False):
        """
        Generator that yields identifiers paired with sequences.
        """
        description = "Reading " + basename(fasta_path)
        with (
            open(fasta_path, "rb") as raw,
            read_progress(raw, description, show_progress) as report_progress,
            self._open(fasta_path, raw) as binary_stream,
            io.TextIOWrapper(binary_stream, encoding="ascii") as text_stream,
        ):
            for record in iter_fasta_records(
                text_stream,
                source_name=str(fasta_path),
                ignore_comments=True,
                allow_empty_header=True,
                allow_pre_header=True,
                on_header=report_progress,
            ):
                identifier = parse_fasta_header_id(
                    record.header.rstrip(), split_pipe=True
                )
                if not identifier:
                    logger.warning(
                        "Unable to parse ID from header line: >%s",
                        record.header,
                    )
                    continue
                if not record.sequence:
                    logger.warning("No sequence data for '%s'", identifier)
                    continue
                yield identifier, record.sequence

    def _open(self, fasta_path, raw):
        """
        Read an open text or gzip-compressed file as a stream of bytes.
        """
        lower = os.fspath(fasta_path).lower()
        if lower.endswith("gz") or lower.endswith("gzip"):
            return GzipFile(fileobj=raw, mode="rb")
        return nullcontext(raw)


def parse_fasta_dictionary(fasta_path, show_progress=False):
    """
    Given a path to a FASTA (or compressed FASTA) file, returns a dictionary
    mapping its sequence identifiers to sequences.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    show_progress : bool, optional
        Display a progress bar while reading.

    Returns dictionary from string identifiers to string sequences.
    """
    parser = FastaParser()
    return parser.read_file(fasta_path, show_progress=show_progress)


def first_fasta_id(fasta_path):
    """Identifier of the first record in a FASTA file, or None if it has none."""
    return datacache_first_fasta_id(
        fasta_path, split_pipe=True, encoding="ascii"
    )

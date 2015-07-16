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

from __future__ import print_function, division, absolute_import

from memoized_property import memoized_property

from .biotypes import is_valid_biotype
from .common import memoize
from .locus import Locus


class Transcript(Locus):
    """
    Transcript encompasses the locus, exons, and sequence of an Ensembl
    transcript.

    Lazily fetches sequence in case we"re constructing many Transcripts
    and not using the sequence, avoid the memory/performance overhead
    of fetching and storing sequences from a FASTA file.
    """
    def __init__(
            self,
            transcript_id,
            transcript_name,
            contig,
            start,
            end,
            strand,
            biotype,
            gene_id,
            ensembl,
            require_valid_biotype=True):
        if require_valid_biotype and not is_valid_biotype(biotype):
            raise ValueError(
                "Invalid biotype '%s' for transcript with ID=%s, name=%s" % (
                    biotype, transcript_id, transcript_name))

        Locus.__init__(self, contig, start, end, strand)
        self.id = transcript_id
        self.name = transcript_name
        self.ensembl = ensembl
        self.db = ensembl.db
        self.biotype = biotype
        self.gene_id = gene_id

    def __str__(self):
        return (
            "Transcript(id=%s,"
            " name=%s,"
            " gene_id=%s,"
            " gene_name=%s,"
            " biotype=%s,"
            " location=%s:%d-%d)") % (
                    self.id,
                    self.name,
                    self.gene.id,
                    self.gene.name,
                    self.biotype,
                    self.contig,
                    self.start,
                    self.end)

    def __repr__(self):
        return str(self)

    def __len__(self):
        """
        Length of a transcript is the sum of its exon lengths
        """
        return sum(len(exon) for exon in self.exons)

    def __eq__(self, other):
        return (
            other.__class__ is Transcript and
            self.id == other.id and
            self.ensembl == other.ensembl)

    def __hash__(self):
        return hash(self.id)

    @memoized_property
    def gene(self):
        return self.ensembl.gene_by_id(self.gene_id)

    @memoized_property
    def exons(self):
        # need to look up exon_number alongside ID since each exon may
        # appear in multiple transcripts and have a different exon number
        # in each transcript
        columns = ["exon_number", "exon_id"]
        exon_numbers_and_ids = self.db.query(
            columns,
            filter_column="transcript_id",
            filter_value=self.id,
            feature="exon")

        # fill this list in its correct order (by exon_number) by using
        # the exon_number as a 1-based list offset
        exons = [None] * len(exon_numbers_and_ids)

        for exon_number, exon_id in exon_numbers_and_ids:
            exon = self.ensembl.exon_by_id(exon_id)
            exon_number = int(exon_number)
            assert exon_number >= 1, "Invalid exon number: %s" % exon_number
            assert exon_number <= len(exons), \
                "Invalid exon number: %s (max expected = %d)" % (
                    exon_number, len(exons))

            # exon_number is 1-based, convert to list index by subtracting 1
            exons[exon_number - 1] = exon

        assert all(exon is not None for exon in exons), \
            "Missing exons %s for transcript %s" % (
                [
                    i
                    for i, maybe_exon
                    in enumerate(exons) if maybe_exon is None
                ],
                self.id)
        return exons

    # possible annotations associated with transcripts
    _TRANSCRIPT_FEATURES = {"start_codon", "stop_codon", "UTR", "CDS"}

    @memoize
    def _transcript_feature_position_ranges(self, feature, required=True):
        """
        Find start/end chromosomal position range of features
        (such as start codon) for this transcript.
        """
        if feature not in self._TRANSCRIPT_FEATURES:
            raise ValueError("Invalid transcript feature: %s" % feature)

        results = self.db.query(
            select_column_names=["start", "end"],
            filter_column="transcript_id",
            filter_value=self.id,
            feature=feature)

        if required and len(results) == 0:
            raise ValueError(
                "Transcript %s does not contain feature %s" % (
                    self.id, feature))
        return results

    @memoize
    def _transcript_feature_positions(self, feature):
        """
        Get unique positions for feature, raise an error if feature is absent.
        """
        ranges = self._transcript_feature_position_ranges(
            feature, required=True)
        results = []
        # a feature (such as a stop codon), maybe be split over multiple
        # contiguous ranges. Collect all the nucleotide positions into a
        # single list.
        for (start, end) in ranges:
            # since Ensembl ranges are [inclusive, inclusive] and
            # Python ranges are [inclusive, exclusive) we have to increment
            # the end position
            for position in range(start, end + 1):
                assert position not in results, \
                    "Repeated position %d for %s" % (position, feature)
                results.append(position)
        return results

    @memoize
    def _codon_positions(self, feature):
        """
        Parameters
        ----------
        feature : str
            Possible values are "start_codon" or "stop_codon"

        Returns list of three chromosomal positions.
        """
        results = self._transcript_feature_positions(feature)
        if len(results) != 3:
            raise ValueError(
                "Expected 3 positions for %s of %s but got %d" % (
                    feature,
                    self.id,
                    len(results)))
        return results

    @memoized_property
    def contains_start_codon(self):
        """
        Does this transcript have an annotated start_codon entry in Ensembl?
        """
        start_codons = self._transcript_feature_position_ranges(
            "start_codon", required=False)
        return len(start_codons) > 0

    @memoized_property
    def contains_stop_codon(self):
        """
        Does this transcript have an annotated stop_codon entry in Ensembl?
        """
        stop_codons = self._transcript_feature_position_ranges(
            "stop_codon", required=False)
        return len(stop_codons) > 0

    @memoized_property
    def start_codon_positions(self):
        """
        Chromosomal positions of nucleotides in start codon.
        """
        return self._codon_positions("start_codon")

    @memoized_property
    def stop_codon_positions(self):
        """
        Chromosomal positions of nucleotides in stop codon.
        """
        return self._codon_positions("stop_codon")

    @memoized_property
    def exon_intervals(self):
        """List of (start,end) tuples for each exon of this transcript,
        in the order specified by the 'exon_number' column of the Ensembl
        exon table.
        """
        results = self.db.query(
            select_column_names=["exon_number", "start", "end"],
            filter_column="transcript_id",
            filter_value=self.id,
            feature="exon")
        sorted_intervals = [None] * len(results)
        for (exon_number, start, end) in results:
            sorted_intervals[int(exon_number) - 1] = (start, end)
        return sorted_intervals

    def spliced_offset(self, position):
        """
        Convert from an absolute chromosomal position to the offset into
        this transcript"s spliced mRNA.

        Position must be inside some exon (otherwise raise exception).
        """
        # this code is performance sensitive, so switching from
        # typechecks.require_integer to a simpler assertion
        assert type(position) == int, \
            "Position argument must be an integer, got %s : %s" % (
                position, type(position))

        if position < self.start or position > self.end:
            raise ValueError(
                "Invalid position: %d (must be between %d and %d)" % (
                    position,
                    self.start,
                    self.end))

        # offset from beginning of unspliced transcript (including introns)
        unspliced_offset = self.offset(position)
        total_spliced_offset = 0

        # traverse exons in order of their appearance on the strand
        # Since absolute positions may decrease if on the negative strand,
        # we instead use unspliced offsets to get always increasing indices.
        #
        # Example:
        #
        # Exon Name:                exon 1                exon 2
        # Spliced Offset:           123456                789...
        # Intron vs. Exon: ...iiiiiieeeeeeiiiiiiiiiiiiiiiieeeeeeiiiiiiiiiii...
        for exon in self.exons:
            exon_unspliced_start, exon_unspliced_end = self.offset_range(
                exon.start, exon.end)
            # If the relative position is not within this exon, keep a running
            # total of the total exonic length-so-far.
            #
            # Otherwise, if the relative position is within an exon, get its
            # offset into that exon by subtracting the exon"s relative start
            # position from the relative position. Add that to the total exonic
            # length-so-far.
            if exon_unspliced_start <= unspliced_offset <= exon_unspliced_end:
                # all offsets are base 0, can be used as indices into
                # sequence string
                exon_offset = unspliced_offset - exon_unspliced_start
                return total_spliced_offset + exon_offset
            else:
                exon_length = len(exon)  # exon_end_position - exon_start_position + 1
                total_spliced_offset += exon_length
        raise ValueError(
            "Couldn't find position %d on any exon of %s" % (
                position, self.id))

    @memoized_property
    def start_codon_unspliced_offsets(self):
        """
        Offsets from start of unspliced pre-mRNA transcript
        of nucleotides in start codon.
        """
        return [
            self.position_offset(position)
            for position
            in self.start_codon_positions
        ]

    @memoized_property
    def stop_codon_unspliced_offsets(self):
        """
        Offsets from start of unspliced pre-mRNA transcript
        of nucleotides in stop codon.
        """
        return [
            self.position_offset(position)
            for position
            in self.stop_codon_positions
        ]

    def _contiguous_offsets(self, offsets):
        """
        Sorts the input list of integer offsets,
        ensures that values are contiguous.
        """
        offsets.sort()
        for i in range(len(offsets) - 1):
            assert offsets[i] + 1 == offsets[i + 1], \
                "Offsets not contiguous: %s" % (offsets,)
        return offsets

    @memoized_property
    def start_codon_spliced_offsets(self):
        """
        Offsets from start of spliced mRNA transcript
        of nucleotides in start codon.
        """
        offsets = [
            self.spliced_offset(position)
            for position
            in self.start_codon_positions
        ]
        return self._contiguous_offsets(offsets)

    @memoized_property
    def stop_codon_spliced_offsets(self):
        """
        Offsets from start of spliced mRNA transcript
        of nucleotides in stop codon.
        """
        offsets = [
            self.spliced_offset(position)
            for position
            in self.stop_codon_positions
        ]
        return self._contiguous_offsets(offsets)

    @memoized_property
    def coding_sequence_position_ranges(self):
        """
        Return absolute chromosome position ranges for CDS fragments
        of this transcript
        """
        return self._transcript_feature_position_ranges("CDS")

    @memoized_property
    def complete(self):
        """
        Consider a transcript complete if it has start and stop codons and
        a coding sequence whose length is divisible by 3
        """
        return (
            self.contains_start_codon and
            self.contains_stop_codon and
            self.coding_sequence is not None and
            len(self.coding_sequence) % 3 == 0
        )

    @memoized_property
    def sequence(self):
        """
        Spliced cDNA sequence of transcript
        (includes 5" UTR, coding sequence, and 3" UTR)
        """
        return self.ensembl.transcript_sequences.get(self.id)

    @memoized_property
    def first_start_codon_spliced_offset(self):
        """
        Offset of first nucleotide in start codon into the spliced mRNA
        (excluding introns)
        """
        start_offsets = self.start_codon_spliced_offsets
        return min(start_offsets)

    @memoized_property
    def last_stop_codon_spliced_offset(self):
        """
        Offset of last nucleotide in stop codon into the spliced mRNA
        (excluding introns)
        """
        stop_offsets = self.stop_codon_spliced_offsets
        return max(stop_offsets)

    @memoized_property
    def coding_sequence(self):
        """
        cDNA coding sequence (from start codon to stop codon, without
        any introns)
        """
        start = self.first_start_codon_spliced_offset
        end = self.last_stop_codon_spliced_offset
        # If start codon is the at nucleotide offsets [3,4,5] and
        # stop codon is at nucleotide offsets  [20,21,22]
        # then start = 3 and end = 22.
        #
        # Adding 1 to end since Python uses non-inclusive ends in slices/ranges.

        # pylint: disable=invalid-slice-index
        # TODO(tavi) Figure out pylint is not happy with this slice
        return self.sequence[start:end + 1]

    @memoized_property
    def five_prime_utr_sequence(self):
        """
        cDNA sequence of 5' UTR
        (untranslated region at the beginning of the transcript)
        """
        # pylint: disable=invalid-slice-index
        # TODO(tavi) Figure out pylint is not happy with this slice
        return self.sequence[:self.first_start_codon_spliced_offset]

    @memoized_property
    def three_prime_utr_sequence(self):
        """
        cDNA sequence of 3' UTR
        (untranslated region at the end of the transcript)
        """
        return self.sequence[self.last_stop_codon_spliced_offset + 1:]

    @memoized_property
    def protein_id(self):
        result_tuple = self.db.query_one(
            select_column_names=["protein_id"],
            filter_column="transcript_id",
            filter_value=self.id,
            feature="CDS",
            distinct=True,
            required=False)
        if result_tuple:
            return result_tuple[0]
        else:
            return None

    @memoized_property
    def protein_sequence(self):
        if self.protein_id:
            return self.ensembl.protein_sequences.get(self.protein_id)
        else:
            return None

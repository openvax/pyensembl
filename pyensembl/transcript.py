from __future__ import print_function, division, absolute_import

from .biotypes import is_valid_biotype
from .exon import Exon
from .locus import Locus

from memoized_property import memoized_property
from .type_checks import require_integer, require_string

class Transcript(Locus):
    """
    Transcript encompasses the locus, exons, and sequence of an Ensembl
    transcript.

    Lazily fetches sequence in case we're constructing many Transcripts
    and not using the sequence, avoid the memory/performance overhead
    of fetching and storing sequences from a FASTA file.
    """
    def __init__(self, transcript_id, db, reference):
        """
        Parameters
        ----------
        transcript_id : str

        db : pyensembl.Database
            Database of annotations

        reference : ReferenceTranscripts
        """
        require_string(transcript_id, "transcript ID")

        self.id = transcript_id
        self.db = db
        self.reference = reference

        columns = [
            'transcript_name',
            'transcript_biotype',
            'seqname',
            'start',
            'end',
            'strand',
            'gene_name',
            'gene_id',
        ]
        name, biotype, contig, start, end, strand, gene_name, gene_id = \
            self.db.query_one(
                select_column_names=columns,
                filter_column='transcript_id',
                filter_value=transcript_id,
                feature='transcript',
                distinct=True)

        Locus.__init__(self, contig, start, end, strand)

        if not name:
            raise ValueError(
                "Missing name for transcript with ID=%s" % transcript_id)
        else:
            self.name = name

        if not biotype:
            raise ValueError(
                "Missing biotype for transcript with ID=%s, name=%s" % (
                    transcript_id, name))
        elif not is_valid_biotype(biotype):
            raise ValueError(
                "Invalid biotype '%s' for transcript with ID=%s, name=%s" % (
                    biotype, transcript_id, name))
        else:
            self.biotype = biotype

        if gene_name is None:
            raise ValueError(
                "Missing gene name for transcript with ID=%s" % transcript_id)
        else:
            self.gene_name = gene_name

        if gene_id is None:
            raise ValueError(
                "Missing gene ID for transcript with ID=%s" % transcript_id)
        else:
            self.gene_id = gene_id


    def __str__(self):
        return "Transcript(id=%s, name=%s, gene_name=%s)" % (
                    self.id, self.name, self.gene_name)

    def __repr__(self):
        return str(self)

    def __len__(self):
        """
        Length of a transcript is the sum of its exon lengths
        """
        return sum(len(exon) for exon in self.exons)

    def __eq__(self, other):
        return (
            isinstance(other, Transcript) and
            self.id == other.id and
            self.db == other.db and
            self.reference == other.reference)

    def __hash__(self):
        return hash(self.id)

    @memoized_property
    def exons(self):
        columns = ['exon_number', 'exon_id']
        results = self.db.query(
            columns,
            filter_column='transcript_id',
            filter_value=self.id,
            feature='exon')

         # fill this list in its correct order (by exon_number) by using
         # the exon_number as a 1-based list offset
        exons = [None] * len(results)

        for exon_number, exon_id in results:
            exon = Exon(exon_id, self.db)
            exon_number = int(exon_number)
            assert exon_number >= 1, "Invalid exon number: %s" % exon_number
            assert exon_number <= len(exons), \
                "Invalid exon number: %s (max expected = %d)" % (
                    exon_number, len(exons))

            # exon_number is 1-based, convert to list index by subtracting 1
            exons[exon_number - 1] = exon

        assert all(exon is not None for exon in exons), \
            "Missing exons %s for transcript %s" % (
                [i for i, exon in enumerate(exons) if exon is None],
                self.id
            )
        return exons

    # possible annotations associated with transcripts
    _TRANSCRIPT_FEATURES = {'start_codon', 'stop_codon', 'UTR', 'CDS'}

    def _transcript_feature_position_ranges(self, feature, required=True):
        """
        Find start/end chromosomal position range of features
        (such as start codon) for this transcript.
        """

        if feature not in self._TRANSCRIPT_FEATURES:
            raise ValueError("Invalid transcript feature: %s" % feature)

        results = self.db.query(
            select_column_names=['start', 'end'],
            filter_column='transcript_id',
            filter_value=self.id,
            feature=feature)

        if required and len(results) == 0:
            raise ValueError(
                "Transcript %s does not contain feature %s" % (
                    self.id, feature))
        return results

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
            for position in range(start, end+1):
                assert position not in results, \
                    "Repeated position %d for %s" % (position, feature)
                results.append(position)
        return results

    def _codon_positions(self, feature):
        """
        Parameters
        ----------
        feature : str
            Possible values are 'start_codon' or 'stop_codon'

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
            'start_codon', required=False)
        return len(start_codons) > 0

    @memoized_property
    def contains_stop_codon(self):
        """
        Does this transcript have an annotated stop_codon entry in Ensembl?
        """
        stop_codons = self._transcript_feature_position_ranges(
            'stop_codon', required=False)
        return len(stop_codons) > 0

    @memoized_property
    def start_codon_positions(self):
        """
        Chromosomal positions of nucleotides in start codon.
        """
        return self._codon_positions('start_codon')

    @memoized_property
    def stop_codon_positions(self):
        """
        Chromosomal positions of nucleotides in stop codon.
        """
        return self._codon_positions('stop_codon')

    def spliced_offset(self, position):
        """
        Convert from an absolute chromosomal position to the offset into
        this transcript's spliced mRNA.

        Position must be inside some exon (otherwise raise exception).
        """
        require_integer(position, "position")
        if position < self.start or position > self.end:
            raise ValueError(
                "Invalid position: %d (must be between %d and %d)" % (
                    position,
                    self.start,
                    self.end))

        # offset from beginning of unspliced transcript (including introns)
        unspliced_offset = self.position_offset(position)
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
            # offset into that exon by subtracting the exon's relative start
            # position from the relative position. Add that to the total exonic
            # length-so-far.
            if exon_unspliced_start <= unspliced_offset <= exon_unspliced_end:
                # all offsets are base 0, can be used as indices into
                # sequence string
                exon_offset = unspliced_offset - exon_unspliced_start
                return total_spliced_offset + exon_offset
            else:
                total_spliced_offset += len(exon)
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
        for i in range(len(offsets)-1):
            assert offsets[i] + 1 == offsets[i+1], \
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
        return self._transcript_feature_position_ranges('CDS')

    @memoized_property
    def complete(self):
        """
        Consider a transcript complete if it has start and stop codons and
        a coding sequence whose length is divisible by 3
        """
        return (
            self.contains_start_codon and
            self.contains_stop_codon and
            len(self.coding_sequence) % 3 == 0
        )

    @memoized_property
    def sequence(self):
        """
        Spliced cDNA sequence of transcript
        (includes 5' UTR, coding sequence, and 3' UTR)
        """
        try:
            # fetch this transcript's Sequence from the attached
            # ReferenceTranscripts object
            return self.reference.transcript_sequence(self.id)
        except KeyError:
            raise ValueError(
                "No sequence for transcript %s in reference %s" % (
                    self.id, self.reference))

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
        return self.sequence[start:end+1]

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
        return self.sequence[self.last_stop_codon_spliced_offset+1:]


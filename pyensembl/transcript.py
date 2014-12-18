from locus import Locus
from exon import Exon

from pyfaidx import Sequence

class Transcript(Locus):
    """
    Transcript encompasses the locus, exons, and sequence of an Ensembl
    transcript.

    Lazily fetches sequence in case we're constructing many Transcripts
    and not using the sequence, avoid the memory/performance overhead
    of fetching and storing sequences from a FASTA file.
    """
    def __init__(self, transcript_id, db, reference):
        if not isinstance(transcript_id, (unicode, str)):
            raise TypeError(
                "Expected transcript ID to be string, got %s : %s" % (
                transcript_id, type(transcript_id)))

        self.id = transcript_id
        self.db = db
        self.reference = reference

        columns = [
            'transcript_name',
            'seqname',
            'start',
            'end',
            'strand',
            'gene_name',
            'gene_id'
        ]
        transcript_name, contig, start, end, strand, gene_name, gene_id = \
            self.db.query_one(
                select_column_names=columns,
                filter_column='transcript_id',
                filter_value=transcript_id,
                feature='transcript',
                distinct=True)

        Locus.__init__(self, contig, start, end, strand)

        if not transcript_name:
            raise ValueError(
                "Missing name for transcript with ID = %s" % transcript_name)
        self.name = transcript_name

        if gene_name is None:
            raise ValueError(
                "Missing gene name for transcript with ID = %s" % transcript_id)
        self.gene_name = gene_name

        if gene_id is None:
            raise ValueError(
                "Missing gene ID for transcript with ID = %s" % transcript_id)
        self.gene_id = gene_id

        self._sequence = None

    def __str__(self):
        return "Transcript(id=%s, name=%s, gene_name=%s)" % (
                    self.id, self.name, self.gene_name)

    def __repr__(self):
        return str(self)

    def __len__(self):
        """
        Length of a transcript is the sum of its exon lengths
        """
        # cache the length once you compute it
        if not hasattr(self, "_length"):
            self._length = sum(len(exon) for exon in self.exons)
        return self._length

    @property
    def exons(self):
        if not hasattr(self, "_exons"):
            columns = ['exon_number', 'exon_id']
            results = self.db.query(
                columns,
                filter_column='transcript_id',
                filter_value=self.id,
                feature='exon')

             # fill this list in its correct order (by exon_number) by using
             # the exon_number as a 1-based list offset
            exons = [None] * len(results)

            for entry in results:
                exon_number, exon_id = entry

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
                    self.transcript_name
                )
            self._exons = exons

        return self._exons


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

    def _transcript_feature_positions(self, feature, require_full_codon=False):
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
            # Python ranges are [inclusive, excsluve) we have to increment
            # the end position
            range_positions = range(start, end+1)
            results.extend(range_positions)
        if require_full_codon:
            if len(results) != 3:
                raise ValueError(
                    "Expected 3 positions for %s of %s but got %d" % (
                        feature,
                        self.id,
                        len(results)))
        return results

    @property
    def contains_start_codon(self):
        """
        Does this transcript have an annotated start_codon entry in Ensembl?
        """
        start_codons = _transcript_feature_position_ranges(
            'start_codon', required=False)
        return len(start_codons) > 0

    @property
    def contains_stop_codon(self):
        """
        Does this transcript have an annotated stop_codon entry in Ensembl?
        """
        stop_codons = _transcript_feature_position_ranges(
            'stop_codon', required=False)
        return len(stop_codons) > 0

    @property
    def start_codon_positions(self):
        """
        Chromosomal positions of nucleotides in start codon.
        """
        return self._transcript_feature_positions(
            'start_codon', require_full_codon=True)

    @property
    def stop_codon_positions(self):
        """
        Chromosomal positions of nucleotides in stop codon.
        """
        return self._transcript_feature_positions(
            'stop_codon', require_full_codon=True)

    @property
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

    def spliced_offset(self, position):
        """
        Convert from an absolute chromosomal position to the offset into
        this transcript's spliced mRNA.

        Position must be inside some exon (otherwise raise exception).
        """
        if not isinstance(position, (int, long)):
            raise TypeError(
                "Expected position to be int, got %s : %s" % (
                    position, type(position)))
        elif position < self.start or position > self.end:
            raise ValueError(
                "Invalid position: %d (must be between %d and %d)" % (
                    position,
                    self.start,
                    self.end))

        # offset from beginning of unspliced transcript (including introns)
        unspliced_offset = self.position_offset(position)
        total_spliced_offset = 0
        # traverse exons in order of their appearance on the strand
        # Since bsolute positions may decrease if on the negative strand,
        # we instead use unspliced offsets.
        for exon in self.exons:
            exon_start_unspliced_offset, exon_end_unspliced_offset = \
                self.offset_range(exon.start, exon.end)
            if unspliced_offset <= exon_end_unspliced_offset:
                exon_offset = unspliced_offset - exon_start_unspliced_offset
                return total_spliced_offset + exon_offset
            else:
                total_spliced_offset += len(exon)
        raise ValueError(
            "Couldn't find position %d on any exon of %s" % (
                position, self.id))

    @property
    def start_codon_spliced_offsets(self):
        """
        Offsets from start of spliced mRNA transcript
        of nucleotides in start codon.
        """
        return [
            self.spliced_offset(position)
            for position
            in self.start_codon_positions
        ]

    @property
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

    @property
    def stop_codon_spliced_offsets(self):
        """
        Offsets from start of spliced mRNA transcript
        of nucleotides in stop codon.
        """
        return [
            self.spliced_offset(position)
            for position
            in self.stop_codon_positions
        ]

    @property
    def coding_sequence_position_ranges(self):
        """
        Return absolute chromosome position ranges for CDS fragments
        of this transcript
        """
        return self._transcript_feature_position_ranges('CDS')

    @property
    def coding_sequence_length(self):
        total = 0
        for (start, stop) in self.coding_sequence_position_ranges:
            assert start <= stop, \
                "Invalid offset range [%d..%d] for transcript %s" % (
                    start, stop, self.id)
            total += stop - start + 1
        return total

    @property
    def complete(self):
        """
        Consider a transcript complete if it has start and stop codons and
        a coding sequence whose length is divisible by 3
        """
        return (
            self.contains_start_codon and
            self.contains_stop_codon and
            self.coding_sequence_length % 3 == 0
        )

    @property
    def sequence(self):
        """
        Spliced cDNA sequence of transcript
        (includes 5' UTR, coding sequence, and 3' UTR)
        """
        if self._sequence is None:
            if self.id not in self.reference:
                raise ValueError(
                    "No sequence for transcript %s in reference %s" % (
                        self.id, self.reference))
            # fetch this transcript's Sequence from the attached
            # ReferenceTranscripts object
            self._sequence = self.reference.transcript_sequence(self.id)
        return self._sequence

    @property
    def first_start_codon_spliced_offset(self):
        """
        Offset of first nucleotide in start codon into the spliced mRNA
        (excluding introns)
        """
        start_offsets = self.start_codon_spliced_offsets
        return min(start_offsets)


    @property
    def last_stop_codon_spliced_offset(self):
        """
        Offset of last nucleotide in stop codon into the spliced mRNA
        (excluding introns)
        """
        stop_offsets = self.stop_codon_spliced_offsets
        return max(stop_offsets)


    @property
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
        return self.sequence[start:end+1]

    @property
    def five_prime_utr_sequence(self):
        """
        cDNA sequence of 5' UTR
        (untranslated region at the beginning of the transcript)
        """
        return self.sequence[:self.first_start_codon_spliced_offset]

    @property
    def three_prime_utr_sequence(self):
        """
        cDNA sequence of 3' UTR
        (untranslated region at the end of the transcript)
        """
        return self.sequence[self.last_stop_codon_spliced_offset+1:]


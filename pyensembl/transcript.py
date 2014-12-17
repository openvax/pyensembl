from locus import Locus
from exon import Exon

from pyfaidx import Sequence

class Transcript(Locus):
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

        # Lazily fetch sequence for this Transcript.
        # Doing this in case we're constructing many Transcripts
        # and not using the sequence, avoid the memory/performance overhead
        # of fetching and storing sequences from a FASTA file.
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

    def _transcript_feature_position_range(self, feature):
        """
        Get unique start and end positions for feature,
        raise an error if feature is absent or has multiple entries
        for this transcript.
        """
        ranges = self._transcript_feature_position_ranges(
            feature, required=True)
        if len(ranges) > 1:
            raise ValueError(
                "Expected %s to be unique for %s but got %d entries" % (
                    feature, self.id, len(ranges)))
        return ranges[0]

    @property
    def contains_start_codon(self):
        start_codons = _transcript_feature_position_ranges('start_codon')
        return len(start_codons) > 0

    @property
    def contains_stop_codon(self):
        stop_codons = _transcript_feature_position_ranges('stop_codon')
        return len(stop_codons) > 0

    @property
    def start_codon_position_range(self):
        """
        Smallest and largest chromosomal positions of nucleotides in start codon
        """
        return self._transcript_feature_position_range('start_codon')

    @property
    def stop_codon_position_range(self):
        """
        Smallest and largest chromosomal positions of nucleotides in stop codon
        """
        return self._transcript_feature_position_range('stop_codon')

    @property
    def start_codon_offset_range(self):
        start, end = self.start_codon_position_range
        return self.offset_range(start, end)

    @property
    def stop_codon_offset_range(self):
        start, end = self.stop_codon_position_range
        return self.offset_range(start, end)

    @property
    def coding_sequence_position_ranges(self):
        """
        Return absolute chromosome position ranges for CDS fragments
        of this transcript
        """
        return self._transcript_feature_position_ranges('CDS')

    @property
    def coding_sequence_offset_ranges(self):
        """
        Return offsets from start of this transcript for CDS fragments
        """
        ranges = self._transcript_feature_position_ranges('CDS')
        return [self.offset_range(r) for r in ranges]

    @property
    def coding_sequence_length(self):
        total = 0
        for (start, stop) in self.coding_sequence_offset_ranges:
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
    def coding_sequence(self):
        start_offset, _ = self.start_codon_offset_range
        _, stop_offset = self.stop_codon_offset_range
        return self.sequence[start_offset:stop_offset]

    @property
    def five_prime_utr_sequence(self):
        start_offset, _ = self.start_codon_offset_range
        return self.sequence[:start_offset]

    @property
    def three_prime_utr_sequence(self):
        _, end_offset = self.stop_codon_offset_range
        return self.sequence[end_offset:]


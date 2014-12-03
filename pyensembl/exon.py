from locus import Locus

import numpy as np

class Exon(Locus):
    def __init__(self, exon_id, db):

        if not isinstance(exon_id, (unicode, str)):
            raise TypeError(
                "Expected exon ID to be string, got %s : %s" % (
                transcript_id, type(transcript_id)))

        self.id = exon_id
        self.db = db

        columns = [
            'seqname',
            'start',
            'end',
            'strand',
            'frame',
            'exon_number',
            'gene_name',
            'gene_id',
            'transcript_name',
            'transcript_id'
        ]
        columns_str = ", ".join(columns)

        query = """
            SELECT %s
            FROM ensembl
            WHERE exon_id = ?
            AND feature='exon'
        """ % columns_str

        cursor = db.execute(query, [exon_id])
        result = cursor.fetchone()

        if result is None:
            raise ValueError("Exon ID not found: %s" % exon_id)

        result_dict = {}
        for i, column_name in columns:
            result_dict[column_name] = result[i]

        seqname = result['seqname']
        start = int(result['start'])
        end = int(result['end'])
        strand = result['strand']

        Locus.__init__(self, seqname, start, end, strand)

        self.frame = int(result['frame'])

        assert 0 <= self.frame <= 2, "Invalid exon frame: %d" % self.frame

        self.exon_number = int(result['exon_number'])
        assert self.exon_number >= 0

        self.transcript_id = result['transcript_id']
        self.transcript_name = result['transcript_name']
        self.gene_name = result['gene_name']
        self.gene_id = result['gene_id']


    def __str__(self):
        return "Exon(id=%s, transcript_id=%s, gene=%s, frame=%s, number=%d)" % (
                    self.id,
                    self.transcript_id,
                    self.gene_name,
                    self.frame,
                    self.exon_number)

    def __repr__(self):
        return str(self)

    def _query_exon_feature_position(feature):
        query = """
            SELECT seqname, start, strand
            WHERE feature=%s
            AND transcript_id = ?
            AND exon_number = ?
        """ % (feature,)
        cursor = db.execute(query, [self.transcript_id, self.exon_number])
        return cursor.fetchall()


    @property
    def contains_start_codon(self):
        """
        Does this exon contain a start codon?
        """
        if not hasattr(self, "_contains_start_codon"):
            results = self._query_exon_feature_position('start_codon')
            self._contains_start_codon = len(results) > 0
        return self._contains_start_codon

    @property
    def contains_stop_codon(self):
        """
        Does this exon contain a stop codon ?
        """
        if not hasattr(self, "_contains_stop_codon"):
            results = self._query_exon_feature_position('stop_codon')
            self._contains_stop_codon = len(results) > 0
        return self._contains_stop_codon


    def _local_feature_position(feature):
        """
        Return the position of feature in coordinates relative
        to the start of this exon.
        """
        if feature not in {'start_codon', 'stop_codon', 'UTR'}:
            raise ValueError("Invalid exon feature: %s" % feature)

        results = self._query_exon_feature_position(feature)

        if len(results) == 0:
            raise ValueError(
                "Exon %d of %s does not contain %s" % (
                    self.exon_number, self.transcript_name, feature))

        # in case there are multiple results, choose the
        # first, which is either a higher or lower position depending
        # on the strand
        earliest_position = np.inf if self.strand == "+" else 0
        for entry in start_codon_results:
            seqname, position, strand = entry
            assert seqname == self.contig, \
                "Wrong contig for exon: %s (should be %s)" % (
                    seqname, self.contig)
            assert strand in ("+", "-"), "Invalid strand: %s" % strand
            assert strand == self.strand

            assert isinstance(position, (int,long)), \
                "Invalid type %s for position %s" % (type(position), position)

            if position < earliest_position:
                earliest_position = position

        local_position = self.position_offset(earliest_position)
        if local_position < 0:
            raise ValueError(
                "%s starts before exon %d of %s" % (
                    feature, self.exon_number, self.transcript_name))
        return local_position

    @property
    def start_codon_offset(self):
        """
        How many bases from the beginning of the exon (starting from 0)
        is the first base of the start codon?
        """
        return self._local_feature_position('start_codon')


    @property
    def stop_codon_offset(self):
        """
        How many bases from the beginning of the exon (starting from 0)
        is the first base of the stop codon?
        """
        return self._local_feature_position("stop_codon")


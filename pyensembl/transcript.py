from locus import Locus
from exon import Exon

class Transcript(Locus):
    def __init__(self, transcript_id, db):

        if not isinstance(transcript_id, (unicode, str)):
            raise TypeError(
                "Expected transcript ID to be string, got %s : %s" % (
                transcript_id, type(transcript_id)))

        self.id = transcript_id
        self.db = db
        query = """
            SELECT
                transcript_name,
                seqname, start, end, strand,
                gene_name, gene_id
            FROM ensembl
            WHERE transcript_id = ?
            AND feature='transcript'
        """
        cursor = db.execute(query, [transcript_id])

        result = cursor.fetchone()
        if result is None:
            raise ValueError("Transcript ID not found: %s" % transcript_id)

        transcript_name, contig, start, end, strand, gene_name, gene_id = result

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



    def __str__(self):
        return "Transcript(id=%s, name=%s, gene_name=%s)" % (
                    self.id, self.name, self.gene_name)

    def __repr__(self):
        return str(self)

    @property
    def exons(self):
        if not hasattr(self, "_exons"):
            exon_ids_query = """
                SELECT exon_number, exon_id
                FROM ensembl
                WHERE transcript_id = ?
                AND feature='exon'
            """
            cursor = self.db.execute(exon_ids_query, [self.id])
            results = cursor.fetchall()

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

    def _transcript_feature_positions(self, feature, required=True):
        """
        Find features (such as start codon) of this transcript.
        """

        if feature not in self._TRANSCRIPT_FEATURES:
            raise ValueError("Invalid transcript feature: %s" % feature)

        query = """
            SELECT DISTINCT start, end
            FROM ensembl
            WHERE feature = ?
            AND transcript_id = ?
        """
        query_params = [
            feature,
            self.id
        ]
        cursor = self.db.execute(query, query_params)
        results = cursor.fetchall()
        if required and len(results) == 0:
            raise ValueError(
                "Transcript %s does not contain feature %s" % (
                    self.id, feature))
        return results

    def _transcript_feature_positions(self, feature):
        """
        Get unique start and end positions for feature,
        raise an error if feature is absent or has multiple entries
        for this transcript.
        """
        results = self._transcript_feature_positions(feature, required=True)
        if len(results) > 1:
            raise ValueError(
                "Expected %s to be unique for %s but got %d entries" % (
                    feature, self.id, len(results)))
        return results[0]


    @property
    def start_codon_locus(self):
        return Locus(self.contig, start, end, self.strand)

    @property
    def start_codon_offset(self):

        results = self._query_transcript_feature_positions(
            feature, required=True)


    def stop_codon_offset(self):
        pass

    @property
    def coding_sequence_length(self):
        pass

    @property
    def contains_start_codon(self):
        pass

    @property
    def contains_stop_codon(self):
        pass

    @property
    def complete(self):
        """
        Consider a transcript complete if it has start and stop codons and
        a coding sequence whose length is divisible by 3
        """
        return self.contains_start_codon and self.contains_stop_codon



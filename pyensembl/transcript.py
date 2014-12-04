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


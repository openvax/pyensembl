from locus import Locus
from transcript import Transcript

class Gene(Locus):

    def __init__(self, gene_id, db):
        if not isinstance(gene_id, (unicode, str)):
            raise TypeError(
                "Expected gene ID to be string, got %s : %s" % (
                    gene_id, type(gene_id))
            )

        self.id = gene_id
        self.db = db
        query = """
            SELECT gene_name, seqname, start, end, strand, gene_biotype
            FROM ensembl
            WHERE gene_id = ?
            AND feature = 'gene'
        """
        cursor = db.execute(query, [gene_id])
        result = cursor.fetchone()
        if not result:
            raise ValueError("Gene ID not found: %s" % gene_id)

        gene_name, contig, start, end, strand, biotype = result
        if not gene_name:
            raise ValueError("Missing name for gene with ID = %s" % gene_id)
        self.name = gene_name

        Locus.__init__(self, contig, start, end, strand)

        if not biotype:
            raise ValueError(
                "Missing gene_biotype for gene with ID = %s" % gene_id
            )
        self.biotype = biotype

    def __str__(self):
        return "Gene(id=%s, name=%s)" % (self.id, self.name)

    def __repr__(self):
        return str(self)

    @property
    def transcripts(self):
        """
        Property which dynamically construct transcript objects for all
        transcript IDs associated with this gene.
        """
        if not hasattr(self, "_transcripts"):
            transcript_ids_query = """
                SELECT transcript_id
                FROM ensembl
                WHERE gene_id = ?
                AND feature = 'transcript'
            """
            cursor = self.db.execute(transcript_ids_query, [self.id])
            results = cursor.fetchall()

            # We're doing a SQL query for each transcript ID to fetch
            # its particular information, might be more efficient if we
            # just get all the columns here, but how do we keep that modular?
            transcripts = [
                Transcript(result[0], self.db)
                for result in results
            ]
            self._transcripts = transcripts
        return self._transcripts

    @property
    def exons(self):
        if not hasattr(self, "_exons"):
            exons_dict = {}
            for transcript in self.transcripts:
                for exon in transcript.exons:
                    if exon.id not in exons_dict:
                        exons_dict[exon.id] = exon
            self._exons = list(exons_dict.values())
        return self._exons


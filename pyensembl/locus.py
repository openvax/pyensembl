
def normalize_chromosome(c):

    if isinstance(c, (int, long)):
        if c == 0:
            raise ValueError("Contig cannot be 0")
        c = str(c)

    if not isinstance(c, (str, unicode)):
        raise TypeError(
            "Expected contig name to be str, got %s : %s" % (c, type(c))
        )

    if len(c) == 0:
        raise ValueError("Contig name cannot be empty string")

    # only strip off lowercase chr since some of the non-chromosomal contigs
    # start with "CHR"
    if c.startswith("chr"):
        c = c[3:]

    # standardize mitochondrial genome to be "MT"
    if c == "M":
        return "MT"
    # just in case someone is being lazy, capitalize "X" and "Y"
    elif c == "x":
        return "X"
    elif c == "y":
        return "Y"
    else:
        return c

def normalize_strand(strand):
    if strand == 1:
        return "+"
    elif strand == -1:
        return "-"

    if not isinstance(strand, (str, unicode)):
        raise TypeError("Expected strand to be string, got %s : %s" % (
            strand, type(strand)))
    elif len(strand) == 0:
        raise ValueError("Strand cannot be empty string")
    elif len(strand) > 1:
        raise ValueError("Invalid strand: %s" % (strand,))

    return strand

class Locus(object):
    """
    Base class for any entity which can be localized at a range of positions
    on a particular chromosome/contig.
    """

    def __init__(
            self,
            contig,
            start,
            end,
            strand,
            gene_id=None,
            gene_name=None,
            transcript_id=None,
            transcript_name=None,
            exon_id=None):
        """
        contig : str
            Chromosome or other sequence name in the reference assembly

        start : int
            Start position of locus on the contig

        end : int
            Inclusive end position on the contig

        strand : str
            Should we read the locus forwards ('+') or backwards ('-')?
        """
        self.contig = normalize_chromosome(contig)

        self.start = int(start)
        self.end = int(end)

        if self.start == 0:
            raise ValueError("Expected start > 0 (using base 1 coordinates)")
        elif self.end == 0:
            raise ValueError("Expected end > 0 (using base 1 coordinates)")
        elif self.end < self.start:
            raise ValueError(
                "Expected start <= end, got start = %d, end = %d" % (
                    self.start, self.end)
            )

        self.strand = normalize_strand(strand)


    @property
    def on_forward_strand(self):
        return self.strand == "+"

    @property
    def on_positive_strand(self):
        """
        Alias for forward strand
        """
        return self.on_forward_strand

    @property
    def on_backward_strand(self):
        return self.strand == "-"

    @property
    def on_negative_strand(self):
        """
        Alias for backward strand
        """
        return self.on_backward_strand

    @property
    def length(self):
        if self.on_forward_strand:
            return self.end - self.start + 1
        else:
            return self.start - self.end + 1

    def position_offset(self, position):
        if self.on_forward_strand:
            return position - self.start
        else:
            return self.start - position

    def on_contig(self, contig):
        return normalize_chromosome(contig) == self.contig

    def on_strand(self, strand):
        return normalize_strand(strand) == self.strand

    def can_overlap(self, contig, strand=None):
        """
        Is this locus on the same contig and (optionally) on the same strand?
        """
        return (
            self.on_contig(contig)
            and
            (strand is None or self.on_strand(strand)))

    def overlaps(self, contig, start, end, strand=None):
        """
        Does this locus overlap with a given range of positions?

        Since locus position ranges are inclusive, we should make sure
        that e.g. chr1:10-10 overlaps with chr1:10-10
        """
        return (
            self.can_overlap(contig, strand)
            and
            end >= self.start
            and
            start <= self.end)

    def overlaps_locus(self, other_locus):
        return self.overlaps(
            other_locus.contig,
            other_locus.start,
            other_locus.end,
            other_locus.strand)

    def contains(self, contig, start, end, strand=None):
        return (
            self.can_overlap(contig, strand)
            and
            start >= self.start
            and
            end <= self.end)

    def contains_locus(self, other_locus):
        return self.contains(
            other_locus.contig,
            other_locus.start,
            other_locus.end,
            other_locus.strand)

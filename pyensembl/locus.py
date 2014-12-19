
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
            strand):
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
        self.strand = normalize_strand(strand)

        start = int(start)
        end = int(end)

        if start == 0:
            raise ValueError("Expected start > 0 (using base 1 coordinates)")
        elif end == 0:
            raise ValueError("Expected end > 0 (using base 1 coordinates)")

        if end < start:
            raise ValueError(
                "Expected start <= end, got start = %d, end = %d" % (
                    start, end))
        self.start = start
        self.end = end

    def __str__(self):
        return "Locus(contig=%s, start=%s, end=%s, strand=%s)" % (
            self.contig, self.start, self.end, self.strand)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.end - self.start + 1

    def __eq__(self, other):
        return (
            self.contig == other.contig and
            self.start == other.start and
            self.end == other.end and
            self.strand == other.strand
        )

    @property
    def length(self):
        return self.end - self.start + 1

    def position_offset(self, position):
        if position > self.end or position < self.start:
            raise ValueError(
                "Position %d outside valid range %d..%d of %s" % (
                    position, self.start, self.end, self))
        elif self.on_forward_strand:
            return position - self.start
        else:
            return self.end - position

    def offset_range(self, start, end):
        """
        Database start/end entries are always ordered such that
        start < end. This makes computing a relative position (e.g. of a stop
        codon relative to its transcript) complicated since the "end"
        position of a backwards locus is actually earlir on the strand.
        This function correctly selects a start vs. end value depending
        on this locuses's strand and determines that position's offset from
        the earliest position in this locus.
        """
        assert start <= end, \
            "Locations should always have start < end, got start=%d, end=%d" % (
                start, end)

        if start < self.start or end > self.end:
            raise ValueError("Range (%d, %d) falls outside %s" % (
                start, end, self))

        if self.on_forward_strand:
            return (start - self.start, end - self.start)

        else:
            return (self.end - end, self.end - start)

    def on_contig(self, contig):
        return normalize_chromosome(contig) == self.contig

    def on_strand(self, strand):
        return normalize_strand(strand) == self.strand

    @property
    def on_forward_strand(self):
        return self.on_strand("+")

    @property
    def on_positive_strand(self):
        return self.on_forward_strand

    @property
    def on_backward_strand(self):
        return self.on_strand("-")

    @property
    def on_negative_strand(self):
        return self.on_backward_strand

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

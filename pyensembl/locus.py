def normalize_chromosome(c):
    if isinstance(c, (int, long)):
        assert c > 0, "Invalid chromosome %d" % c
        c = str(c)

    assert isinstance(c, str), "Invalid chromosome %s : %s" % (c, type(c))

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

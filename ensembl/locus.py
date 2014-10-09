def normalize_chromosome(c):
    if isinstance(c, (int, long)):
        assert c > 0
        return c

    assert isinstance(c, str), "Invalid chromosome %s : %s" % (c, type(c))
    if c.startswith("chr"):
        c = c[3:]

    # if after trimming 'chr' prefix we get a number, return it as an int
    try:
        return int(c)
    except:
        pass

    # in case human-entered chromosomes come in as e.g. 'y', convert to 'Y'
    c = c.upper()

    # standardize mitochondrial genome to be "MT"
    if c == "M":
        c = "MT"

    return c

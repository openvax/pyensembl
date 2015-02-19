"""
Helper functions for searching over collections of PyEnsembl objects
"""

def find_nearest_locus(start, end, loci):
    """
    Finds nearest locus (object with method `distance_to_interval`) to the
    interval defined by the given `start` and `end` positions.
    Returns the distance to that locus, along with the locus object itself.
    """
    best_distance = float("inf")
    best_locus = None
    for locus in loci:
        distance = locus.distance_to_interval(start, end)

        if best_distance > distance:
            best_distance = distance
            best_locus = locus

    return best_distance, best_locus

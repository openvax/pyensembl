"""Compare complete fresh SQLite indexing with pandas and Polars parsing.

Run from the checkout: python scripts/benchmark_gtf_indexing.py /path/to/file.gtf.gz
Each result includes parsing, feature synthesis, SQL loading and indexing.
"""

import argparse
import json
import logging
from statistics import median
from tempfile import TemporaryDirectory
from time import perf_counter
from unittest.mock import patch

import gtfparse
from pyensembl.database import Database


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf")
    parser.add_argument("--runs", type=int, default=2)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be positive")
    logging.getLogger().setLevel(logging.WARNING)
    reader = gtfparse.read_gtf
    timings = {"polars": [], "pandas": []}
    for run in range(args.runs):
        order = ("polars", "pandas") if run % 2 == 0 else ("pandas", "polars")
        for result_type in order:
            def read(*a, **kw):
                return reader(*a, **{**kw, "result_type": result_type})

            with TemporaryDirectory(prefix="pyensembl-index-benchmark-") as directory:
                with patch.object(gtfparse, "read_gtf", read), Database(
                    args.gtf, cache_directory_path=directory
                ) as database:
                    start = perf_counter()
                    database.create()
                    seconds = perf_counter() - start
                    timings[result_type].append(seconds)
                    print(json.dumps({"run": run + 1, "result_type": result_type,
                                      "seconds": seconds}), flush=True)
    print(json.dumps({"median_seconds": {name: median(values)
                                         for name, values in timings.items()}}))


if __name__ == "__main__":
    main()

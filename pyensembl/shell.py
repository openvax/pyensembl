#!/usr/bin/env python

# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Manipulate pyensembl's local cache.

    %(prog)s {install,download,index} [--release XXX ...]

To install the latest Ensembl release:
    %(prog)s install

To install particular Ensembl release(s):
    %(prog)s install --release 75 77

To install any Genomic database:
    %(prog)s install --gtf-url URL --transcript-fasta-url URL --protein-fasta-url URL

"""

from __future__ import absolute_import
import argparse
from .ensembl_release import EnsemblRelease
from .genome import Genome
from .genome_source import GenomeSource
from .release_info import MAX_ENSEMBL_RELEASE

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--gtf-url",
                        type=str,
                        nargs=1,
                        help="URL for a GTF file containing annotations.")
    parser.add_argument("--transcript-fasta-url",
                        type=str,
                        nargs=1,
                        default=None,
                        help="URL for a FASTA file containing transcript data.")
    parser.add_argument("--protein-fasta-url",
                        type=str,
                        nargs=1,
                        default=None,
                        help="URL for a FASTA file containing protein data.")
    parser.add_argument("--release",
        type=int,
        nargs="+",
        default=[MAX_ENSEMBL_RELEASE],
        help="Ensembl release. Defaults to latest release %(default)s. "
             "Multiple releases may be specified.")
    parser.add_argument("action", choices=("install", "download", "index"),
        help=""install" will download and index any data that is  not "
        "currently downloaded or indexed. \"download\" will download data, "
        "regardless of whether it is already downloaded. \"index\" will index "
        "data, regardless of whether it is already indexed, and will raise "
        "an error if the data is not already downloaded.")

    args = parser.parse_args()
    genomes = []

    # If specific genome source URLs are provided, use those
    if "gtf_url" in args:
        genome_source = GenomeSource(
            gtf_url=args.gtf_url,
            transcript_fasta_url=args.transcript_fasta_url,
            protein_fasta_url=args.protein_fasta_url)
        genomes.append(Genome(genome_source=genome_source))
    # Otherwise, use Ensembl release information
    else:
        for release in args.release:
            genomes.append(EnsemblRelease(release))

    for genome in genomes:
        assert args.action in ["install", "download", "index"], \
            "Invalid action: %s" % args.action
        if args.action == "install":
            genome.install()
        elif args.action == "download":
            genome.download()
        else:
            genome.index()


if __name__ == "__main__":
    run()

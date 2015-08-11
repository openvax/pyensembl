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

To install any genome:
    %(prog)s install --reference_name "GRCh38" --gtf-path_or_url URL_OR_PATH --transcript-fasta-path_or_url URL_OR_PATH --protein-fasta-path-or-url URL_OR_PATH

"""

from __future__ import absolute_import
import argparse
from .ensembl_release import EnsemblRelease
from .genome import Genome
from .release_info import MAX_ENSEMBL_RELEASE

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        "--force",
        default=False,
        action="store_true",
        help="Force download and indexing even if files already exist locally")

    root_group = parser.add_mutually_exclusive_group()

    release_group = root_group.add_argument_group()
    release_group.add_argument(
        "--release",
        type=int,
        nargs="+",
        default=None,
        help=("Ensembl release. Defaults to latest release %(default)s. "
              "Multiple releases may be specified."))

    release_group.add_argument(
        "--species",
        default="human",
        help="Which species to download Ensembl data for (default=%(default)s.")

    path_group = root_group.add_argument_group()
    path_group.add_argument(
        "--reference-name",
        type=str,
        default=None,
        help="Name of the reference, e.g. GRCh38")
    path_group.add_argument(
        "--annotation-name",
        default=None,
        help="Name of annotation source (e.g. refseq)")
    path_group.add_argument(
        "--annotation-version",
        default=None,
        help="Version of annotation database")
    path_group.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="URL or local path to a GTF file containing annotations.")
    path_group.add_argument(
        "--transcript-fasta",
        type=str,
        default=None,
        help="URL or local path to a FASTA file containing transcript data.")
    path_group.add_argument(
        "--protein-fasta",
        type=str,
        default=None,
        help="URL or local path to a FASTA file containing protein data.")

    parser.add_argument("action", choices=("install", "download", "index"),
        help="\"install\" will download and index any data that is  not "
        "currently downloaded or indexed. \"download\" will download data, "
        "regardless of whether it is already downloaded. \"index\" will index "
        "data, regardless of whether it is already indexed, and will raise "
        "an error if the data is not already downloaded.")

    args = parser.parse_args()

    genomes = []
    # If specific genome source URLs are provided, use those
    if args.gtf or args.transcript_fasta or args.protein_fasta:
        if args.release:
            raise ValueError(
                "An Ensembl release cannot be specified if "
                "specific paths are also given")
        if not args.reference_name:
            raise ValueError("Must specify a reference name")
        if not args.annotation_name:
            raise ValueError("Must specify the name of the annotation source")
        genomes.append(Genome(
            reference_name=args.reference_name,
            annotation_name=args.annotation_name,
            gtf_path_or_url=args.gtf,
            transcript_fasta_path_or_url=args.transcript_fasta,
            protein_fasta_path_or_url=args.protein_fasta,
            auto_download=True,
            force_download=args.force))
    # Otherwise, use Ensembl release information
    else:
        versions = args.release if args.release else [MAX_ENSEMBL_RELEASE]
        for version in versions:
                genomes.append(
                    EnsemblRelease(
                        version,
                        species=args.species,
                        auto_download=True))

    for genome in genomes:
        if args.action not in ["download"]:
            "Invalid action: %s" % args.action
        if args.action == "install":
            genome.install()
        elif args.action == "download":
            genome.download()
        else:
            genome.index()


if __name__ == "__main__":
    run()

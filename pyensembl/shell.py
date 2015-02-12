#!/usr/bin/env python
"""
A shell wrapper around various PyEnsembl commands.
"""
import argparse

import ensembl_release
from release_info import MAX_ENSEMBL_RELEASE


def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    subparsers = parser.add_subparsers(dest='program')
    download_parser = subparsers.add_parser(
        'update',
        usage=('retrieve all Ensembl data for a release, some of '
               'which is installed into a local SQLite database'))
    download_parser.add_argument(
        'release',
        nargs='?',
        default=MAX_ENSEMBL_RELEASE,
        help=('specify the Ensembl release that you would like to '
              'download and install (defaults to latest release)'))
    download_parser.add_argument(
        '-a', '--annotations-only',
        action='store_true',
        dest='annotations_only',
        help=('retrieve only Ensembl annotation data, and install it '
              'into a local SQLite database'))
    download_parser.add_argument(
        '-t', '--transcripts-only',
        action='store_true',
        dest='transcripts_only',
        help='retrieve only Ensembl transcript data')

    args = parser.parse_args()
    if args.program == 'update':
        data = ensembl_release.EnsemblRelease(args.release)
        if args.annotations_only:
            data.download_annotations()
        elif args.transcripts_only:
            data.download_transcripts()
        else:
            data.download_all()


if __name__ == '__main__':
    run()

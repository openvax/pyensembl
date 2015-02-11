#!/usr/bin/env python
"""
This tool is a shell wrapper around various PyEnsembl commands.
"""
import argparse

import ensembl_release
from release_info import MAX_ENSEMBL_RELEASE


def download_annotations(release=MAX_ENSEMBL_RELEASE):
    data = ensembl_release.EnsemblRelease(release)
    data.download_annotations()


def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    subparsers = parser.add_subparsers(dest='program')
    download_parser = subparsers.add_parser(
        'download', help ='retrieve and install Ensembl data')
    download_parser.add_argument(
        'release',
        help=('Specify the Ensembl release that you would like to'
              'download and install (defaults to latest release)'))
    args = parser.parse_args()
    if args.program == 'download':
        download_annotations(args.release)


if __name__ == '__main__':
    run()

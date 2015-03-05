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
A shell wrapper around various PyEnsembl commands.

Example:

    pyensembl 75 77 install
"""
import argparse

import ensembl_release
from release_info import MAX_ENSEMBL_RELEASE


def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        'release',
        nargs='*',
        default=MAX_ENSEMBL_RELEASE,
        help=('specify the Ensembl release(s) that you would like to '
              'download and install (defaults to latest release)'))
    subparsers = parser.add_subparsers(dest='program')
    subparsers.add_parser(
        'install',
        help=('download and index any data for this release that '
               'is not yet downloaded and/or indexed'))
    subparsers.add_parser(
        'download',
        help=('download all data for this release, regardless of '
              'whether it is already downloaded'))
    subparsers.add_parser(
        'index',
        help=('index all data for this release, regardless of '
              'whether it is already indexed, and raises an error if '
              'data is not yet downloaded'))

    args = parser.parse_args()
    releases = args.release
    if type(releases) == int:
        releases = [releases]
    for release in releases:
        data = ensembl_release.EnsemblRelease(release)
        if args.program == 'install':
            data.install()
        elif args.program == 'download':
            data.download()
        elif args.program == 'index':
            data.index()


if __name__ == '__main__':
    run()

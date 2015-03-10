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

To install the latest ensembl release:
    %(prog)s install

To install particular release(s):
    %(prog)s install --release 75 77

"""
from __future__ import absolute_import
import argparse
from . import ensembl_release
from .release_info import MAX_ENSEMBL_RELEASE

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--release',
        type=int,
        nargs="+",
        default=[MAX_ENSEMBL_RELEASE],
        help='Ensembl release. Defaults to latest release %(default)s. '
             'Multiple releases may be specified.')
    parser.add_argument('action', choices=('install', 'download', 'index'),
        help="'install' will download and index any data that is  not "
        "currently downloaded or indexed. 'download' will download data, "
        "regardless of whether it is already downloaded. 'index' will index "
        "data, regardless of whether it is already indexed, and will raise "
        "an error if the data is not already downloaded.")

    args = parser.parse_args()
    for release in args.release:
        ensembl = ensembl_release.EnsemblRelease(release)
        if args.action == 'install':
            ensembl.install()
        elif args.action == 'download':
            ensembl.download()
        elif args.action == 'index':
            ensembl.index()
        else:
            assert False, "shouldn't get here"

if __name__ == '__main__':
    run()

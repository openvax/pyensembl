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

from __future__ import print_function, division, absolute_import

MIN_ENSEMBL_RELEASE = 54
MAX_ENSEMBL_RELEASE = 85

def check_release_number(release):
    """
    Check to make sure a release is in the valid range of
    Ensembl releases.
    """
    try:
        release = int(release)
    except:
        raise ValueError("Invalid Ensembl release: %s" % release)

    if release < MIN_ENSEMBL_RELEASE or release > MAX_ENSEMBL_RELEASE:
        raise ValueError(
            "Invalid Ensembl releases %d, must be between %d and %d" % (
                release, MIN_ENSEMBL_RELEASE, MAX_ENSEMBL_RELEASE))
    return release

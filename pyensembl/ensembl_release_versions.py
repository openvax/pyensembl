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

from .config import (MAX_ENSEMBL_RELEASE, MAX_ENSEMBLGENOME_RELEASE,
                     MIN_ENSEMBL_RELEASE, MIN_ENSEMBLGENOME_RELEASE)


def check_release_number(release, database=None):
    """
    Check to make sure a release is in the valid range of Ensembl releases.
    """
    if release is None:
        return MAX_ENSEMBL_RELEASE if database is None else MAX_ENSEMBLGENOME_RELEASE
    try:
        release = int(release)
    except ValueError:
        raise ValueError("Invalid Ensembl release: %s" % release)
    if database is None:
        min_release = MIN_ENSEMBL_RELEASE
    else:
        min_release = MIN_ENSEMBLGENOME_RELEASE
    if release < min_release:
        raise ValueError(
            "Invalid Ensembl releases %d, must be greater than %d"
            % (release, min_release)
        )
    return release

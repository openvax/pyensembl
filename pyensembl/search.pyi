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

from typing import Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from .locus import Locus

def find_nearest_locus(start: int, end: int, loci: Iterable[Locus]):
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

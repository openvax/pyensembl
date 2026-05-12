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

from serializable import Serializable


class Protein(Serializable):
    """
    Lightweight view object exposing the protein identity for a transcript.
    Accessed via :attr:`Transcript.protein`.
    """

    def __init__(self, protein_id, protein_version=None):
        self.protein_id = protein_id
        self.protein_version = protein_version

    @property
    def id(self):
        """Alias for :attr:`protein_id`."""
        return self.protein_id

    @property
    def version(self):
        """Alias for :attr:`protein_version`."""
        return self.protein_version

    @property
    def versioned_protein_id(self):
        """``protein_id.protein_version`` when available, else ``protein_id``."""
        if self.protein_version is None:
            return self.protein_id
        return "%s.%d" % (self.protein_id, self.protein_version)

    @property
    def versioned_id(self):
        """Alias for :attr:`versioned_protein_id`."""
        return self.versioned_protein_id

    def __eq__(self, other):
        return (
            other.__class__ is Protein
            and self.protein_id == other.protein_id
            and self.protein_version == other.protein_version
        )

    def __hash__(self):
        return hash((self.protein_id, self.protein_version))

    def __str__(self):
        return "Protein(protein_id='%s', protein_version=%s)" % (
            self.protein_id,
            self.protein_version,
        )

    def __repr__(self):
        return str(self)

    def to_dict(self):
        return {
            "protein_id": self.protein_id,
            "protein_version": self.protein_version,
        }

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

from .config import SPECIES_DATA

# TODO: replace Serializable with data class


class Species(Serializable):
    """
    Container for combined information about a species name, its synonyn names
    and which reference to use for this species in each Ensembl release.
    """

    # as species instances get created, they get registered in these
    # dictionaries
    _latin_names_to_species = {}
    _common_names_to_species = {}
    _reference_names_to_species = {}

    @classmethod
    def register(
        cls, latin_name, synonyms, reference_assemblies, database=None
    ):
        """
        Create a Species object from the given arguments and enter into all the
        dicts used to look the species up by its fields.
        """
        species = Species(
            latin_name=latin_name,
            synonyms=synonyms,
            reference_assemblies=reference_assemblies,
            database=database,
        )
        cls._latin_names_to_species[species.latin_name] = species
        for synonym in synonyms:
            if synonym in cls._common_names_to_species:
                raise ValueError(
                    "Can't use synonym '%s' for both %s and %s"
                    % (synonym, species, cls._common_names_to_species[synonym])
                )
            cls._common_names_to_species[synonym] = species
        for reference_name in reference_assemblies:
            if reference_name in cls._reference_names_to_species:
                raise ValueError(
                    "Can't use reference '%s' for both %s and %s"
                    % (
                        reference_name,
                        species,
                        cls._reference_names_to_species[reference_name],
                    )
                )
            cls._reference_names_to_species[reference_name] = species
        return species

    @classmethod
    def all_registered_latin_names(cls):
        """
        Returns latin name of every registered species.
        """
        return list(cls._latin_names_to_species.keys())

    @classmethod
    def all_species_release_pairs(cls):
        """
        Generator which yields (species, release) pairs for all possible
        combinations.
        """
        for species_name in cls.all_registered_latin_names():
            species = cls._latin_names_to_species[species_name]
            for _, release_range in species.reference_assemblies.items():
                for release in range(release_range[0], release_range[1] + 1):
                    yield species_name, release

    def __init__(
        self, latin_name, synonyms=[], reference_assemblies={}, database=None
    ):
        """
        Parameters
        ----------
        latin_name : str

        synonyms : list of strings

        reference_assemblies : dict
            Mapping of names of reference genomes onto inclusive ranges of
            Ensembl releases Example: {"GRCh37": (54, 75)}
        """
        self.latin_name = latin_name.lower().replace(" ", "_")
        self.synonyms = synonyms
        self.reference_assemblies = reference_assemblies
        self.database = database
        self._release_to_genome = {}
        for genome_name, (start, end) in self.reference_assemblies.items():
            for i in range(start, end + 1):
                if i in self._release_to_genome:
                    raise ValueError(
                        "Ensembl release %d already has an associated genome"
                        % i
                    )
                self._release_to_genome[i] = genome_name

    def which_reference(self, ensembl_release):
        if ensembl_release not in self._release_to_genome:
            raise ValueError(
                "No genome for %s in Ensembl release %d"
                % (self.latin_name, ensembl_release)
            )
        return self._release_to_genome[ensembl_release]

    def __str__(self):
        return (
            "Species(latin_name='%s', synonyms=%s, reference_assemblies=%s, database=%s)"
            % (
                self.latin_name,
                self.synonyms,
                self.reference_assemblies,
                self.database,
            )
        )

    def __eq__(self, other):
        return (
            other.__class__ is Species
            and self.latin_name == other.latin_name
            and self.synonyms == other.synonyms
            and self.reference_assemblies == other.reference_assemblies
            and self.database == other.database
        )

    def to_dict(self):
        return {"latin_name": self.latin_name}

    @classmethod
    def from_dict(cls, state_dict):
        return cls._latin_names_to_species[state_dict["latin_name"]]

    def __hash__(self):
        return hash(
            (
                self.latin_name,
                tuple(self.synonyms),
                frozenset(self.reference_assemblies.items()),
                self.database,
            )
        )


def normalize_species_name(name):
    """
    If species name was "Homo sapiens" then replace spaces with underscores and
    return "homo_sapiens".

    Also replace common names like "human" with "homo_sapiens".
    """
    lower_name = name.lower().strip()

    # if given a common name such as "human", look up its latin equivalent
    if lower_name in Species._common_names_to_species:
        return Species._common_names_to_species[lower_name].latin_name

    return lower_name.replace(" ", "_")


def find_species_by_name(species_name):
    latin_name = normalize_species_name(species_name)
    if latin_name not in Species._latin_names_to_species:
        raise ValueError(
            "Species not found: %s, for non-Ensembl data see https://github.com/openvax/pyensembl#non-ensembl-data"
            % (species_name,)
        )
    return Species._latin_names_to_species[latin_name]


def check_species_object(species_name_or_object):
    """
    Helper for validating user supplied species names or objects.

    Return `Species` Object
    """
    if isinstance(species_name_or_object, Species):
        return species_name_or_object
    elif isinstance(species_name_or_object, str):
        return find_species_by_name(species_name_or_object)
    else:
        raise ValueError(
            "Unexpected type for species: %s : %s"
            % (species_name_or_object, type(species_name_or_object))
        )


for data in SPECIES_DATA:
    globals()[data["synonyms"][0]] = Species.register(
        latin_name=data["latin_name"],
        synonyms=data["synonyms"],
        reference_assemblies=data["reference_assemblies"],
        database=data.get("database", None),
    )

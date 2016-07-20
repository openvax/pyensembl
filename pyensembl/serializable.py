# Copyright (c) 2016. Mount Sinai School of Medicine
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

import json

class Serializable(object):
    """
    Base class for all PyEnsembl objects which provides default
    methods such as to_json, from_json, __reduce__, and from_dict
    and relies only on (1) a user-defined to_dict method and (2) that the
    keys of to_dict match the arguments to __init__.
    """

    def __str__(self):
        fields_str = ", ".join(
            "%s=%s" % (k, v) for (k, v) in sorted(self.to_dict().items()))
        return "%s(%s)" % (self.__class__.__name__, fields_str)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if self.__class__ is not other.__class__:
            return False
        self_dict = self.to_dict()
        other_dict = other.to_dict()
        if len(self_dict) != other_dict(other_dict):
            return False
        self_keys = set(self_dict.keys())
        other_keys = set(other_dict.keys())
        if self_keys != other_keys:
            return False
        for key in self_keys:
            if self_dict[key] != other_dict[key]:
                return False
        return True

    def to_dict(self):
        raise NotImplementedError("Method to_dict() not implemented for %s" % (
            self.__class__.__name__,))

    @classmethod
    def _tuple_to_class(cls, module_and_class_name):
        """
        Given the name of a module and a class it contains, imports that module
        and gets the class object from it.
        """
        module_string, class_name = module_and_class_name
        module = __import__(module_string)
        return getattr(module, class_name)

    @classmethod
    def _class_to_tuple(cls):
        """
        Given an class, return two strings:
            - fully qualified import path for its module
            - name of the class

        The class can be reconstructed from these two strings by calling
        _tuple_to_class.
        """
        return cls.__module__, cls.__name__

    def _object_to_primitive_types(self):
        """
        Given an instance of a Python object, returns a tuple whose
        first element is a primitive representation of the class and whose
        second element is a dictionary of instance data.
        """
        class_representation = self._class_to_tuple()
        state_dict = self.to_dict()
        return (class_representation, state_dict)

    @classmethod
    def _object_from_primitive_types(cls, class_and_data_pair):
        """
        Given a primitive representation of some object, reconstructs
        the class from its module and class names and then instantiates. Returns
        instance object.

        It's confusing that there are *two* class variables here:
            `cls` corresponds to Serializable
            `subclass` is the dynamically constructed subclass which we're
            trying to deserialize (and is presumably a subclass of Serializable)
        """
        class_representation, state_dict = class_and_data_pair
        subclass = cls._tuple_to_class(class_representation)
        return subclass.from_dict(state_dict)

    @classmethod
    def _initialize_nested_objects_in_state_dict(cls, state_dict):
        """
        Nested serializable objects will be represented as dictionaries so we
        allow manual reconstruction of those objects in this method.
        """
        return state_dict

    @classmethod
    def from_dict(cls, state_dict):
        """
        Given a dictionary of flattened fields (result of calling to_dict()),
        returns an instance.
        """
        state_dict = cls._initialize_nested_objects_in_state_dict(state_dict)
        return cls(**state_dict)

    def to_json(self):
        """
        Returns a string containing a JSON representation of this Genome.
        """
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_string):
        """
        Reconstruct an instance from a JSON string.
        """
        state_dict = json.loads(json_string)
        return cls.from_dict(state_dict)

    def _to_pairs(self):
        return tuple(sorted(self.to_dict.items()))

    def __hash__(self):
        return hash(self._to_pairs())

    def __getstate__(self):
        return self.to_dict()

    def __setstate__(self, state_dict):
        state_dict = self._initialize_nested_objects_in_state_dict(state_dict)
        self.__init__(**state_dict)

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

    Relies on the following condition:
         (1) a user-defined to_dict method
         (2) the keys of to_dict() must match the arguments to __init__
    """

    def __str__(self):
        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join("%s=%s" % (k, v) for (k, v) in self.to_dict().items()))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.__class__ is other.__class__ and self.to_dict() == other.to_dict()

    def to_dict(self):
        """
        Derived classes must implement this method and return a dictionary
        whose keys match the parameters to __init__. The values must be
        primitive types (string, int, float, tuple, dict, list).
        """
        raise NotImplementedError("Method to_dict() not implemented for %s" % (
            self.__class__.__name__,))

    @classmethod
    def _reconstruct_nested_objects(cls, state_dict):
        """
        Nested serializable objects will be represented as dictionaries so we
        allow manual reconstruction of those objects in this method.

        By default just returns the state dictionary unmodified.
        """
        return state_dict

    @classmethod
    def from_dict(cls, state_dict):
        """
        Given a dictionary of flattened fields (result of calling to_dict()),
        returns an instance.
        """
        state_dict = cls._reconstruct_nested_objects(state_dict)
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

    def write_json_file(self, path):
        """
        Serialize this VariantCollection to a JSON representation and write it
        out to a text file.
        """
        with open(path, "w") as f:
            f.write(self.to_json())

    @classmethod
    def read_json_file(cls, path):
        """
        Construct a VariantCollection from a JSON file.
        """
        with open(path, 'r') as f:
            json_string = f.read()
        return cls.from_json(json_string)

    def __hash__(self):
        return hash(tuple(sorted(self.to_dict().items())))

    def __reduce__(self):
        """
        Overriding this method directs the default pickler to reconstruct
        this object using our from_dict method.
        """

        # Due to differences between Python 2.7 and Python 3.x I have to go
        # through some acrobatics to return a pickle-able function
        return self.from_dict.__func__, (self.__class__, self.to_dict(),)

def class_from_serializable_representation(class_repr):
    """
    Given the name of a module and a class it contains, imports that module
    and gets the class object from it.
    """
    module_string, class_name = class_repr
    module = __import__(module_string)
    return getattr(module, class_name)

def class_to_serializable_representation(cls):
    """
    Given a class, return two strings:
        - fully qualified import path for its module
        - name of the class

    The class can be reconstructed from these two strings by calling
    class_from_serializable_representation.
    """
    return (cls.__module__, cls.__name__)

def function_from_serializable_representation(fn_repr):
    """
    Given the name of a module and a function it contains, imports that module
    and gets the class object from it.
    """
    if fn_repr is None:
        return None

    module_string, class_name = fn_repr
    module = __import__(module_string)
    return getattr(module, class_name)


def function_to_serializable_representation(fn):
    """
    Given a function, return two strings:
        - fully qualified import path for its module
        - name of the class

    The original function can be reconstructed from these two strings by
    calling class_from_serializable_representation. Function with closure
    variables will not work. If fn is None then None is returned as its
    representation.
    """
    if fn is None:
        return None

    if fn.__closure__ is not None:
        raise ValueError("No serializable representation for closure %s" % (fn,))
    return (fn.__module__, fn.__name__)

def object_to_serializable_representation(obj):
    """
    Given an instance of a Python object, returns a tuple whose
    first element is a primitive representation of the class and whose
    second element is a dictionary of instance data.
    """
    if obj is None:
        return None

    if not hasattr(obj, 'to_dict'):
        raise ValueError("Expected %s to have method to_dict()" % (obj,))

    state_dict = obj.to_dict()
    class_representation = class_to_serializable_representation(obj.__class__)
    return (class_representation, state_dict)


def object_from_serializable_representation(obj_repr):
    """
    Given a primitive representation of some object, reconstructs
    the class from its module and class names and then instantiates. Returns
    instance object.
    """
    if obj_repr is None:
        return None

    class_repr, state_dict = obj_repr
    subclass = class_from_serializable_representation(class_repr)
    return subclass.from_dict(state_dict)

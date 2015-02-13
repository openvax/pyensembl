"""
Routines for testing or asserting that a value is certain type, in a way that
works in both Python 2.7 and Python 3.4.
"""
from __future__ import print_function, division, absolute_import

# In Python 2, "unicode" is a type. In Python 3, the str type replaces it.
try:
    unicode
except NameError:
    # Python 3
    unicode = str

# Similarly, in Python 2 "long" is a type. In Python 3, int replaces it.
try:
    long
except NameError:
    # Python 3
    long = int

def is_string(obj):
    """
    Is the object an instance of str or unicde?
    """
    return isinstance(obj, (str, unicode))
    
def is_integer(obj):
    """
    Is the object an instance of int or long?
    """
    return isinstance(obj, (int, long))

def assert_string(obj, name=None, nonempty=False):
    """
    Raise an exception if the obj is not of type str or unicode.
    
    If name is provided it is used in the exception message.

    If nonempty=True, then an exception is raised if the object is the empty
    string.
    """
    if not is_string(obj):
        raise TypeError(
            (("%s: " % name) if name else "") +
            "expected string, got: '%s' of type '%s'"
                % (name, str(obj), type(obj)))
    if nonempty and not obj:
        raise ValueError(
            (("%s: " % name) if name else "") +
            "string must be nonempty.")

def assert_integer(obj, name=None):
    """
    Raise an exception if the obj is not of type int or long.
    
    If name is provided it is used in the exception message.
    """
    if not is_integer(obj):
        raise TypeError(
            (("%s: " % name) if name else "") +
            "expected int, got: '%s' of type '%s'"
                % (name, str(obj), type(obj)))



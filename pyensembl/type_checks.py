"""
Routines for testing or asserting that a value is certain type, in a way that
works in both Python 2.7 and Python 3.4.
"""
from __future__ import print_function, division, absolute_import

# In Python 3, the "unicode" and "long" types went away.
string_types = (str,)  # Python 2 and 3.
try: string_types += (unicode,) # Python 2.
except NameError: pass  # Python 3.

integer_types = (int,)  # Python 2 and 3.
try: integer_types += (long,)  # Python 2.
except NameError: pass  # Python 3.

def is_string(obj):
    """
    Is the object a string?
    """
    return isinstance(obj, string_types)
    
def is_integer(obj):
    """
    Is the object an integer?
    """
    return isinstance(obj, integer_types)

def require_string(obj, name=None, nonempty=False):
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
                % (str(obj), type(obj)))
    if nonempty and not obj:
        raise ValueError(
            (("%s: " % name) if name else "") +
            "string must be nonempty.")

def require_integer(obj, name=None):
    """
    Raise an exception if the obj is not of type int or long.
    
    If name is provided it is used in the exception message.
    """
    if not is_integer(obj):
        raise TypeError(
            (("%s: " % name) if name else "") +
            "expected int, got: '%s' of type '%s'"
                % (str(obj), type(obj)))



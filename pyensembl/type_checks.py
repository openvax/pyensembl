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
    require_instance(obj, string_types, name, "string")
    if nonempty and not obj:
        raise ValueError(
            (("%s: " % name) if name else "") +
            "string must be nonempty.")

def require_integer(obj, name=None):
    """
    Raise an exception if the obj is not of type int or long.
    
    If name is provided it is used in the exception message.
    """
    require_instance(obj, integer_types, name, "integer")

def require_instance(obj, types=None, name=None, type_name=None, truncate_at=80):
    """
    Raise an exception if obj is not an instance of one of the specified types.
    
    Similarly to isinstance, 'types' may be either a single type or a tuple of
    types.

    If name or type_name is provided, it is used in the exception message.
    The object's string representation is also included in the message,
    truncated to 'truncate_at' number of characters.
    """
    if not isinstance(obj, types):
        obj_string = str(obj)
        if len(obj_string) > truncate_at:
            obj_string = obj_string[:truncate_at - 3] + "..."
        if type_name is None:
            try:
                type_name = "one of " + ", ".join(str(t) for t in types)
            except TypeError:
                type_name = str(types)
        raise TypeError(
            (("%s: " % name) if name else "") +
            "expected %s. Got: '%s' of type '%s'"
                % (type_name, obj_string, type(obj)))

def require_iterable_of(objs, types, name=None, type_name=None, truncate_at=80):
    """
    Raise an exception if objs is not an iterable with each element an instance
    of one of the specified types.

    See `require_instance` for descriptions of the other parameters.
    """
    # Fast pass for common case where all types are correct.
    # This avoids the more expensive loop below. A typical speedup from this
    # optimization is 6.6 sec -> 1.7 sec, for testing a list of size 10,000,000.
    try:
        if all(isinstance(obj, types) for obj in objs):
            return
    except TypeError:
        # We don't require that objs is a list in this function, just that it's
        # iterable. We specify 'list' below as a convenient way to throw the
        # desired error.
        require_instance(objs, list, name, "iterable", truncate_at)

    # Some type isn't correct. We reuse the require_instance function to raise
    # the exception.
    prefix = ("%s: " % name) if name else ""
    for (i, obj) in enumerate(objs):
        element_name = prefix + ("element at index %d" % i)
        require_instance(obj, types, element_name, type_name, truncate_at)
    assert False, "Shouldn't reach here."


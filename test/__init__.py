"""
Utility functions for tests.
"""

import os

def data_path(name):
    """
    Return the absolute path to a file in the varcode/test/data directory.
    The name specified should be relative to varcode/test/data.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)

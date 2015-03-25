# Copyright (c) 2015. Mount Sinai School of Medicine
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

from functools import wraps
from typechecks import is_string

CACHE_SUBDIR = "ensembl"

def memoize(fn):
    """Simple memoization decorator for functions and methods,
    assumes that all arguments to the function can be hashed and
    compared.
    """
    memoized_values = {}

    @wraps(fn)
    def wrapped_fn(*args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        if key not in memoized_values:
            memoized_values[key] = fn(*args, **kwargs)
        return memoized_values[key]

    return wrapped_fn

def is_valid_ensembl_id(ensembl_id):
    """Is the argument a valid ID for any Ensembl feature?"""
    return is_string(ensembl_id) and ensembl_id.startswith("ENS")

def require_ensembl_id(ensembl_id):
    if not is_valid_ensembl_id(ensembl_id):
        raise ValueError("Invalid Ensembl ID '%s'" % ensembl_id)

def is_valid_human_transcript_id(transcript_id):
    """Is the argument a valid identifier for human Ensembl transcripts?"""
    return is_string(transcript_id) and transcript_id.startswith("ENST")

def require_human_transcript_id(transcript_id):
    if not is_valid_human_transcript_id(transcript_id):
        raise ValueError("Invalid transcript ID '%s'" % transcript_id)

def is_valid_human_protein_id(protein_id):
    """Is the argument a valid identifier for human Ensembl proteins?"""
    return is_string(protein_id) and protein_id.startswith("ENSP")

def require_human_protein_id(protein_id):
    if not is_valid_human_protein_id(protein_id):
        raise ValueError("Invalid protein ID '%s'" % protein_id)

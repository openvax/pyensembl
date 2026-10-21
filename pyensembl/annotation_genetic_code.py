"""Read sparse transcript/CDS genetic-code declarations without cache writes."""

from contextlib import closing
import gzip
import json
from os import fspath
from os.path import join, isfile, abspath
import re

import datacache

from .database import (
    DATABASE_SCHEMA_VERSION, GENETIC_CODE_SOURCE_TABLE, _file_state, annotation_db_filename,
)
from .genetic_code import translation_table


_COLUMNS = ("transl_table", "translation_table")
_DECLARATION = re.compile(r"(?:^|;)\s*(?:transl_table|translation_table)\b")
_ATTRIBUTE = re.compile(
    r'''(?:^|;)\s*([^\s;]+)[ \t]+(?:"([^"\n]*)"|'([^'\n]*)'|([^;\n]*))\s*(?=;|$)'''
)


def _add(result, identifier, values):
    # gtfparse joins repeated attributes with commas. Retain every declaration
    # so agreement/conflicts have identical semantics in new and old indexes.
    values = [item.strip() for value in values if value is not None
              for item in str(value).split(",") if item.strip()]
    if values:
        if not identifier:
            raise ValueError("Translation-table declaration has no transcript_id")
        result.setdefault(identifier, set()).update(values)


class AnnotationGeneticCodes:
    """One in-memory sparse map, invalidated when the source/index changes.

    Complete indexes with a matching source fingerprint avoid reading the
    GTF body. Others use a streaming scan, never a reindex or sidecar.
    Transcript IDs here are their stored forms, already resolved by Genome.
    """

    def __init__(self, path, source, cache_directory):
        self.path = fspath(path)
        source = fspath(source)
        self.source = source if "://" in source else abspath(source)
        self.index = join(cache_directory, annotation_db_filename(self.path, self.source))
        self._state = None
        self._values = {}

    def _file_state(self):
        try:
            source = _file_state(self.path)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                "Cannot check source translation-table declarations: installed GTF is missing: " + self.path
            ) from error
        return source, _file_state(self.index) if isfile(self.index) else None

    def _read_index(self):
        if not isfile(self.index):
            return None
        connection = datacache.connect_if_correct_version(
            self.index, DATABASE_SCHEMA_VERSION, read_only=True,
        )
        if connection is None:
            return None
        with closing(connection):
            marker = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                (GENETIC_CODE_SOURCE_TABLE,),
            ).fetchone()
            if marker is None:
                return None
            row = connection.execute('SELECT state FROM "' + GENETIC_CODE_SOURCE_TABLE + '"').fetchone()
            if row is None or json.loads(row[0]) != _file_state(self.path):
                return None
            tables = []
            for feature in ("transcript", "CDS"):
                columns = {row[1] for row in connection.execute('PRAGMA table_info("%s")' % feature)}
                if not columns:
                    continue
                if not {"transcript_id", *_COLUMNS}.issubset(columns):
                    return None
                tables.append(feature)
            if not tables:
                return None
            result = {}
            for feature in tables:
                for row in connection.execute(
                    'SELECT transcript_id, transl_table, translation_table FROM "%s" '
                    "WHERE COALESCE(transl_table, '') != '' OR COALESCE(translation_table, '') != ''" % feature
                ):
                    _add(result, row[0], row[1:])
            return result

    def _read_gtf(self):
        with open(self.path, "rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
        opener = gzip.open if compressed else open
        result = {}
        with opener(self.path, "rt", encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t", 8)
                if len(fields) != 9 or fields[2] not in ("transcript", "CDS"):
                    continue
                attributes = fields[8]
                if not any(column in attributes for column in _COLUMNS):
                    continue
                parsed = list(_ATTRIBUTE.finditer(attributes))
                def value(match):
                    return next(group for group in match.groups()[1:] if group is not None)
                identifiers = {value(match) for match in parsed if match[1] == "transcript_id"}
                values = [value(match) for match in parsed if match[1] in _COLUMNS]
                if _DECLARATION.search(_ATTRIBUTE.sub("", attributes)) or (values and len(identifiers) != 1):
                    raise ValueError("Malformed translation-table declaration at %s:%d" % (self.path, number))
                if values:
                    _add(result, identifiers.pop(), values)
        return result

    def get(self, transcript_id):
        state = self._file_state()
        if state != self._state:
            values = self._read_index()
            if values is None:
                values = self._read_gtf()
            if self._file_state() != state:
                raise OSError("Annotation changed while reading translation-table declarations")
            self._values, self._state = values, state
        codes = set()
        for value in self._values.get(transcript_id, ()):
            if not re.fullmatch(r"[1-9][0-9]*", value):
                raise ValueError("Invalid translation table for transcript %s: %r" % (transcript_id, value))
            codes.add(translation_table(int(value)))
        if len(codes) > 1:
            raise ValueError("Conflicting translation tables for transcript %s: %s" % (transcript_id, sorted(codes)))
        return next(iter(codes), None)

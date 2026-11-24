"""Read sparse transcript/CDS genetic-code declarations without cache writes."""

from contextlib import closing
import json
from os import fspath
from os.path import join, isfile, abspath
import re

import datacache

from .common import file_fingerprint
from .database import (
    DATABASE_SCHEMA_VERSION, GENETIC_CODE_SOURCE_TABLE, annotation_db_filename,
)
from .genetic_code import translation_table


_COLUMNS = ("transl_table", "translation_table")


def add_translation_table_declarations(result, identifier, values):
    """Record raw GTF translation-table attribute values for a transcript."""
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

    def file_state(self):
        """Return ``(gtf_fingerprint, index_fingerprint)`` for cache invalidation."""
        try:
            source = file_fingerprint(self.path)
        except FileNotFoundError as error:
            raise FileNotFoundError(
                "Cannot check source translation-table declarations: installed GTF is missing: " + self.path
            ) from error
        return source, file_fingerprint(self.index) if isfile(self.index) else None

    def read_index(self):
        """Read sparse transcript/CDS translation-table declarations from SQLite, or None if unavailable."""
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
            if row is None or json.loads(row[0]) != file_fingerprint(self.path):
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
                    add_translation_table_declarations(result, row[0], row[1:])
            return result

    def read_gtf(self):
        """Read transcript/CDS translation-table declarations directly from the GTF via gtfparse."""
        # Imported lazily so plain `import pyensembl` does not load pandas/PyArrow.
        from gtfparse import read_gtf

        frame = read_gtf(
            self.path,
            result_type="pandas",
            features={"transcript", "CDS"},
            usecols=["transcript_id", *_COLUMNS],
            cast_version_columns=False,
        )
        result = {}
        for row in frame.to_dict("records"):
            identifier = row.get("transcript_id")
            values = [row.get(column) for column in _COLUMNS]
            if any(v is not None and str(v).strip() for v in values):
                if not identifier or "," in str(identifier):
                    raise ValueError("Malformed translation-table declaration in %s" % self.path)
            add_translation_table_declarations(result, identifier, values)
        return result

    def get(self, transcript_id):
        state = self.file_state()
        if state != self._state:
            values = self.read_index()
            if values is None:
                values = self.read_gtf()
            if self.file_state() != state:
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

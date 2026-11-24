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

"""Logging configuration for the ``pyensembl`` command-line tool."""

import logging


class ConciseFormatter(logging.Formatter):
    """Plain progress formatter that prefixes only warnings and errors."""

    def format(self, record):
        """Format ``record``, prepending the lowercase level name for warnings and above.

        Parameters
        ----------
        record : logging.LogRecord
            Log record to format.

        Returns
        -------
        str
            Formatted log line.
        """
        message = super().format(record)
        if record.levelno >= logging.WARNING:
            return "%s: %s" % (record.levelname.lower(), message)
        return message


CLI_LOGGING_STATE = {"handler": None}


def reset_cli_logging():
    """Clear the remembered CLI logging handler reference.

    Useful in tests that snapshot and restore logger handlers directly.
    """
    CLI_LOGGING_STATE["handler"] = None


def configure_logging(verbose=False):
    """Configure concise stderr progress logging for ``pyensembl`` and ``datacache``.

    Only the CLI entry point calls this function, so importing ``pyensembl`` as
    a library never alters the host application's logging configuration.

    Parameters
    ----------
    verbose : bool, optional
        When ``True``, emit ``DEBUG`` logs with timestamps and logger names.
        When ``False`` (default), emit plain ``INFO`` lines for ``pyensembl``
        and ``WARNING`` lines for ``datacache``.
    """
    previous = CLI_LOGGING_STATE["handler"]
    handler = logging.StreamHandler()  # stderr; stdout carries results
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
        if verbose
        else ConciseFormatter()
    )
    for name, level in (
        ("pyensembl", logging.DEBUG if verbose else logging.INFO),
        # datacache reports each download step and SQL statement.
        ("datacache", logging.DEBUG if verbose else logging.WARNING),
    ):
        package_logger = logging.getLogger(name)
        if previous is not None:
            package_logger.removeHandler(previous)
        package_logger.setLevel(level)
        package_logger.addHandler(handler)
        # Our handler prints; a host application's root handler must not too.
        package_logger.propagate = False
    CLI_LOGGING_STATE["handler"] = handler

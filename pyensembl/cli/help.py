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

"""Readable argparse help and table styling with optional terminal color."""

import argparse
import os
import re
import shutil
import sys

ANSI_RESET = "\033[0m"
ANSI_BOLD = "\033[1m"
ANSI_DIM = "\033[2m"
ANSI_YELLOW = "\033[33m"
ANSI_CYAN = "\033[36m"
ANSI_BOLD_CYAN = "\033[1;36m"

OPTION_PATTERN = re.compile(r"(?<![\w-])(--[a-zA-Z][\w-]*|-[a-zA-Z])\b")

DIM_STATUS_VALUES = frozenset({
    "-",
    "?",
    "absent",
    "n/a",
    "none",
    "not applicable",
    "not checked",
    "not configured",
    "unrecorded",
    "unsupported",
})

WARNING_STATUS_VALUES = frozenset({
    "changed",
    "corrupt",
    "corrupted",
    "error",
    "inaccessible",
    "incomplete",
    "invalid",
    "invalid index",
    "invalid reference",
    "invalid-receipt",
    "mismatch",
    "missing",
    "needs index",
    "not indexed",
    "stale",
})


def color_enabled(stream):
    """Return whether ``stream`` supports optional ANSI terminal styling.

    Styling is enabled only when ``stream`` is an interactive TTY and neither
    ``NO_COLOR`` nor ``TERM=dumb`` is set in the environment.

    Parameters
    ----------
    stream : file-like or None
        Destination output stream (for example ``sys.stdout``).

    Returns
    -------
    bool
        ``True`` if ANSI styling should be applied, ``False`` otherwise.
    """
    return bool(
        getattr(stream, "isatty", lambda: False)()
        and not os.environ.get("NO_COLOR")
        and os.environ.get("TERM") != "dumb"
    )


def style_help(text, stream):
    """Style argparse help text when writing to a color-capable terminal.

    Applies one subtle cyan hue: section headings in bold cyan, option flags
    in cyan, and ``usage:`` in bold. Descriptions, help text, metavars, and
    defaults remain plain.

    Parameters
    ----------
    text : str
        Plain help text returned by ``ArgumentParser.format_help()``.
    stream : file-like or None
        Destination output stream checked with :func:`color_enabled`.

    Returns
    -------
    str
        Styled help text for color-capable terminals, or ``text`` unchanged.
    """
    if not color_enabled(stream):
        return text

    lines = []
    for line in text.splitlines(keepends=True):
        if line.rstrip().endswith(":") and not line[0].isspace():
            lines.append(ANSI_BOLD_CYAN + line.rstrip("\n") + ANSI_RESET + "\n")
        else:
            line = OPTION_PATTERN.sub(
                lambda match: ANSI_CYAN + match[0] + ANSI_RESET, line
            )
            if line.startswith("usage:"):
                line = line.replace("usage:", ANSI_BOLD + "usage:" + ANSI_RESET, 1)
            lines.append(line)
    return "".join(lines)


def style_table_header(text, stream=None, use_color=None):
    """Style the first line of an aligned table in bold cyan on color terminals.

    Parameters
    ----------
    text : str
        Aligned table string (single-line header or multi-line table).
    stream : file-like, optional
        Destination stream checked when ``use_color`` is ``None``. Defaults
        to ``sys.stdout``.
    use_color : bool, optional
        Explicitly enable or disable ANSI styling, overriding TTY detection.

    Returns
    -------
    str
        ``text`` with its first line wrapped in bold cyan when styling is
        enabled, or ``text`` unchanged.
    """
    enabled = (
        color_enabled(sys.stdout if stream is None else stream)
        if use_color is None
        else bool(use_color)
    )
    if not text or not enabled:
        return text
    lines = text.splitlines(keepends=True)
    first = lines[0]
    ending = "\n" if first.endswith("\n") else ""
    lines[0] = ANSI_BOLD_CYAN + first.rstrip("\n") + ANSI_RESET + ending
    return "".join(lines)


def classify_status(text):
    """Classify a table status cell as ``"dim"``, ``"warning"``, or ``"plain"``.

    Absent / not-applicable values map to ``"dim"``, warning states (including
    compound DNA strings such as ``"local dna.fa, missing"``) map to
    ``"warning"``, and normal values (including ``"OK"``, ``"indexed"``, and
    ``"verified"``) map to ``"plain"``.

    Parameters
    ----------
    text : object
        Status cell value.

    Returns
    -------
    str
        One of ``"dim"``, ``"warning"``, or ``"plain"``.
    """
    value = str(text).strip().lower()
    if value in DIM_STATUS_VALUES:
        return "dim"
    if value in WARNING_STATUS_VALUES:
        return "warning"
    if "," in value and value.rsplit(",", 1)[-1].strip() in WARNING_STATUS_VALUES:
        return "warning"
    return "plain"


def style_status(text, stream=None, level=None, use_color=None):
    """Style a table status value using dim or yellow only when color is enabled.

    Uses dim (``\\033[2m``) for absent / not-applicable values and yellow
    (``\\033[33m``) for warnings. Normal values, including ``"OK"``, stay plain
    (never green or red).

    Parameters
    ----------
    text : object
        Status value to format.
    stream : file-like, optional
        Destination stream checked when ``use_color`` is ``None``. Defaults
        to ``sys.stdout``.
    level : str, optional
        Explicit status category: ``"dim"`` (or ``"absent"``), ``"warning"``,
        or ``"plain"`` (or ``"normal"`` / ``"ok"``). When ``None`` (default),
        inferred via :func:`classify_status`.
    use_color : bool, optional
        Explicitly enable or disable ANSI styling, overriding TTY detection.

    Returns
    -------
    str
        Styled status string when color is enabled and the status is dim or a
        warning; plain string otherwise.

    Raises
    ------
    ValueError
        If ``level`` is not a recognized status category.
    """
    rendered = str(text)
    if level is None:
        category = classify_status(rendered)
    else:
        normalized = str(level).strip().lower()
        if normalized in ("dim", "absent"):
            category = "dim"
        elif normalized == "warning":
            category = "warning"
        elif normalized in ("plain", "normal", "ok"):
            category = "plain"
        else:
            raise ValueError("Unknown status level: %r" % (level,))

    enabled = (
        color_enabled(sys.stdout if stream is None else stream)
        if use_color is None
        else bool(use_color)
    )
    if not enabled or not rendered or category == "plain":
        return rendered
    code = ANSI_DIM if category == "dim" else ANSI_YELLOW
    return code + rendered + ANSI_RESET


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Limit reading width and separate options without adding table borders."""

    def __init__(self, prog):
        """Initialize the formatter with a bounded terminal width.

        Parameters
        ----------
        prog : str
            Program name displayed in usage output.
        """
        width = max(40, min(100, shutil.get_terminal_size().columns - 2))
        super().__init__(prog, max_help_position=32, width=width)

    def _split_lines(self, text, width):
        # RST literals are useful in generated docs but noisy in terminal help.
        return super()._split_lines(text.replace("``", ""), width)

    def _format_action(self, action):
        text = super()._format_action(action)
        return text + "\n" if text else text


class HelpArgumentParser(argparse.ArgumentParser):
    """Keep ``format_help()`` plain; add color only when printing to a terminal."""

    def __init__(self, *args, **kwargs):
        """Construct the parser using :class:`HelpFormatter` and uniform color rules.

        On Python 3.14+, disables argparse's built-in multi-hue theme so all
        supported Python versions share the same subtle cyan styling.
        """
        kwargs.setdefault("formatter_class", HelpFormatter)
        if sys.version_info >= (3, 14):
            # Apply the same styling on all supported Python versions.
            kwargs["color"] = False
        super().__init__(*args, **kwargs)

    def print_help(self, file=None):
        """Print help text, styling headings and flags on color-capable terminals.

        Parameters
        ----------
        file : file-like, optional
            Destination stream; defaults to ``sys.stdout``.
        """
        stream = sys.stdout if file is None else file
        self._print_message(style_help(self.format_help(), stream), stream)

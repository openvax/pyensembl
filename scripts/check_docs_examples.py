"""Run tutorial Python blocks and compare their actual output with the page.

Requires human Ensembl release 93 installed. Does not download missing data.
"""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pyensembl import EnsemblRelease  # noqa: E402


def check_examples():
    with EnsemblRelease(93, species="human") as data:
        if not data.installed():
            raise SystemExit("First run: pyensembl install --release 93 --species human")
    executed = 0
    for filename in ("docs/index.md", "docs/guides/features.md"):
        text = (ROOT / filename).read_text()
        blocks = re.findall(r"^```([^\n]*)\n(.*?)^```[ \t]*$", text,
                            re.MULTILINE | re.DOTALL)
        examples = [(code, output)
                    for (language, code), (next_language, output)
                    in zip(blocks, blocks[1:])
                    if language.strip() == "python" and next_language.strip() == "text"]
        if not examples:
            raise SystemExit("No checked examples found in %s" % filename)
        namespace = {}
        try:
            for code, expected in examples:
                stream = StringIO()
                with redirect_stdout(stream):
                    exec(compile(code, filename, "exec"), namespace)
                if stream.getvalue().strip() != expected.strip():
                    raise SystemExit("%s output differs:\nexpected:\n%s\nactual:\n%s" %
                                     (filename, expected, stream.getvalue()))
                executed += 1
        finally:
            if "data" in namespace:
                namespace["data"].close()
    print("Executed %d examples against human GRCh38 / Ensembl 93" % executed)


if __name__ == "__main__":
    check_examples()

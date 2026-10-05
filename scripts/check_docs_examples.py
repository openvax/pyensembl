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
    text = (ROOT / "docs/getting-started.md").read_text()
    examples = re.findall(r"```python\n(.*?)```\s*```text\n(.*?)```", text, re.DOTALL)
    if not examples:
        raise SystemExit("No tutorial examples found")
    namespace = {}
    try:
        for code, expected in examples:
            stream = StringIO()
            with redirect_stdout(stream):
                exec(compile(code, "docs/getting-started.md", "exec"), namespace)
            if stream.getvalue().strip() != expected.strip():
                raise SystemExit("Tutorial output differs:\nexpected:\n%s\nactual:\n%s" %
                                 (expected, stream.getvalue()))
    finally:
        if "data" in namespace:
            namespace["data"].close()
    print("Executed %d tutorial examples against human GRCh38 / Ensembl 93" % len(examples))


if __name__ == "__main__":
    check_examples()

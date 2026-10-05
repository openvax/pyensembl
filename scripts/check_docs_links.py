"""Validate rendered links/anchors, README destinations and public API coverage.

Run after `mkdocs build --strict` from the repository root. No network required.
"""

import ast
from html.parser import HTMLParser
import inspect
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit

import markdown

# Direct execution places scripts/ first; inspect this checkout's interfaces.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pyensembl  # noqa: E402


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        if tag == "a":
            if "name" in attrs:
                self.ids.add(attrs["name"])
            if "href" in attrs:
                self.links.append(attrs["href"])


def check_links():
    site = ROOT / "site"
    if not (site / "index.html").exists():
        raise SystemExit("Build the site first: mkdocs build --strict")
    pages = {path.resolve(): Page(path.read_text()) for path in site.rglob("*.html")}
    failures = []
    links = 0
    for path, page in pages.items():
        for href in page.links:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            if url.path.startswith("/"):
                target = site / unquote(url.path).removeprefix("/pyensembl/").lstrip("/")
            else:
                target = path.parent / unquote(url.path) if url.path else path
            if target.is_dir():
                target /= "index.html"
            target = target.resolve()
            links += 1
            if not target.is_file():
                failures.append("%s: missing %s" % (path.relative_to(site), href))
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                failures.append("%s: missing anchor %s" % (path.relative_to(site), href))
    readme = Page(markdown.markdown((ROOT / "README.md").read_text(),
                                   extensions=["tables", "toc", "fenced_code"]))
    for href in readme.links:
        url = urlsplit(href)
        repo_prefix = "/openvax/pyensembl/blob/main/"
        if url.netloc == "github.com" and url.path.startswith(repo_prefix):
            target = ROOT / unquote(url.path.removeprefix(repo_prefix))
        elif not url.scheme and not url.netloc and url.path:
            target = ROOT / unquote(url.path)
        else:
            continue
        if not target.is_file():
            failures.append("README: missing %s" % href)
    # Signatures and own public methods must actually appear in the rendered
    # reference. This catches decorators or doc generation omitting interfaces.
    reference_ids = set().union(*(page.ids for path, page in pages.items()
                                 if "reference" in path.parts))
    for name in set(pyensembl.__all__):
        value = getattr(pyensembl, name)
        if not (inspect.isclass(value) or inspect.isfunction(value)):
            continue
        qualified = "pyensembl." + name
        if qualified not in reference_ids:
            failures.append("API reference omits %s" % qualified)
        if inspect.isclass(value):
            tree = ast.parse(inspect.getsource(value))
            for method in tree.body[0].body:
                if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef)) and not method.name.startswith("_"):
                    if qualified + "." + method.name not in reference_ids:
                        failures.append("API reference omits %s.%s" % (qualified, method.name))
    if failures:
        raise SystemExit("\n".join(failures))
    print("Checked %d pages, %d internal links, README destinations and public API coverage" %
          (len(pages), links))


if __name__ == "__main__":
    check_links()

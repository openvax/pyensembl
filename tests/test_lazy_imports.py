"""``import pyensembl`` and the CLI module stay cheap and do not read the cache (#510, #514)."""

import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pyensembl

REPOSITORY = Path(pyensembl.__file__).resolve().parents[1]

# Modules that make an import slow. Importing the package or the command-line
# module must load none of them; commands import them when they run.
HEAVY_MODULES = (
    "pandas", "numpy", "pyarrow", "polars", "gtfparse", "pyfaidx", "sqlite3",
    "datacache", "serializable", "urllib.request", "http.client",
    "pyensembl.database", "pyensembl.genome", "pyensembl.ensembl_release",
    "pyensembl.download_cache", "pyensembl.reference_name", "pyensembl.species",
    "pyensembl.sequence_data", "pyensembl.genome_fasta", "pyensembl.genome_fasta_cache",
    "pyensembl.dated_releases", "pyensembl.release_catalogue",
)


def run_python(code, cache):
    """Run Python code in a fresh interpreter that imports this checkout.

    Parameters
    ----------
    code : str
        Source passed to ``python -c``.
    cache : pathlib.Path
        Value for ``PYENSEMBL_CACHE_DIR``.

    Returns
    -------
    str
        What the code printed.
    """
    environment = dict(os.environ, PYENSEMBL_CACHE_DIR=str(cache), PYTHONPATH=str(REPOSITORY))
    result = subprocess.run([sys.executable, "-c", code], env=environment,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_importing_package_and_cli_loads_no_heavy_modules(tmp_path):
    loaded = run_python(
        "import json, sys, pyensembl, pyensembl.shell; "
        "print(json.dumps(sorted(m for m in %r if m in sys.modules)))" % (HEAVY_MODULES,),
        tmp_path,
    )
    assert json.loads(loaded) == []


def test_import_does_not_read_the_genome_cache(tmp_path):
    (tmp_path / "pyensembl" / "GRCh38" / "ensembl93").mkdir(parents=True)
    touched = run_python(textwrap.dedent("""
        import json, os, sys
        root = os.environ["PYENSEMBL_CACHE_DIR"]
        events = []

        def record_cache_reads(event, args):
            if event in ("open", "os.listdir", "os.scandir", "sqlite3.connect") and args:
                path = args[0]
                if isinstance(path, (str, bytes, os.PathLike)) and os.fsdecode(path).startswith(root):
                    events.append([event, os.fsdecode(path)])

        sys.addaudithook(record_cache_reads)
        import pyensembl, pyensembl.shell, pyensembl.reference_name
        print(json.dumps(events))
    """), tmp_path)
    assert json.loads(touched) == []


def test_convenience_genome_is_chosen_once_on_first_use(tmp_path):
    result = run_python(
        "import json, pyensembl\n"
        "import pyensembl.reference_name as reference_name\n"
        "real = reference_name.genome_for_reference_name\n"
        "calls = []\n"
        "def counting(name):\n"
        "    calls.append(name)\n"
        "    return real(name)\n"
        "reference_name.genome_for_reference_name = counting\n"
        "first = pyensembl.ensembl_grch38\n"
        "from pyensembl import ensembl_grch38\n"
        "from pyensembl.reference_name import ensembl_grch38 as again\n"
        "print(json.dumps(dict(\n"
        "    calls=calls, same=first is ensembl_grch38 is again,\n"
        "    expected=first == real('grch38'), release=first.release,\n"
        "    listed='ensembl_grch37' in dir(reference_name)\n"
        ")))\n",
        tmp_path,
    )
    # With an empty cache the newest supported GRCh38 release is chosen.
    assert json.loads(result) == dict(
        calls=["grch38"], same=True, expected=True,
        release=pyensembl.MAX_ENSEMBL_RELEASE, listed=True,
    )


def test_public_names_resolve_and_are_listed(tmp_path):
    result = run_python(
        "import json, pyensembl\n"
        "unlisted = [name for name in pyensembl.__all__ if name not in dir(pyensembl)]\n"
        "missing = [name for name in pyensembl.__all__ if getattr(pyensembl, name, None) is None]\n"
        "namespace = {}\n"
        "exec('from pyensembl import *', namespace)\n"
        "print(json.dumps(dict(unlisted=unlisted, missing=missing,\n"
        "                      star=sorted(set(pyensembl.__all__) - set(namespace)))))\n",
        tmp_path,
    )
    assert json.loads(result) == dict(unlisted=[], missing=[], star=[])


def test_submodules_are_attributes_after_a_bare_import(tmp_path):
    # Before 2.33 the package imported most submodules; code still uses them so.
    result = run_python(
        "import pyensembl\n"
        "print(pyensembl.species.Species.__name__, pyensembl.download_cache.DownloadCache.__name__)\n"
        "try:\n"
        "    pyensembl.no_such_name\n"
        "except AttributeError as error:\n"
        "    print(error)\n",
        tmp_path,
    )
    assert result.splitlines() == [
        "Species DownloadCache", "module 'pyensembl' has no attribute 'no_such_name'",
    ]


def test_public_objects_come_from_their_defining_modules():
    assert pyensembl.EnsemblRelease.__module__ == "pyensembl.ensembl_release"
    assert pyensembl.Genome.__module__ == "pyensembl.genome"
    assert pyensembl.cached_release is pyensembl.ensembl_release.cached_release
    assert pyensembl.MAX_ENSEMBL_RELEASE == pyensembl.ensembl_versions.MAX_ENSEMBL_RELEASE

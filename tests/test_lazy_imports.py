"""``import pyensembl`` and the CLI module stay cheap and do not read the cache (#510, #514)."""

import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pyensembl

REPOSITORY = Path(pyensembl.__file__).resolve().parents[1]

# Modules that make an import slow. Importing the package must load none of
# them (including ``datacache``); importing the CLI modules may load the
# stdlib-only ``datacache.cli_help`` helper (and lazy ``datacache``), but must
# load no heavy runtime or download/database submodules.
HEAVY_MODULES = (
    "pandas", "numpy", "pyarrow", "polars", "gtfparse", "pyfaidx", "sqlite3",
    "datacache.download", "datacache.cache", "datacache.database_helpers",
    "requests", "filelock", "serializable", "urllib.request", "http.client",
    "pyensembl.database", "pyensembl.genome", "pyensembl.ensembl_release",
    "pyensembl.download_cache", "pyensembl.reference_name", "pyensembl.species",
    "pyensembl.sequence_data", "pyensembl.fasta_cache",
    "pyensembl.genome_fasta", "pyensembl.genome_fasta_cache",
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
    existing = os.environ.get("PYTHONPATH")
    pythonpath = os.pathsep.join([str(REPOSITORY), existing]) if existing else str(REPOSITORY)
    environment = dict(os.environ, PYENSEMBL_CACHE_DIR=str(cache), PYTHONPATH=pythonpath)
    result = subprocess.run([sys.executable, "-c", code], env=environment,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_importing_package_and_cli_loads_no_heavy_modules(tmp_path):
    package_loaded = run_python(
        "import json, sys, pyensembl; "
        "print(json.dumps(sorted(m for m in %r if m in sys.modules)))"
        % (HEAVY_MODULES + ("datacache",),),
        tmp_path,
    )
    assert json.loads(package_loaded) == []

    loaded = run_python(
        "import json, sys, pyensembl, pyensembl.shell, pyensembl.cli, "
        "pyensembl.cli.help, pyensembl.cli.parser; "
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


def test_import_creates_no_ensembl_releases_until_convenience_genome_accessed(tmp_path):
    """Importing pyensembl or reference_name builds zero EnsemblRelease objects (#510)."""
    (tmp_path / "pyensembl" / "GRCh38" / "ensembl93").mkdir(parents=True)
    result = run_python(textwrap.dedent("""
        import json, os, sys

        root = os.environ["PYENSEMBL_CACHE_DIR"]
        cache_events = []

        def record_cache_reads(event, args):
            if event in ("open", "os.listdir", "os.scandir", "sqlite3.connect") and args:
                path = args[0]
                if isinstance(path, (str, bytes, os.PathLike)) and os.fsdecode(path).startswith(root):
                    cache_events.append([event, os.fsdecode(path)])

        sys.addaudithook(record_cache_reads)

        import pyensembl
        assert "pyensembl.ensembl_release" not in sys.modules
        assert "pyensembl.reference_name" not in sys.modules

        # Patching EnsemblRelease.__init__ before importing reference_name proves
        # importing reference_name creates zero EnsemblRelease instances.
        import pyensembl.ensembl_release as er
        constructed = []
        real_init = er.EnsemblRelease.__init__

        def counting_init(self, *args, **kwargs):
            constructed.append((args, kwargs))
            return real_init(self, *args, **kwargs)

        er.EnsemblRelease.__init__ = counting_init

        import pyensembl.reference_name as reference_name
        real_select = reference_name.genome_for_reference_name
        reference_name.genome_for_reference_name = lambda *a, **kw: (_ for _ in ()).throw(
            AssertionError("genome_for_reference_name must not run at import")
        )
        # Non-convenience helpers work without triggering convenience genomes.
        assert reference_name.normalize_reference_name("grch38") == "GRCh38"
        before_access = dict(
            constructed=len(constructed),
            cache_size=len(er.EnsemblRelease._genome_cache),
            cache_events=list(cache_events),
            grch36_bound="ensembl_grch36" in vars(reference_name),
            grch37_bound="ensembl_grch37" in vars(reference_name),
            grch38_bound="ensembl_grch38" in vars(reference_name),
        )

        # Now restore genome_for_reference_name and access only ensembl_grch38.
        reference_name.genome_for_reference_name = real_select
        genome = pyensembl.ensembl_grch38
        after_access = dict(
            release=genome.release,
            constructed_gt_zero=len(constructed) > 0,
            grch36_bound="ensembl_grch36" in vars(reference_name),
            grch37_bound="ensembl_grch37" in vars(reference_name),
            grch38_bound="ensembl_grch38" in vars(reference_name),
        )
        print(json.dumps(dict(before=before_access, after=after_access)))
    """), tmp_path)
    data = json.loads(result)
    assert data["before"] == dict(
        constructed=0,
        cache_size=0,
        cache_events=[],
        grch36_bound=False,
        grch37_bound=False,
        grch38_bound=False,
    )
    assert data["after"] == dict(
        release=pyensembl.MAX_ENSEMBL_RELEASE,
        constructed_gt_zero=True,
        grch36_bound=False,
        grch37_bound=False,
        grch38_bound=True,
    )

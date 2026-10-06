# Contributing to PyEnsembl

[PyEnsembl](https://github.com/openvax/pyensembl) is open source software under
the Apache 2.0 license, and we welcome contributions. Contributed code is
assumed to use the same license.

## Filing issues

Check the [open issues](https://github.com/openvax/pyensembl/issues) first,
then [open a new issue](https://github.com/openvax/pyensembl/issues/new) for a
bug or feature request. Include your PyEnsembl and Python versions. If the
problem involves a particular gene, transcript or locus, name it and the
release, e.g. "Missing transcript sequence for BRCA1-002 in Ensembl release 74".

## Pull requests

- Start a new feature with an issue explaining its scope and rationale, and
  reference the issue in the PR, e.g. "Closes #123".
- Follow [PEP 8](https://peps.python.org/pep-0008/); `./lint.sh` runs ruff.
- Accompany new code with unit tests.
- Support Python 3.9 and later.
- Bump the version in `pyensembl/version.py` in every PR, including
  documentation-only changes.

## Development setup

```sh
git clone https://github.com/openvax/pyensembl.git
cd pyensembl
pip install -e '.[dev]'
./lint.sh
./test.sh
```

The `dev` extra installs pytest, pytest-cov, ruff and build. Most tests need
Ensembl data installed first; `.github/workflows/tests.yml` lists the releases
CI installs. Tests use exactly those releases, so other genomes in your cache
do not change what they check.

Species assembly ranges are checked against Ensembl's archive. After raising
`MAX_ENSEMBL_RELEASE`, recheck every assembly boundary on the live FTP servers:

```sh
PYENSEMBL_NETWORK_TESTS=1 ./test.sh tests/test_species_assemblies.py
```

Timed benchmarks are opt-in because wall-clock limits depend on the machine
and its load. Run them on an otherwise idle machine:

```sh
PYENSEMBL_BENCHMARKS=1 ./test.sh tests/test_timings.py -s
```

## Documentation

Follow the [documentation style guide](https://openvax.github.io/pyensembl/dev/documentation-style/),
then build and check the site:

```sh
pip install -e '.[dev,docs]'
./docs.sh  # strict build plus link, anchor and API checks
python scripts/check_docs_examples.py  # needs human release 93
```

Preview with `mkdocs serve` and review changed pages at desktop and narrow
widths.

## Releasing

Maintainers merge through a PR and publish from a clean main; see
[releasing](https://openvax.github.io/pyensembl/dev/releasing/).

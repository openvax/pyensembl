# Contributing to PyEnsembl

[PyEnsembl](http://www.github.com/openvax/pyensembl) is open source software and
we welcome your contributions. This document should help you get started
contributing to PyEnsembl.

## Filing Issues

If you find any bugs or problems while using PyEnsembl or have any feature requests, please feel free to file an issue against the project. When doing so, please follow the guidelines below:

To report any bugs, issues, or feature requests, please [open an issue](https://github.com/openvax/pyensembl/issues)
Please check the [current open issues](https://github.com/openvax/pyensembl/issues) to see if the request already exists
If you are filing a bug report, please describe the version of PyEnsembl and Python you are using. If your problem involves a particular gene, transcript, or genomic locus, please include that information (e.g. "Missing transcript sequence for BRCA1-002 for Ensembl release 74").

## Coding Guidelines

- PyEnsembl is written in Python and adheres to the [PEP8](https://www.python.org/dev/peps/pep-0008/)
  style guidelines.
- Contributions should come in the form of GitHub pull requests.
- New features should start with a GitHub issue explaining their scope and rationale.
- If the work is based on an existing issue, please reference the issue in the PR.
- All new code should be accompanied by comprehensive unit tests.
- If the PR fixes or implements an issue, please state "Closes #XYZ" or "Fixes #XYZ", where XYZ is the issue number.
- Please ensure that your code works under Python >= 3.9.

## Licensing

PyEnsembl is licensed under the Apache 2.0 license. Your code is assumed to be as well.

## Development setup

For development, install PyEnsembl in editable mode with development dependencies:

```sh
git clone https://github.com/openvax/pyensembl.git
cd pyensembl
pip install -e .[dev]
```

This installs the package in development mode along with tools for testing, linting, and building:
- `pytest` for running tests
- `ruff` for code linting
- `pytest-cov` for coverage reporting
- `build` for package building

Run lint and tests with:
```sh
./lint.sh
./test.sh
```

Most tests need Ensembl data installed first; `.github/workflows/tests.yml`
lists the releases CI installs. Tests use exactly those releases, so other
genomes in your cache do not change what they check.

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

Follow the [documentation style guide](documentation-style.md). Install `.[docs]` and run `./docs.sh`; with human release 93 installed, run `python scripts/check_docs_examples.py`. Review first-use and reference pages at desktop and narrow widths. Every PR includes a version bump and must pass `./lint.sh` and `./test.sh`.

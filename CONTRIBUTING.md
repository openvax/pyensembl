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

## Documentation

Follow the [documentation style guide](docs/dev/documentation-style.md): short
entry pages, installation before queries, pinned examples with interpreted
results, task guides and complete reference. Use ordinary linked names for
concepts and libraries; inline code is for literal syntax. Keep scientific
qualifications, source links and existing linked anchors.

Install `.[docs]` and run `./docs.sh`; run
`python scripts/check_docs_examples.py` with human release 93 installed. Review
first-use and reference pages at desktop and narrow widths. Run `./lint.sh`
and `./test.sh` for documentation changes too. Every PR bumps the version.

## Licensing

PyEnsembl is licensed under the Apache 2.0 license. Your code is assumed to be as well.

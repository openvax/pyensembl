# Publish a release <a id="releasing-pyensembl"></a>

Every PR includes a version bump in `pyensembl/version.py`, including
documentation changes. The Release version check compares a PR with current
main, existing tags, and PyPI. It is required on an up-to-date branch; if
another PR merges first, update the version and rerun CI before merging.
Once lint, tests and review have passed:

1. Merge the PR into main.
2. Check out main, pull the merged commit, and ensure the working tree is clean.
3. Run `./deploy.sh`. It reruns lint and tests, builds the package, checks that
   main is clean and current and the version is still new, uploads to PyPI,
   and pushes the version tag.

The release is shipped after the upload succeeds. Documentation changes also
publish through the Documentation workflow; check the live site after it runs.

Source archives include the test package and its tracked fixtures. When adding
a fixture, add its exact path to `MANIFEST.in`; avoid recursive fixture globs,
which can package generated caches or developer-local data. CI and `deploy.sh`
extract the sdist and run offline tests before upload. To check locally:

```sh
python -m build --sdist
python scripts/check_sdist.py dist/pyensembl-*.tar.gz
```

From a Git checkout this also verifies all tracked test inputs against the
archive and rejects extra test files. The extracted archive can run its own
`./lint.sh` and `./test.sh`; the full suite still needs the Ensembl datasets
listed in the Tests workflow.

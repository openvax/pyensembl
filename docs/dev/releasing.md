# Publish a release <a id="releasing-pyensembl"></a>

Every PR includes a version bump, including documentation changes. Once lint,
tests and review have passed:

1. Merge the PR into main.
2. Check out main, pull the merged commit, and ensure the working tree is clean.
3. Run `./deploy.sh`. It reruns lint and tests, builds the package, uploads to
   PyPI, and pushes the version tag.

The release is shipped after the upload succeeds. Documentation changes also
publish through the Documentation workflow; check the live site after it runs.

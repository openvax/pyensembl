"""Tests for the PR and deployment release-version guard."""

import pytest

from scripts import check_release_version as release


@pytest.mark.parametrize(
    ("proposed", "base", "tags", "published", "valid"),
    [
        ("2.27.7", "2.27.6", ["2.27.6"], "2.27.6", True),
        ("2.27.5", "2.27.5", ["2.27.4"], "2.27.4", False),
        ("2.27.5", "2.27.6", ["2.27.4"], "2.27.4", False),
        ("2.27.5", "2.27.4", ["2.27.5"], "2.27.4", False),
        ("2.27.5", "2.27.4", ["2.27.4"], "2.27.5", False),
        ("2.28.0", "2.27.9", ["2.27.9"], "2.27.9", True),
    ],
)
def test_version_must_exceed_every_existing_source(
    proposed, base, tags, published, valid
):
    if valid:
        release.check_version(proposed, base=base, tags=tags, published=published)
    else:
        with pytest.raises(ValueError, match="must be newer"):
            release.check_version(proposed, base=base, tags=tags, published=published)


def test_invalid_version_file_is_rejected():
    with pytest.raises(ValueError, match="Missing stable"):
        release.version_from_text('__version__ = "2.27.7rc1"')


def test_pr_check_uses_current_main_and_tags(monkeypatch, capsys):
    calls = []

    def fake_git(*args):
        calls.append(args)
        if args[0] == "show":
            return '__version__ = "2.27.7"'
        if args[0] == "tag":
            return "v2.27.6\nv2.27.7"
        raise AssertionError(args)

    monkeypatch.setattr(release, "git", fake_git)
    monkeypatch.setattr(release, "published_version", lambda: "2.27.6")
    assert release.main(["--base-ref", "origin/main"]) == 1
    assert ("show", "origin/main:pyensembl/version.py") in calls
    assert "current main" in capsys.readouterr().err


def test_pr_check_fails_closed_when_pypi_unavailable(monkeypatch, capsys):
    monkeypatch.setattr(
        release, "git", lambda *args: '__version__ = "2.27.6"'
        if args[0] == "show" else "v2.27.6"
    )

    def unavailable():
        raise OSError("PyPI unavailable")

    monkeypatch.setattr(release, "published_version", unavailable)
    assert release.main(["--base-ref", "origin/main"]) == 1
    assert "PyPI unavailable" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("branch", "status", "head", "origin", "message"),
    [
        ("feature", "", "a", "a", "main branch"),
        ("main", " M file", "a", "a", "clean working tree"),
        ("main", "", "a", "b", "current origin/main"),
    ],
)
def test_deploy_requires_clean_current_main(
    monkeypatch, capsys, branch, status, head, origin, message
):
    values = {
        ("branch", "--show-current"): branch,
        ("status", "--porcelain"): status,
        ("rev-parse", "HEAD"): head,
        ("rev-parse", "origin/main"): origin,
    }
    monkeypatch.setattr(release, "git", lambda *args: values[args])
    monkeypatch.setattr(release, "published_version", lambda: "2.27.6")
    assert release.main(["--deploy"]) == 1
    assert message in capsys.readouterr().err

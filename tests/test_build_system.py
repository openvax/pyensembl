"""
Tests for build system configuration and pyproject.toml.

This module verifies that the package build system is properly configured
with pyproject.toml and that the build process works correctly.
"""

import subprocess
import sys
import os
from pathlib import Path


def test_pyproject_toml_exists():
    """Test that pyproject.toml exists in the project root."""
    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"
    assert pyproject_path.exists(), "pyproject.toml not found in project root"


def test_pyproject_toml_is_valid():
    """Test that pyproject.toml is valid TOML syntax."""
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib  # Python < 3.11

    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    with open(pyproject_path, "rb") as f:
        config = tomllib.load(f)

    # Verify required sections exist
    assert "build-system" in config, "build-system section missing"
    assert "project" in config, "project section missing"


def test_build_system_backend():
    """Test that build-system uses setuptools backend."""
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib

    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    with open(pyproject_path, "rb") as f:
        config = tomllib.load(f)

    assert config["build-system"]["build-backend"] == "setuptools.build_meta"


def test_dependencies_correct():
    """Test that runtime dependencies match specification."""
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib

    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    with open(pyproject_path, "rb") as f:
        config = tomllib.load(f)

    expected_deps = {
        "typechecks>=0.0.2,<1.0.0",
        "datacache>=1.4.0,<2.0.0",
        "memoized-property>=1.0.2",
        "tinytimer>=0.0.0,<1.0.0",
        "gtfparse>=2.6.0,<3.0.0",
        "serializable>=0.2.1,<1.0.0",
        "numpy>=2.0.0,<3.0.0",
    }

    actual_deps = set(config["project"]["dependencies"])

    assert actual_deps == expected_deps, (
        f"Dependencies mismatch.\n"
        f"Expected: {expected_deps}\n"
        f"Actual: {actual_deps}\n"
        f"Missing: {expected_deps - actual_deps}\n"
        f"Extra: {actual_deps - expected_deps}"
    )


def test_no_pylint_in_runtime_deps():
    """
    Test that pylint is not in runtime dependencies.

    As per specification, pylint should only be in dev dependencies if at all,
    not in the main runtime dependencies list.
    """
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib

    project_root = Path(__file__).parent.parent
    pyproject_path = project_root / "pyproject.toml"

    with open(pyproject_path, "rb") as f:
        config = tomllib.load(f)

    runtime_deps = config["project"]["dependencies"]

    # Check that no dependency starts with "pylint"
    pylint_deps = [dep for dep in runtime_deps if dep.lower().startswith("pylint")]
    assert not pylint_deps, f"pylint found in runtime dependencies: {pylint_deps}"

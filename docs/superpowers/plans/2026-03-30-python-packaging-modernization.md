# Python Packaging Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize pyensembl Python packaging from deprecated pkg_resources and setup.py to modern pyproject.toml with importlib.resources.

**Architecture:** Replace deprecated pkg_resources resource loading with Python 3.9+ native importlib.resources, consolidate all package metadata and dependencies into pyproject.toml, and update CI for modern installation patterns.

**Tech Stack:** Python 3.9+ importlib.resources, setuptools build backend, pyproject.toml configuration

---

### Task 1: Create pyproject.toml configuration

**Files:**
- Create: `pyproject.toml`

- [ ] **Step 1: Write failing build test**

```bash
# Test that we can build the package
python -m build --version || echo "build not installed, will install later"
rm -rf dist/ build/ *.egg-info/
python -m pip install build
python -m build
```

Expected: FAIL - no pyproject.toml found

- [ ] **Step 2: Create basic pyproject.toml with build system**

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"
```

- [ ] **Step 3: Add project metadata from setup.py**

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "pyensembl"
version = "2.4.0"
description = "Python interface to Ensembl reference genome metadata"
authors = [
    {name = "Alex Rubinsteyn", email = "alex.rubinsteyn@unc.edu"}
]
license = {text = "Apache-2.0"}
readme = "README.md"
requires-python = ">=3.9"
classifiers = [
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
```

- [ ] **Step 4: Add runtime dependencies from requirements.txt**

```toml
dependencies = [
    "typechecks>=0.0.2,<1.0.0",
    "datacache>=1.4.0,<2.0.0",
    "memoized-property>=1.0.2",
    "tinytimer>=0.0.0,<1.0.0",
    "gtfparse>=2.5.0,<3.0.0",
    "serializable>=0.2.1,<1.0.0",
    "numpy<2",
]
```

- [ ] **Step 5: Add entry point and optional dependencies**

```toml
[project.scripts]
pyensembl = "pyensembl.shell:run"

[project.optional-dependencies]
dev = [
    "flake8",
    "pytest",
    "pytest-cov",
    "ruff",
    "coveralls",
    "build",
]
```

- [ ] **Step 6: Add setuptools package configuration**

```toml
[tool.setuptools]
packages = ["pyensembl"]

[tool.setuptools.package-data]
pyensembl = ["logging.conf"]
```

- [ ] **Step 7: Test build with new pyproject.toml**

```bash
rm -rf dist/ build/ *.egg-info/
python -m build
ls dist/
```

Expected: PASS - creates pyensembl-2.4.0.tar.gz and pyensembl-2.4.0-py3-none-any.whl

- [ ] **Step 8: Commit pyproject.toml**

```bash
git add pyproject.toml
git commit -m "feat: add pyproject.toml configuration

Replaces setup.py with modern Python packaging configuration.
Includes all metadata, dependencies, and build system setup.
"
```

### Task 2: Update resource loading to use importlib.resources

**Files:**
- Modify: `pyensembl/shell.py:43,52`

- [ ] **Step 1: Write test to verify current resource loading works**

```bash
# Test current CLI with logging
python -c "
import sys
sys.path.insert(0, '.')
from pyensembl.shell import run
print('Current resource loading works')
"
```

Expected: PASS - should work with pkg_resources

- [ ] **Step 2: Update import statement**

Replace line 43 in `pyensembl/shell.py`:
```python
# OLD:
import pkg_resources

# NEW:
from importlib import resources
```

- [ ] **Step 3: Update resource loading call**

Replace line 52 in `pyensembl/shell.py`:
```python
# OLD:
logging.config.fileConfig(pkg_resources.resource_filename(__name__, "logging.conf"))

# NEW:
logging.config.fileConfig(resources.files(__name__) / "logging.conf")
```

- [ ] **Step 4: Test updated resource loading**

```bash
# Test new CLI with logging still works
python -c "
import sys
sys.path.insert(0, '.')
from pyensembl.shell import run
print('New resource loading works')
"
```

Expected: PASS - should work with importlib.resources

- [ ] **Step 5: Verify no pkg_resources imports remain**

```bash
grep -r "pkg_resources" pyensembl/
```

Expected: No matches found

- [ ] **Step 6: Run existing tests to ensure nothing broke**

```bash
# Quick smoke test of imports
python -c "import pyensembl; print('Import successful')"
```

Expected: PASS

- [ ] **Step 7: Commit resource loading update**

```bash
git add pyensembl/shell.py
git commit -m "feat: replace pkg_resources with importlib.resources

Modernizes resource loading to use Python 3.9+ standard library
instead of deprecated pkg_resources.
"
```

### Task 3: Update CI workflow for modern installation

**Files:**
- Modify: `.github/workflows/tests.yml:29-34`

- [ ] **Step 1: Update Python version matrix**

In `.github/workflows/tests.yml`, update the matrix:
```yaml
# OLD:
python-version: ["3.9", "3.10", "3.11"]

# NEW:
python-version: ["3.9", "3.10", "3.11", "3.12"]
```

- [ ] **Step 2: Update installation commands**

Replace the installation section (around lines 29-34):
```yaml
# OLD:
- name: Install dependencies
  run: |
    python -m pip install --upgrade pip
    python -m pip install flake8 pytest pytest-cov ruff coveralls
    pip install -r requirements.txt
    pip install .

# NEW:
- name: Install dependencies
  run: |
    python -m pip install --upgrade pip
    pip install .[dev]
```

- [ ] **Step 3: Test CI changes locally**

```bash
# Simulate CI installation process
python -m venv test_env
source test_env/bin/activate
pip install --upgrade pip
pip install .[dev]
pyensembl --version
deactivate
rm -rf test_env
```

Expected: PASS - installs successfully and CLI works

- [ ] **Step 4: Commit CI workflow updates**

```bash
git add .github/workflows/tests.yml
git commit -m "feat: modernize CI installation process

Updates GitHub Actions workflow to use pyproject.toml-based
installation with optional development dependencies.
Adds Python 3.12 to test matrix.
"
```

### Task 4: Remove legacy configuration files

**Files:**
- Delete: `setup.py`, `requirements.txt`

- [ ] **Step 1: Verify build still works without setup.py**

```bash
rm -rf dist/ build/ *.egg-info/
python -m build
```

Expected: PASS - builds successfully using pyproject.toml

- [ ] **Step 2: Test installation from wheel without requirements.txt**

```bash
pip install dist/pyensembl-2.4.0-py3-none-any.whl --force-reinstall
pyensembl --version
```

Expected: PASS - installs and runs correctly

- [ ] **Step 3: Remove setup.py**

```bash
rm setup.py
```

- [ ] **Step 4: Remove requirements.txt**

```bash
rm requirements.txt
```

- [ ] **Step 5: Update MANIFEST.in if it references removed files**

Check if MANIFEST.in references the removed files:
```bash
grep -E "(setup\.py|requirements\.txt)" MANIFEST.in || echo "No references found"
```

If references found, remove them from MANIFEST.in.

- [ ] **Step 6: Test build after cleanup**

```bash
rm -rf dist/ build/ *.egg-info/
python -m build
ls dist/
```

Expected: PASS - still builds correctly

- [ ] **Step 7: Commit removal of legacy files**

```bash
git rm setup.py requirements.txt
git commit -m "feat: remove legacy setup.py and requirements.txt

All configuration now consolidated in pyproject.toml.
Modern Python packaging no longer needs these files.
"
```

### Task 5: Integration testing and verification

**Files:**
- Test: All components working together

- [ ] **Step 1: Clean build and install test**

```bash
# Full clean build and install
rm -rf dist/ build/ *.egg-info/
python -m build
pip install dist/pyensembl-2.4.0-py3-none-any.whl --force-reinstall
```

Expected: PASS - clean build and install

- [ ] **Step 2: Test CLI functionality with logging**

```bash
# Test that logging configuration loads properly
pyensembl --help 2>&1 | head -5
```

Expected: PASS - help output appears without import errors

- [ ] **Step 3: Test development installation**

```bash
# Test editable development install
pip install -e .[dev]
pyensembl --version
```

Expected: PASS - development install works

- [ ] **Step 4: Run linting to ensure code quality**

```bash
# Test that dev dependencies work
ruff check pyensembl/ || echo "Ruff found issues (expected for old code)"
flake8 --version
pytest --version
```

Expected: PASS - dev tools are available

- [ ] **Step 5: Test package imports**

```bash
python -c "
import pyensembl
from pyensembl import Genome
from pyensembl.shell import run
print('All imports successful')
print(f'Version: {pyensembl.version.__version__}')
"
```

Expected: PASS - all imports work correctly

- [ ] **Step 6: Verify no deprecated warnings**

```bash
python -W error::DeprecationWarning -c "
import pyensembl.shell
print('No deprecation warnings')
" 2>&1 || echo "Found deprecation warnings - needs investigation"
```

Expected: PASS - no deprecation warnings from pkg_resources

- [ ] **Step 7: Final commit with verification**

```bash
git add .
git commit -m "feat: complete Python packaging modernization

Successfully migrated from deprecated pkg_resources and setup.py
to modern pyproject.toml configuration with importlib.resources.

- Replaced pkg_resources with importlib.resources
- Consolidated all configuration in pyproject.toml
- Updated CI for modern installation patterns
- Removed legacy setup.py and requirements.txt
- Verified CLI and package functionality

Tested on Python 3.9+ with no breaking changes to public API.
"
```
# Python Packaging Modernization Design

**Date:** 2026-03-30
**Project:** pyensembl
**Type:** Infrastructure Modernization

## Overview

Modernize pyensembl's Python packaging stack to work with current Python tooling by eliminating deprecated `pkg_resources` dependency and migrating from `setup.py` to `pyproject.toml` configuration.

## Problem Statement

The current packaging setup has several issues:
- Uses deprecated `pkg_resources.resource_filename()` in `pyensembl/shell.py:52`
- Relies on legacy `setup.py` configuration (has TODO to migrate to pyproject.toml)
- Dependencies split between `setup.py` and `requirements.txt`
- Missing modern development dependency management
- Build process uses outdated patterns

## Architecture

### Core File Changes

1. **Replace `pkg_resources` usage** in `pyensembl/shell.py` with `importlib.resources`
2. **Create `pyproject.toml`** as single source of truth for project configuration
3. **Remove `setup.py`** and `requirements.txt` to eliminate confusion
4. **Update CI workflows** to use modern installation patterns

### Package Structure

**No changes to package layout** - maintains current structure:
- `pyensembl/` package directory unchanged
- All imports remain identical
- Entry points preserved exactly
- No migration to `src/` layout

### File Operations Summary

- **Add:** `pyproject.toml`
- **Modify:** `pyensembl/shell.py` (1 line change), `.github/workflows/tests.yml`
- **Delete:** `setup.py`, `requirements.txt`

## Dependencies & Build System

### Build Backend
- **Backend:** `setuptools` (familiar, stable, handles current package structure)
- **Frontend:** Modern `python -m build` instead of direct setuptools calls
- **Requirements:** `["setuptools>=61.0", "wheel"]`

### Python Version Support
- **Minimum:** Python 3.9 (matches current CI matrix)
- **Rationale:** Python 3.8 EOL, scientific ecosystem moved to 3.9+, enables native `importlib.resources`
- **Testing:** Python 3.9, 3.10, 3.11, 3.12

### Dependency Migration

**Current runtime dependencies** (from `requirements.txt`):
```
typechecks>=0.0.2,<1.0.0
datacache>=1.4.0,<2.0.0
memoized-property>=1.0.2
tinytimer>=0.0.0,<1.0.0
gtfparse>=2.5.0,<3.0.0
serializable>=0.2.1,<1.0.0
numpy<2
```

**Development dependencies** (move to optional-dependencies):
- `dev` group: ruff, flake8, pytest, pytest-cov, coveralls
- Remove `pylint>=2.17.2,<3.0.0` from runtime deps (development-only)

## Code Changes

### Resource Access Modernization

**File:** `pyensembl/shell.py`
**Lines:** 43, 52

```python
# BEFORE:
import pkg_resources
logging.config.fileConfig(pkg_resources.resource_filename(__name__, "logging.conf"))

# AFTER:
from importlib import resources
logging.config.fileConfig(resources.files(__name__) / "logging.conf")
```

**Benefits:**
- Uses Python 3.9+ standard library
- No deprecated warnings
- Path-like interface more intuitive
- Better type hinting support

## CI/CD Updates

### GitHub Actions Changes

**File:** `.github/workflows/tests.yml`

**Install process:**
```yaml
# BEFORE:
pip install -r requirements.txt
pip install .

# AFTER:
pip install .[dev]
```

**Python matrix:**
- Add Python 3.12 to existing 3.9, 3.10, 3.11
- Keep same test process and Ensembl data installation

### Development Workflow

**Local development:**
- `pip install -e .[dev]` - editable install with dev tools
- `python -m build` - create wheel/sdist for verification
- Single command replaces requirements file management

## Error Handling & Testing

### Backwards Compatibility
- **Package imports:** No changes required in user code
- **CLI interface:** `pyensembl` command unchanged
- **Installation:** Standard `pip install pyensembl` unchanged

### Testing Strategy
- Verify resource loading works in tests
- Build verification in CI with `python -m build`
- Test installation from wheel/sdist
- Ensure dev dependencies properly isolated

### Migration Safety
- No breaking changes to public API
- Package structure preserved
- Entry points identical
- Import paths unchanged

## Implementation Order

1. **Create `pyproject.toml`** with all metadata and dependencies
2. **Update `pyensembl/shell.py`** to use `importlib.resources`
3. **Update CI workflows** for modern install patterns
4. **Remove legacy files** (`setup.py`, `requirements.txt`)
5. **Test build and installation** process
6. **Verify CLI functionality** with resource loading

## Success Criteria

- [ ] No `pkg_resources` imports anywhere in codebase
- [ ] Single `pyproject.toml` configuration file
- [ ] CI passes on Python 3.9-3.12
- [ ] `pyensembl` CLI loads logging.conf correctly
- [ ] Local development with `pip install -e .[dev]` works
- [ ] Package builds successfully with `python -m build`
- [ ] No breaking changes to user-facing API

## Risks & Mitigations

**Risk:** Resource loading fails in some environments
**Mitigation:** Test thoroughly, `importlib.resources` is well-established in 3.9+

**Risk:** Build process changes break CI
**Mitigation:** Update CI incrementally, test each step

**Risk:** Missing development dependencies
**Mitigation:** Comprehensive dev dependency list from current tools

This design provides a clean migration path to modern Python packaging standards while maintaining full backwards compatibility for end users.
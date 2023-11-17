# Compatibility for us old-timers.

# Note: This makefile include remake-style target comments.
# These comments before the targets start with #:
# remake --tasks to shows the targets and the comments

PHONY=check clean dist distclean test
GIT2CL ?= git2cl
PYTHON ?= python
PYTHON3 ?= python3

#: the default target - same as running "check"
all: check

check:
	./lint.sh
	cd test && pytest

#: Clean up temporary files
clean:
	$(PYTHON) ./setup.py $@

#: Create source (tarball) and binary (egg) distribution
dist:
	$(PYTHON) ./setup.py sdist bdist

#: Create source tarball
sdist:
	$(PYTHON) ./setup.py sdist

upload: sdist
	 twine upload dist/*

#: Create binary egg distribution
bdist_egg:
	$(PYTHON) ./setup.py bdist_egg

# It is too much work to figure out how to add a new command to distutils
# to do the following. I'm sure distutils will someday get there.
DISTCLEAN_FILES = build dist *.egg-info *.pyc *.so py*.py

#: Remove ALL dervied files
distclean: clean
	-rm -fr $(DISTCLEAN_FILES) || true

#: Install package locally
install:
	$(PYTHON) ./setup.py install

#: Same as 'check' target
test: check

rmChangeLog:
	rm ChangeLog || true

#: Create a ChangeLog from git via git log and git2cl
ChangeLog: rmChangeLog
	git log --pretty --numstat --summary | $(GIT2CL) >$@

.PHONY: $(PHONY)

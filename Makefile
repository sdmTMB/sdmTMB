PACKAGE=sdmTMB
VERSION := $(shell sed -n '/^Version: /s///p' DESCRIPTION)
DATE := $(shell sed -n '/^Date: /s///p' DESCRIPTION)
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

# Allow e.g. "make R=R-devel install"
R=R

all:
	make doc-update
	make build-package
	make install

doc-update:
	echo "roxygen2::roxygenize(\".\")" | $(R) --no-echo

build-package:
	$(R) CMD build --no-build-vignettes --no-manual .

install:
	$(R) CMD INSTALL --preclean --no-multiarch --with-keep.source .

# devtools::load_all() compiles src/ at -O0 by default, which makes the TMB
# backend several times slower; build an optimized DLL it will reuse instead
compile-opt:
	echo "pkgbuild::clean_dll(); pkgbuild::compile_dll(debug = FALSE)" | $(R) --no-echo

# Quick tests: skip compile-opt (reuses the existing DLL) and use a terse
# reporter. Override with e.g. "make TEST_REPORTER=progress test-quick"
TEST_REPORTER=check
TEST_NCPUS=ncpus <- parallel::detectCores(); options(Ncpus = if (is.na(ncpus)) 1L else max(1L, as.integer(round(ncpus / 2 - 1))))

test:
	echo "$(TEST_NCPUS); devtools::test(reporter = '$(TEST_REPORTER)')" | $(R) --no-echo

test-rtmb:
	echo "$(TEST_NCPUS); devtools::test(reporter = '$(TEST_REPORTER)')" | SDMTMB_TEST_BACKEND=rtmb $(R) --no-echo

test-optional:
	echo "devtools::load_all(); testthat::test_dir('tests/optional')" | $(R) --no-echo

cran-check:
	echo "devtools::check(\".\")" | $(R) --no-echo

# Reference-fit regression suite (see reference-fits/README.md)
reference-check:
	$(R) --no-echo -f reference-fits/run.R --args check --backend both --cores $(TEST_NCPUS)

reference-record:
	$(R) --no-echo -f reference-fits/run.R --args record --cores $(TEST_NCPUS)

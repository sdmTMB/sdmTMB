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
	echo "roxygen2::roxygenize(\".\")" | $(R) --slave

build-package:
	$(R) CMD build --no-build-vignettes --no-manual .

install:
	$(R) CMD INSTALL --preclean --no-multiarch --with-keep.source .

# devtools::load_all() compiles src/ at -O0 by default, which makes the TMB
# backend several times slower; build an optimized DLL it will reuse instead
compile-opt:
	echo "pkgbuild::clean_dll(); pkgbuild::compile_dll(debug = FALSE)" | $(R) --slave

test: compile-opt
	echo "ncpus <- parallel::detectCores(); options(Ncpus = if (is.na(ncpus)) 1L else max(1L, as.integer(round(ncpus / 2)))); devtools::test()" | $(R) --slave

test-rtmb: compile-opt
	echo "ncpus <- parallel::detectCores(); options(Ncpus = if (is.na(ncpus)) 1L else max(1L, as.integer(round(ncpus / 2)))); devtools::test()" | SDMTMB_TEST_BACKEND=rtmb $(R) --slave

# Quick tests: skip compile-opt (reuses the existing DLL) and use a terse
# reporter. Override with e.g. "make TEST_REPORTER=progress test-quick"
TEST_REPORTER=check
TEST_NCPUS=ncpus <- parallel::detectCores(); options(Ncpus = if (is.na(ncpus)) 1L else max(1L, as.integer(round(ncpus / 2 - 1))))


test-quick:
	echo "$(TEST_NCPUS); devtools::test(reporter = '$(TEST_REPORTER)')" | $(R) --slave

test-rtmb-quick:
	echo "$(TEST_NCPUS); devtools::test(reporter = '$(TEST_REPORTER)')" | SDMTMB_TEST_BACKEND=rtmb $(R) --slave

test-optional: compile-opt
	echo "devtools::load_all(); testthat::test_dir('tests/optional')" | $(R) --slave

cran-check:
	echo "devtools::check(\".\")" | $(R) --slave

# Reference-fit regression suite (see reference-fits/README.md)
REF_CORES ?= 6

reference-check: compile-opt
	$(R) --slave -f reference-fits/run.R --args check --backend both --cores $(REF_CORES)

reference-record: compile-opt
	$(R) --slave -f reference-fits/run.R --args record --cores $(REF_CORES)

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

# Tests with a terse reporter. Override with e.g. "make TEST_REPORTER=progress test"
TEST_REPORTER=check
TEST_NCPUS=ncpus <- parallel::detectCores(); options(Ncpus = if (is.na(ncpus)) 1L else max(1L, as.integer(round(ncpus / 2 - 1))))

test:
	echo "$(TEST_NCPUS); devtools::test(reporter = '$(TEST_REPORTER)')" | $(R) --no-echo

test-optional:
	echo "devtools::load_all(); testthat::test_dir('tests/optional')" | $(R) --no-echo

cran-check:
	echo "devtools::check(\".\")" | $(R) --no-echo

# Knit README.Rmd and losslessly compress its figures (requires optipng)
.PHONY: readme
readme:
	echo "devtools::load_all(quiet = TRUE); rmarkdown::render('README.Rmd', quiet = TRUE)" | $(R) --no-echo
	optipng -o2 -quiet man/figures/README-*.png

# Reference-fit regression suite (see reference-fits/README.md)
REF_CORES=5
reference-check:
	$(R) --no-echo -f reference-fits/run.R --args check --backend rtmb --cores $(REF_CORES)

reference-record:
	$(R) --no-echo -f reference-fits/run.R --args record --cores $(REF_CORES)

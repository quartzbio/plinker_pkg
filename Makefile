PKG=plinker
PARALLEL=TRUE
FILTER=
QUIET=TRUE

LIB=.checks

qbtest:
	Rscript --no-save -e 'suppressPackageStartupMessages(library(qbdev));load_pkg("$(PKG)", quiet = TRUE)'
	LONG_TESTS=FALSE Rscript --no-save -e 'library(qbdev);test_pkg("$(PKG)", parallel=$(PARALLEL), filter="$(FILTER)", quiet = $(QUIET))'

.checks:
	mkdir -p .checks

qbcheck: .checks
	Rscript --no-save -e 'library(qbdev);check_pkg("$(PKG)", "$(LIB)", "$(LIB)", parallel=$(PARALLEL))'

rox:
	Rscript --no-save -e 'devtools::document("$(PKG)")'

test:
	Rscript --no-save -e 'devtools::test("$(PKG)", filter="$(FILTER)")'

check: .checks
	Rscript --no-save -e 'devtools::check("$(PKG)", cleanup = FALSE, check_dir = "$(LIB)")'


build:
	Rscript --no-save -e 'devtools::build("$(PKG)")'

install: rox
	R CMD INSTALL $(PKG)

manual: rox
	rm -f $(PKG).pdf
	R CMD Rd2pdf -o $(PKG).pdf $(PKG)

coverage:
	Rscript --no-save -e 'covr::package_coverage("$(PKG)", quiet = FALSE)'

zero-coverage:
	Rscript -e 'library(covr); zero_coverage(package_coverage("$(PKG)"))'

coverage.html:
	Rscript -e 'library(covr); report(package_coverage("$(PKG)"), file = "$@")'

roxygen/clean:
	@rm -rf $(PKG)/man $(PKG)/NAMESPACE

clean: roxygen/clean 
	@rm -rf .checks *.Rcheck .Rd2*  $(PKG)/src/*.o $(PKG)/src/*.so $(PKG)/src/*.gcda


PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: install

build:
	cd ..;\
	R CMD build $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGSRC)

check: 
	Rscript -e "devtools::check(document = FALSE)"

news:
	sed -e 's/^-/  -/' -e 's/^## *//' -e 's/^# //' <NEWS.md | fmt -80 >NEWS

test:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);options(warn=2);test_dir('tests/tests')"

paralleltest:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);mirtCluster();options(warn=2);test_dir('tests/tests')"

extratest:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);options(warn=2);test_dir('tests/extratests')"

knitdocs:
	sed -i 's/# opts$verbose <- FALSE/opts$verbose <- FALSE/g' R/03-estimation.R
	sed -i 's/dontrun/donttest/g' man/*.Rd
	Rscript -e "devtools::install('mirt');library('knitr',quietly=TRUE);knit_rd('$(PKGNAME)')"
	mkdir html
	mv *.html html/
	rm R.css
	rm -r figure/
	git checkout -- .
	make install

clean:
	$(RM) src/*.o
	$(RM) src/*.so
	$(RM) ../$(PKGNAME)_$(PKGVERS).tar.gz
	$(RM) -r ../$(PKGNAME).Rcheck/



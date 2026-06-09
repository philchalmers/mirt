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
	Rscript -e "devtools::check(document = FALSE, args = '--as-cran')"

news:
	sed -e 's/^-/  -/' -e 's/^## *//' -e 's/^# //' <NEWS.md | fmt -80 >NEWS

test:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);options(warn=2);test_dir('tests/testthat')"

paralleltest:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);mirtCluster();options(warn=2);test_dir('tests/testthat')"

extratest:
	Rscript -e "library('testthat',quietly=TRUE);library('mirt',quietly=TRUE);options(warn=2);test_dir('tests/testthat/extratests')"

# knitdocs:
#	rm -rf html/
#	sed -i 's/# opts$$verbose/opts$$verbose/g' R/03-estimation.R
#	sed -i s/dontrun/donttest/g man/*.Rd
#	make install
#	Rscript -e "library('mirt');mirtCluster(2);library('knitr',quietly=TRUE);knit_rd('$(PKGNAME)')"
#	mkdir html
#	mv *.html html/
#	rm R.css
#	rm -rf figure/
#	git checkout -- .
#	make install
#	make kniterrors

pushdocs:
	mv docs/ ../mirtdocs
	git checkout gh-pages
	cp -r ../mirtdocs/* docs/
	git commit -am "update pkgdown"
	git push
	git checkout main

pkgdownerrors:
	cd docs/reference
	printf "\n#########################\nErrors:\n\n"
	grep -Hrn '>Error:'
	printf "\n#########################\nWarnings:\n\n"
	grep -Hrn '>Warning:'

pkgdown:
	sed -i 's/# opts$$verbose/opts$$verbose/g' R/03-estimation.R
	sed -i s/dontrun/donttest/g man/*.Rd
	make install
	Rscript -e "library('pkgdown',quietly=TRUE);build_site()"
	git checkout -- .
	make install
	make pkgdownerrors

clean:
	$(RM) src/*.o
	$(RM) src/*.so
	$(RM) ../$(PKGNAME)_$(PKGVERS).tar.gz
	$(RM) -r ../$(PKGNAME).Rcheck/



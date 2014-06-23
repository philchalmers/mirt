all: vignettes

vignettes:
	Rscript -e "library('knitr');files=dir();for(file in files[grepl('*.rmd',tolower(files))]) knit2html(file)"
	$(RM) *.md;\
	$(RM) -r figure;\
	mv -f *.html html/

manuscript = Thesis.tex
#includes = grab all lines starting with \input in Thesis.tex
includes = $(shell awk '/^\\input/ {i=index($$0,"{");j=index($$0,"}");print substr($$0,i+1,j-i-1)}' $(manuscript))
subincludes = $(shell awk '/^\\include/ {i=index($$0,"{");j=index($$0,"}");print substr($$0,i+1,j-i-1)".tex"}' $(manuscript))
#figures = grab all lines starting with \includegraphics in all .tex files from above
figures = $(shell awk '/^\\includegraphics/ {i=index($$0,"{");j=index($$0,"}");print substr($$0,i+1,j-i-1)}' $(manuscript) $(includes))
#bibliography = grab line starting with \bibliography in Paper.tex, then add .bib at the end
bibliography = $(addsuffix .bib,$(shell awk '/^\\bibliography{/ {i=index($$0,"{");j=index($$0,"}");print substr($$0,i+1,j-i-1)}' $(manuscript)))
#bibliography = AutoBibliography.bib

tar = $(patsubst %.tex,%.tar.gz,$(manuscript))

#Prevents intermediate files (see http://theory.uwinnipeg.ca/localfiles/infofiles/make/make_94.html#SEC93) from being deleted
.SECONDARY:

all: pdf

#Changes ".tex" to ".pdf"
pdf: $(subst .tex,.pdf,$(manuscript))

ps: $(subst .tex,.ps,$(manuscript))

bib: AutoBibliography.bib

tar: $(tar)

arxiv: $(arxiv)

$(subst .tex,.pdf,$(manuscript)): $(manuscript) $(includes) $(subincludes) $(figures) $(subst .tex,.bbl,$(manuscript))
	pdflatex $<
	pdflatex $<

$(subst .tex,.dvi,$(manuscript)): $(manuscript) $(includes) $(figures) $(subst .tex,.bbl,$(manuscript))
	latex $<
	latex $<

$(patsubst %.tex,%.aux,$(manuscript) $(includes)): $(manuscript) $(includes) $(subincludes)
	latex $<

$(subst .tex,.bbl,$(manuscript)): $(bibliography)
	bibtex $(subst .bbl,,$@)

#Added the first line to generate an .aux where none exists
AutoBibliography.bib: $(manuscript) $(includes) $(subincludes) fetch_bibliography.py
	pdflatex $(manuscript)
	#python fetch_bibliography.py -filename $(manuscript) -addition bib_addition.bib

#	cat $(filter %.aux,$(patsubst %.tex,%.aux,$<)) >$(subst .bib,.aux,$@)
#	ln -sf $< $(subst .bib,.tex,$@)
#	./makebib $(subst .bib,,$@)

$(tar): $(manuscript) $(includes) $(subincludes) $(figures) $(bibliography) Makefile fetch_bibliography.py
	tar czvfh $@ $^

%.ps: %.dvi
	dvips $*

#$(subst .tex,.bbl,$(manuscript))
test:
	@echo $(manuscript) $(includes) $(subincludes) $(figures)

clean:
	rm -rf Thesis.pdf *.aux *.dvi *.bib *.bbl *.blg

################################################################################
# Cheng-Han Yu
# 09/2016
# Description of File:
# Makefile for knitr compiling
################################################################################
all:pdf  # default rule DO NOT EDIT
################################################################################
MAINFILE  := AMS223_HW1
RNWFILES  := AMS223_HW1.Rnw
RFILES    := 
TEXFILES  := 
CACHEDIR  := cache
FIGUREDIR := figure
LATEXMK_FLAGS := 

##### Explicit Dependencies #####
################################################################################
RNWTEX = $(RNWFILES:.Rnw=.tex)
ROUTFILES = $(RFILES:.R=.Rout)
RDATAFILES= $(RFILES:.R=.Rdata)
MAINTEX = $(MAINFILE:=.tex)
MAINPDF = $(MAINFILE:=.pdf)
ALLTEX = $(MAINTEX) $(RNWTEX) $(TEXFILES)

# Dependencies
$(RNWTEX): $(RDATAFILES)
$(MAINTEX): $(RNWTEX) $(TEXFILES)
$(MAINPDF): $(MAINTEX) $(ALLTEX) 

.PHONY: pdf tex clean cleanmore clearcache cleanall
pdf: $(MAINPDF)
tex: $(RDATAFILES) $(ALLTEX) 

$(CACHEDIR):
	mkdir $(CACHEDIR)

$(FIGUREDIR):
	mkdir $(FIGUREDIR)

%.tex:%.Rnw
	Rscript \
	  -e "library(knitr)" \
	  -e "knitr::opts_chunk[['set']](fig.path='$(FIGUREDIR)/$*-')" \
	  -e "knitr::opts_chunk[['set']](cache.path='$(CACHEDIR)/$*-')" \
	  -e "knitr::knit('$<','$@')"


%.R:%.Rnw
	Rscript -e "Sweave('$^', driver=Rtangle())"

%.Rout:%.R
	R CMD BATCH "$^" "$@"

%.pdf: %.tex 
	latexmk -pdf $<
    
clean:
	#rm -f *.{fls,aux,fdb_latexmk,log}
	-latexmk -c -quiet $(MAINFILE).tex

cleanfigure:
	-rm -f $(FIGUREDIR)/*.{aux,out,log,tikz}

cleantex:
	-rm -f $(MAINTEX) $(RNWTEX)

clearcache:
	-rm -rf cache

cleanall: clean cleanfigure cleantex clearcache
	-rm *tikzDictionary
	-rm $(MAINPDF)
	-rm -rf $(FIGUREDIR)





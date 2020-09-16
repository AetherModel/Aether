#
# This file was modified from the mpich-3.3.2/doc/design/Makefile
# Rev0 Date: Sept. 16, 2020
# Rev0 Author: Aaron Ridley
# Modifications:
#

ALL: all-redirect
SHELL = /bin/sh

srcdir = .
.SUFFIXES: .pdf .dvi .tex

VPATH=.:${srcdir}

LATEX      = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; latex
BIBTEX     = BIBINPUTS=".:$(srcdir):" ; export BIBINPUTS ; bibtex
DVIPDFM    = TEXINPUTS=".:$(srcdir):" ; export TEXINPUTS ; dvipdfm
LATEX2HTML = latex2html

.tex.dvi:
	$(LATEX) $< </dev/null >/dev/null 2>&1
	$(LATEX) $< </dev/null >/dev/null 2>&1
	# $(BIBTEX) $* </dev/null
	$(LATEX) $<

.dvi.pdf:
	$(DVIPDFM) $*

all-redirect: standards.pdf

clean:
	rm -f *.toc *.out *.dvi *.aux *.blg *.bbl *.log *.idx
	rm -f standards.pdf

mandoc:

htmldoc:

latexdoc: ALL

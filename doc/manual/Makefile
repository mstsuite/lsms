TEX_FILES = lsms_manual.tex \
            constrained_moment.tex instalation.tex units.tex

all: pdf

pdf: $(TEX_FILES)
	pdflatex lsms_manual.tex
	pdflatex lsms_manual.tex

clean:
	rm -f lsms_manual.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.toc

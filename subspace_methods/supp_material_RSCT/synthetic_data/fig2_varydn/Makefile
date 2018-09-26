# Makefile to generate fig_varydn.pdf
# Reinhard Heckel, 2013


all: fig_varydn.pdf

fig_varydn.pdf: CE_varydn.dat
	pdflatex -halt-on-error -interaction=batchmode fig_varydn.tex

CE_varydn.dat:
	matlab -nojvm < varydn.m

clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux


# Makefile to generate fig_varydn.pdf
# Reinhard Heckel, 2013


all: fig_incomplete_data.pdf

fig_incomplete_data.pdf: erasure_20.dat
	pdflatex -halt-on-error -interaction=batchmode fig_incomplete_data.tex 

erasure_20.dat:
	matlab -nojvm < incomplete_data.m

clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux


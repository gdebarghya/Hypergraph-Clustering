# Makefile to generate fig_intersection.pdf
# Reinhard Heckel, 2013


all: fig_intersection.pdf

fig_intersection.pdf: CE_SSintersect_TSC.dat
	pdflatex fig_intersection.tex

CE_SSintersect_TSC.dat:
	matlab -nojvm < intersection_subspaces.m

clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux


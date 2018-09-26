# Makefile to generate fig_digits.pdf and fig_sigvals.pdf
# Reinhard Heckel, 2013

all: fig_digits.pdf fig_singvals.pdf

fig_digits.pdf: TSC.dat
	pdflatex fig_digits.tex
TSC.dat:
	matlab -nojvm < compare_algorithms.m
##
fig_singvals.pdf: svals.dat
	pdflatex fig_singvals.tex
svals.dat:
	matlab -nojvm < lie_digits_in_subspaces.m

##
clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux


## Generate all the tikz images

TEX_FILES = chi_less_theta.tex chi_more_theta.tex chi_equal_theta.tex
PDF_FILES = chi_less_theta.pdf chi_more_theta.pdf chi_equal_theta.pdf

all : tex clean mv
.PHONY: tex clean mv all

tex:
	$(foreach var,$(TEX_FILES),pdflatex $(var);)

clean:
	rm -f *.log
	rm -f *.aux
	rm -f *.synctex*



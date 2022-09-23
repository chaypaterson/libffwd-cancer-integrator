tmpfile = CPaterson_Documentation
make : documentation.tex
	cp documentation.tex $(tmpfile).tex
	pdflatex $(tmpfile).tex
	bibtex $(tmpfile).aux
	pdflatex $(tmpfile).tex
	pdflatex $(tmpfile).tex
	rm $(tmpfile).tex

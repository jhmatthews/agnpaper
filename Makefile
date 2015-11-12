FILE=draft
FILE2=paper
FILE3=paper_mnras
DIFF_FILE=diffs_mnras

paper:
	python mnras.py paper.tex paper_mnras.tex
	pdflatex ${FILE2}
	bibtex ${FILE2}
	pdflatex ${FILE2}
	pdflatex ${FILE2}
	#dvips -o ${FILE2}.ps ${FILE2}
	#ps2pdf ${FILE2}.ps ${FILE2}.pdf 

	pdflatex ${FILE3}
	bibtex ${FILE3}
	pdflatex ${FILE3}
	pdflatex ${FILE3}
	#
	#cp ${FILE}.pdf ~/Dropbox/Documents/
	#dvips -o ${FILE3}.ps ${FILE3}
	#ps2pdf ${FILE3}.ps ${FILE3}.pdf 	
	open -a preview ${FILE3}.pdf

paperdvi:
	python mnras.py paper.tex paper_mnras.tex
	latex ${FILE2}
	bibtex ${FILE2}
	latex ${FILE2}
	latex ${FILE2}
	dvips -o ${FILE2}.ps ${FILE2}
	ps2pdf ${FILE2}.ps ${FILE2}.pdf 

	latex ${FILE3}
	bibtex ${FILE3}
	latex ${FILE3}
	latex ${FILE3}
	#
	#cp ${FILE}.pdf ~/Dropbox/Documents/
	dvips -o ${FILE3}.ps ${FILE3}
	ps2pdf ${FILE3}.ps ${FILE3}.pdf 	
	open -a preview ${FILE3}.pdf


draft: 
	latex draft
	bibtex draft
	latex draft
	latex draft
	#
	#cp ${FILE}.pdf ~/Dropbox/Documents/
	dvips -o draft.ps draft
	ps2pdf draft.ps draft.pdf 	
	open -a preview draft.pdf
#	gv ${FILE}.ps &
	@echo "WORDCOUNT"
	#wc report.out

ldiff:
	latexdiff --flatten compare.tex paper.tex > diffs.tex
	python mnras.py diffs.tex ${DIFF_FILE}.tex

diffs:
	#latexdiff --flatten compare.tex paper.tex > diffs.tex
	#python mnras.py diffs.tex ${DIFF_FILE}.tex

	pdflatex ${DIFF_FILE}
	bibtex ${DIFF_FILE}
	pdflatex ${DIFF_FILE}
	pdflatex ${DIFF_FILE}

	open -a preview ${DIFF_FILE}.pdf 


	
clean:	
	/bin/rm -f *.aux 



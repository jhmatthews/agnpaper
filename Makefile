FILE=draft

all: 
	pdflatex ${FILE}
	bibtex ${FILE}
	pdflatex ${FILE}
	pdflatex ${FILE}
	open -a preview ${FILE}.pdf
	cp ${FILE}.pdf ~/Dropbox/Documents/
#	dvips -o ${FILE}.ps ${FILE}
#	ps2pdf ${FILE}.ps ${FILE}.pdf 	
#	gv ${FILE}.ps &
	@echo "WORDCOUNT"
	wc report.out
	
clean:	
	/bin/rm -f *.aux 



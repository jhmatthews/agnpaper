FILE=paper

paper:
	latex ${FILE}
	# bibtex ${FILE}
	# latex ${FILE}
	# latex ${FILE}
	dvips -o ${FILE}.ps ${FILE}
	ps2pdf ${FILE}.ps ${FILE}.pdf 	
	open -a preview ${FILE}.pdf

	
clean:	
	/bin/rm -f *.aux 
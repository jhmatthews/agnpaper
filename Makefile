FILE=draft
FILE2=paper

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

paper:
	latex ${FILE2}
	bibtex ${FILE2}
	latex ${FILE2}
	latex ${FILE2}
	#
	#cp ${FILE}.pdf ~/Dropbox/Documents/
	dvips -o ${FILE2}.ps ${FILE2}
	ps2pdf ${FILE2}.ps ${FILE2}.pdf 	
	open -a preview ${FILE2}.pdf
#	gv ${FILE}.ps &
	@echo "WORDCOUNT"
	#wc report.out



	
clean:	
	/bin/rm -f *.aux 



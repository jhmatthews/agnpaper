

import sys
import numpy as np
filename = sys.argv[1]
outfilename = sys.argv[2]


texfile =  open(filename, "r")
outfile = open(outfilename, "w")


outfile.write("\documentclass[useAMS,usenatbib]{mn2e_x}\n")

fig_env = False

for line in texfile:
	if "begin{abstract}" in line:
		outfile.write("\maketitle\n")
		outfile.write(line)

	elif "\\author{" in line:
		outfile.write("\\author[Matthews et al.]{\n")

	elif "indent" in line and "rule" in line:
		if "%" not in line:
			outfile.write("\\noindent\\rule{8cm}{0.4pt}\n")

	elif "fullpage" in line and "begin{figure}" in line:
		if fig_env == False:
			outfile.write("\\begin{figure*}\n")
			fig_env = True 

	elif fig_env and "end{figure}" in line:
		outfile.write("\end{figure*}")
		outfile.write("\n")
		fig_env = False 

	elif "documentclass" not in line and "maketitle" not in line:
		outfile.write(line)



outfile.close()

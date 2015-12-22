#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *

eddmod =[0.2]
mbhmod = [1e9]


mbh, l = np.loadtxt("balqsos_mbh_ledd.dat", unpack=True, usecols=(2,3))

names1, names2 = np.loadtxt("balqsos_mbh_ledd.dat", unpack=True, usecols=(0,1), dtype="str")


names = []
locs = []
for i in range(len(names1)):

	names.append(names1[i] + ' ' + names2[i])

	print names[i]

	#if names[i] == 'PG 1001+054' or names[i] == 'Mkn 231':
	if names[i] == 'Mkn 231':
		locs.append(i)
#names = names1 + names2
#print names

set_pretty()

colors=get_colors()

scatter(10.0**mbh, 10.0**l, marker="x", s=100, c="k")
scatter(mbhmod, eddmod, marker="*", s=100, c="r")


for j in locs:
	print mbh[j], l[j], names[j]
	text(mbh[j], l[j], names[j])

xlabel(r"$M_{BH} (M_\odot)$", fontsize=20)
ylabel(r"$L/L_{edd}$", fontsize=20)

loglog()

savefig("egv.png", dpi=500)





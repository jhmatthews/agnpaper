#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *

def get_linelist(filenames, ionreturn=False):

	lines = np.array([])
	z = np.array([])

	for i in range(len(filenames)):

		x = np.loadtxt(filenames[i], comments="#", usecols=(2,3), unpack=True)

		lines = np.concatenate((lines, x[1]))

		z = np.concatenate((z, x[0]))

	if ionreturn:
		return lines, z
	else:
		return lines

z, ion, wave, f = np.loadtxt("data/atomic78/lines_linked_ver_2.py", comments="#", usecols=(1,2,3,4), unpack=True)

set_pretty()
#select = (z == )
E = 1e-3 * HEV * C / (wave * ANGSTROM)

#zz = np.array([6,7,8,10,26])

zz, ions = np.loadtxt("elems", unpack=True, usecols=(1,2), dtype='string')
#ions = ["C","N","O","Ne","Fe"] 
zz = zz.astype(int)

ions = ions[(zz > 3)]
zz = zz[(zz > 3)]


colors = np.concatenate((get_colors(),get_colors(),get_colors()))

for j in range(len(zz)):
	select = (z == zz[j])
	scatter(E[select], ion[select], s=200*f[select], label=ions[j],c=colors[j])

float_legend()
xlim(0.3,7)
ylim(0,27)
semilogx()
xlabel("Line Energy (eV)")
ylabel("Ion")
savefig("ions.png",dpi=200)
clf()
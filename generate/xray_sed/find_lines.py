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

	'''
	function to read in linelist from a bunch of filenames
	'''

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


# some global plotting params
logx = True
# only select lines above this oscillator threshold
oscillator_threshold = 0.01
xlims = (0.4,3.5)
elim = (0.4,3)
ionlim = (15,27)	# ion limits

# angle we are going to plot at
angle = 75



# read in Python synthetic spectrum
fname = "grid11/nextgen_a05_pre_xray"
s = r.read_spectrum("../specs/"+fname)

nu = s["Freq."]

# converting from fnu to flambda
fl_to_fnu = s["Lambda"] / s["Freq."]

# converting from flux to normalisation
norm = 4.0 * PI * (100.0 * PARSEC) * (100.0 * PARSEC)

#convert to lnu - flambda x fl_to_fnu x norm
lnu = s["A%iP0.50" % angle] * fl_to_fnu * norm




# now get the atomic data
z, ion, wave, f = np.loadtxt("data/atomic78/lines_linked_ver_2.py", comments="#", usecols=(1,2,3,4), unpack=True)

energies = C / (wave * ANGSTROM) * HEV * 1e-3 # energy in keV

# print the strong lines
for i in range(len(f)):
	if f[i] > oscillator_threshold and energies[i] < elim[1] and energies[i] > elim[0]:
		print z[i], ion[i], f[i], energies[i]

# select only strong lines
select2 = (f > oscillator_threshold) * (energies < elim[1]) * (energies > elim[0]) * (z < ionlim[1]) * (z > ionlim[0])

# read in the elements used in the simulation
zz, elem_names = np.loadtxt("elems", unpack=True, usecols=(1,2), dtype='string')

zz = zz.astype(int)
elem_names = elem_names[(zz > ionlim[0]) * (zz < ionlim[1])]
zz = zz[(zz > ionlim[0]) * (zz < ionlim[1])]



figure(figsize=(12,12))

# cosmetics
set_pretty()
big_tick_labels(16)
long_ticks()
colors = np.concatenate((get_colors(),get_colors(),get_colors()))


# first create the first subplot which is our spectrum
subplot(2,1,1)
# put vertical lines where the lines are
vlines(energies[select2], 1e41, 1e45)

# plot the spectrum
plot(nu * HEV * 1e-3, util.smooth(nu*lnu, window_len = 20), c=colors[0], label="Model", linewidth=2)


# set limits
semilogy()
if logx:
	semilogx()
xlim(xlims[0],xlims[1])
ylim(1e41,1e45)

#xlabel(r"$E$ (eV)", fontsize=20)
ylabel(r"$\nu L_\nu$ (erg~s$^{-1}$)", fontsize=20)



subplot(2,1,2)

for j in range(len(zz)):
	select = (z == zz[j]) * select2

	if np.sum(select) > 0:
		vlines(energies[select],ion[select],27)
		scatter(energies[select], ion[select], s=800*f[select], label=elem_names[j],c=colors[j])


text(1.5,7,r"Size~$\propto$~Oscillator Strength", fontsize=16)

# put a legend
float_legend()

# set limits and log
xlim(xlims[0],xlims[1])
ylim(4,ionlim[1])
if logx:
	semilogx()

xlabel("Energy (keV)", fontsize=20)
ylabel("Ionic Stage", fontsize=20)

subplots_adjust(hspace=0, wspace=0)
savefig("ions%ito%i.png" % (ionlim[0],ionlim[1]),dpi=200)
clf()







# f_2kev = 2000.0/HEV
# f_2500 = C/(2500.0 / 1e8)
# vlines([f_2kev, f_2500],0.2e47,1e47,linestyle="-", linewidth=1.2)
# text(f_2kev*1.1,5e46,"$2$~keV")
# text(f_2500*1.1,5e46,"$2500$~\AA")


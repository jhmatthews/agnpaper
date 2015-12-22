#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *
from plot_norm import make_f

set_pretty()
NORM = True
almost_black = '#262626'
#nameb = "nextgen_alpha1_long"
namea = "nextgen_a5"
nameb = "nextgen_opt"
# namec = "asi_thmin70_rmin50_a1_f0p01_mdotw5"
# namec = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"

name2 = "fiducial"
smoda = r.read_spectrum("../specs/"+namea)
smodb = r.read_spectrum("../specs/"+nameb)
# smodc = r.read_spectrum("specs/"+namec)


sdisk = r.read_spectrum("../specs/disk")
sh13 = r.read_spectrum("../specs/h13")


angles = [20,60,72,75,80]

colors = get_colors()

figure(figsize=(8,12))



#run_agnspec

for i in range(len(angles)):

	angle = angles[i]
	subplot(5,1,i+1)
	angstr = "A%iP0.50" % angle
	wwa = smoda["Lambda"]
	wwb = smodb["Lambda"]
	select = (wwa <= 2900)
	select2 = (wwb > 2900)
	# wwc = smodc["Lambda"]


	#Fdisk = sdisk[angstr]
	#Fh13 = sh13[angstr]

	Fmoda = smoda[angstr]
	Fmodb = smodb[angstr]
	#Fmodc = smodc[angstr]



	plot(wwa[select],util.smooth(Fmoda[select]), label=r"Model A, $\alpha = 0.6$", c='k', linewidth=2)
	plot(wwb[select2],util.smooth(Fmodb[select2]), c='k', linewidth=2)
	#plot(wwc,util.smooth(Fmodc), label=r"Model C, $\alpha = 1$", c=colors[2], linewidth=2)
	#plot(sh13["Lambda"],util.smooth(Fh13), label="H13", c=almost_black, linewidth=1)
	#plot(sdisk["Lambda"],util.smooth(Fdisk, window_len=100), label="Disc continuum", c="k", linestyle="--")


	NORM = True
	if NORM:
		f2000 = util.get_flux_at_wavelength(wwa,Fmoda,2000)
	#plot(ww,util.smooth(Fmod/f2000), label="Original Model")

	#read the composites and plot them
	if angles[i]<=70:

		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdss_hst_combined.dat", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="non-BAL composite")

	else:
		#names = ["nonbal","hibal", "lobal","felobal"]
		names = ["nonbal","hibal", "lobal"]
		labels = ["non-BAL composite", "hiBAL composite", "loBAL composite"]

		for j in range(len(names)-1):
			w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.%s.spec" % names[j], unpack=True,usecols=(0,1))
			f2000_comp = util.get_flux_at_wavelength(w,f,2000)
			plot(w,f/f2000_comp*f2000, label=labels[j])

	start, end = gca().get_ylim()
	text(2950,start+(0.5*(end-start)),"$%i^\circ$" % angle, fontsize=16, rotation="vertical")
	xlim(800,6800)
	#ylim(0,6.9)

	semilogy()


	if i<3:
		gca().set_xticklabels([])
	if i == 2:
		#text(600,6,"$F_{\lambda}$ at 100pc (erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)", fontsize=20, rotation="vertical")
		ylabel("$F_{\lambda}$ at 100pc (erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)", fontsize=20)
	if i == 0: float_legend()

subplots_adjust(hspace=0)
xlabel("Wavelength (\AA)", fontsize=20)
savefig("wide.png", dpi=500)
#savefig("../../figures/uvspec.eps", bbox="tight")

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

path = "."

spec_files = [f for f in os.listdir(path) if (f.endswith('.spec') and f!='lines.spec')]


for iname in range(len(spec_files)):



	almost_black = '#262626'
	name1 = "nextgen_alpha1_long"
	name2 = spec_files[iname][:-5]
	specname = name2
	smod = r.read_spectrum(name1)
	sdisk = r.read_spectrum("disk_only")
	sh13 = r.read_spectrum(name2)

	print name2

	angles = [20,40, 75,80]

	figure(figsize=(12,8))

	#run_agnspec

	for i in range(len(angles)):

		angle = angles[i]
		subplot(4,1,i+1)
		angstr = "A%iP0.50" % angle
		ww = smod["Lambda"]


		Fdisk = sdisk[angstr]
		Fmod = smod[angstr]
		Fh13 = sh13[angstr]



		plot(ww,util.smooth(Fmod), label="Full model", c="k", linewidth=2)
		plot(sh13["Lambda"],util.smooth(Fh13), label="for comp", c=almost_black, linewidth=1)
		plot(sdisk["Lambda"],util.smooth(Fdisk, window_len=100), label="Disc continuum", c="k", linestyle="--")



		# if NORM:
		# 	f2000 = util.get_flux_at_wavelength(ww,Fmod,2000)
		# plot(ww,util.smooth(Fmod/f2000), label="Original Model")

		# read the composites and plot them
		# if angles[i]<=70:

		# 	w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.nonbal.spec", unpack=True,usecols=(0,1))
		# 	f2000 = util.get_flux_at_wavelength(w,f,2000)
		# 	plot(w,f/f2000, label="non-BAL composite")

		# else:
		# 	#names = ["nonbal","hibal", "lobal","felobal"]
		# 	names = ["nonbal","hibal", "lobal"]
		# 	labels = ["non-BAL composite", "hiBAL composite", "loBAL composite"]

		# 	for j in range(len(names)-1):
		# 		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.%s.spec" % names[j], unpack=True,usecols=(0,1))
		# 		f2000 = util.get_flux_at_wavelength(w,f,2000)
		# 		plot(w,f/f2000, label=labels[j])

		start, end = gca().get_ylim()
		text(2950,start+(0.5*(end-start)),"$%i^\circ$" % angle, fontsize=16, rotation="vertical")
		xlim(800,2900)
		#ylim(0,6.9)


		if i<3:
			gca().set_xticklabels([])
		if i == 3:
			text(600,6,"$F_{\lambda}$ at 100pc (erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)", fontsize=20, rotation="vertical")

		if i == 0: float_legend()

	subplots_adjust(hspace=0)
	xlabel("Wavelength (\AA)", fontsize=20)
	savefig("uvspec_%s.png" % name2)

	#savefig("../../figures/uvspec.eps", bbox="tight")

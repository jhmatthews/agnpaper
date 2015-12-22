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
namea = "nextgen_uv"
namea = "checkvol_a6"
#namea = "nextgen_a5"
nameb = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"
namec = "asi_thmin70_rmin50_a1_f0p01_mdotw5"
#namec = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"

name2 = "fiducial"
smoda = r.read_spectrum("../specs/"+namea)
#smodb = r.read_spectrum("../specs/"+nameb)
#smodc = r.read_spectrum("../specs/"+namec)


sdisk = r.read_spectrum("../specs/disk")
sh13 = r.read_spectrum("../specs/h13")


angles = [20,60,72,75,78]

colors = get_colors()

figure(figsize=(10,10))

#run_agnspec

for i in range(len(angles)):

	angle = angles[i]
	subplot(5,1,i+1)
	angstr = "A%iP0.50" % angle
	wwa = smoda["Lambda"]
	#wwb = smodb["Lambda"]
	#wwc = smodc["Lambda"]

	NORM_GPC = (1e9*1e9)/(1e4)/1e14


	Fdisk = sdisk[angstr] / NORM_GPC
	Fh13 = sh13[angstr] / NORM_GPC

	Fmoda = smoda[angstr] / NORM_GPC
	#Fmodb = smodb[angstr]
	#Fmodc = smodc[angstr]



	plot(wwa,util.smooth(Fmoda,window_len=5), label=r"Model", c='k', linewidth=2)
	#plot(wwb,util.smooth(Fmodb), label=r"Model B, $\alpha = 0.75$", c=colors[1], linewidth=2)
	#plot(wwc,util.smooth(Fmodc), label=r"Model C, $\alpha = 1$", c=colors[2], linewidth=2)
	#plot(sh13["Lambda"],util.smooth(Fh13,window_len=35), label="H13", c=almost_black, linewidth=1)
	#plot(sdisk["Lambda"],util.smooth(Fdisk, window_len=100), label="Disc continuum", c="k", linestyle="--")

	f2000= util.get_flux_at_wavelength(wwa,Fmoda,2000)



	w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdss_hst_combined.dat", unpack=True,usecols=(0,1))
	f2000_comp = util.get_flux_at_wavelength(w,f,2000)
	plot(w,f/f2000_comp*f2000, label="non-BAL composite",c=colors[2], linewidth=1.5)



	if angle == 72:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.hibal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		bal_select = (w < 3000)
		plot(w[bal_select],f[bal_select]/f2000_comp*f2000, label="HiBAL composite", c=colors[3], linewidth=1.5)

	elif angle == 75:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.lobal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		bal_select = (w < 3000)
		plot(w[bal_select],f[bal_select]/f2000_comp*f2000, label="LoBAL composite", c=colors[4], linewidth=1.5)

	elif angle == 78:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.felobal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		bal_select = (w < 3000)
		plot(w[bal_select],f[bal_select]/f2000_comp*f2000, label="FeLoBAL composite", c=colors[0], linewidth=1.5)


	# else:
	# # 	#names = ["nonbal","hibal", "lobal","felobal"]
	# 	names = ["nonbal","hibal", "lobal"]
	# 	labels = ["non-BAL composite", "hiBAL composite", "loBAL composite"]

	# 	for j in range(len(names)-1):
	# 		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.%s.spec" % names[j], unpack=True,usecols=(0,1))
	# 		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
	# 		plot(w,f/f2000_comp*f2000, label=labels[j], c=colors[2+j])

	start, end = gca().get_ylim()
	text(7100,start+(0.5*(end-start)),"$%i^\circ$" % angle, fontsize=16, rotation="vertical")
	xlim(1000,7000)
	ylim(1e-3,100)
	semilogy()
	#ylim(0,6.9)


	if i<4:
		gca().set_xticklabels([])
	if i == 2:
		ylabel("$F_{\lambda}$ at 10~Gpc ($10^{-14}$~erg s$^{-1}$ cm$^{-3}$ \AA$^{-1}$)", fontsize=20)

	if i != 1: float_legend(loc=4)

subplots_adjust(hspace=0)
xlabel("Wavelength (\AA)", fontsize=20)
savefig("wide_opt.png", dpi=500)
savefig("../../figures/wide_opt.png", dpi=500)
savefig("../../figures/wide_opt.eps", bbox="tight")

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

#set_pretty()
xkcd()
NORM = True
almost_black = '#262626'
#nameb = "nextgen_alpha1_long"
suffix = "_a075_pre_he"


namea = "grid11/nextgen%s" % suffix
namea = "/Users/jmatthews/Documents/runs/qso_he_test/with_estfix/nextgen_a05_pre_he"
#nameb = "../specs/grid11/nextgen_uv"
namec = "asi_thmin70_rmin50_a1_f0p01_mdotw5"
#namec = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"

name2 = "fiducial"
smoda = r.read_spectrum(namea)
#smodb = r.read_spectrum("../specs/"+nameb)
#smodc = r.read_spectrum("../specs/"+namec)


sdisk = r.read_spectrum("../specs/disk")
sh13 = r.read_spectrum("../specs/h13")

lines = [1215, 1240, 1400, 1550, 1640, 1860, 1909, 2800]
textx = [1215, 1400, 1550, 1640, 1860, 2800]
ycoords = [2, 50, 50, 50, 50, 50, 50]
line_labels = [r"Ly $\alpha$", r"N~\textsc{v}", r"Si~\textsc{iv}", r"C~\textsc{iv}", r"He~\textsc{ii}", r"Al~\textsc{iii}", "C~\\textsc{iii}]", r"Mg~\textsc{ii}"]
line_labels = [r"Ly $\alpha$~\&~N~\textsc{v}", r"Si~\textsc{iv}", r"C~\textsc{iv}", r"He~\textsc{ii}", r"Al~\textsc{iii}~\&~C~\textsc{iii}]", r"Mg~\textsc{ii}"]

angles = [20,60,72,75,78]

colors = get_colors()

figure(figsize=(10,8))

#run_agnspec

#for i in range(len(angles)):
for i in range(1):
	angle = angles[i]
	subplot(1,1,i+1)
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
	#plot(wwb,util.smooth(Fmodb), label=r"Model B", c=almost_black, linewidth=2)
	#plot(wwc,util.smooth(Fmodc), label=r"Model C, $\alpha = 1$", c=colors[2], linewidth=2)
	#plot(sh13["Lambda"],util.smooth(Fh13,window_len=35), label="H13", c=almost_black, linewidth=1)
	#plot(sdisk["Lambda"],util.smooth(Fdisk, window_len=100), label="Disc continuum", c="k", linestyle="--")

	f2000= util.get_flux_at_wavelength(wwa,Fmoda,2000)

	if angles[i]<=70:

		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.nonbal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="non-BAL composite",c=colors[2], linewidth=1.5)

	# if NORM:
	# 	f2000 = util.get_flux_at_wavelength(ww,Fmod,2000)
	# plot(ww,util.smooth(Fmod/f2000), label="Original Model")

	# read the composites and plot them
	# if angles[i]<=70:

	# 	w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.nonbal.spec", unpack=True,usecols=(0,1))
	# 	f2000 = util.get_flux_at_wavelength(w,f,2000)
	# 	plot(w,f/f2000, label="non-BAL composite")
	elif angle == 72 or angle == 75:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.hibal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="HiBAL composite", c=colors[3], linewidth=1.5)

	elif angle == 78:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.lobal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="LoBAL composite", c=colors[4], linewidth=1.5)

	elif angle == 80:
		w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.felobal.spec", unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="FeLoBAL composite", c=colors[0], linewidth=1.5)


	else:
	# 	#names = ["nonbal","hibal", "lobal","felobal"]
		names = ["nonbal","hibal", "lobal"]
		labels = ["non-BAL composite", "hiBAL composite", "loBAL composite"]

		for j in range(len(names)-1):
			w,f = np.loadtxt("/Users/jmatthews/Documents/Python/examples/sdssedrqsobal/sdssedrqsobal.%s.spec" % names[j], unpack=True,usecols=(0,1))
			f2000_comp = util.get_flux_at_wavelength(w,f,2000)
			plot(w,f/f2000_comp*f2000, label=labels[j], c=colors[2+j])

	start, end = gca().get_ylim()
	# if i == 0:
	# 	for l in range(len(textx)):
	# 		text(textx[l], ycoords[l], line_labels[l])

	# 	vlines(lines[:2], 8./60.0*end,0.25*end, linewidth=1)
	# 	vlines(lines[2:-1], 35./60.0*end,48/60.0*end, linewidth=1)
	# 	vlines(lines[-1:], 15/60.0*end,25/60.0*end, linewidth=1)
	# #else:
		#vlines(lines[:2], 8/60.0*end,15/60.0*end, linewidth=1)
	#	vlines(lines[2:-1], 35/60.0*end,48/60.0*end, linewidth=1)
	#	vlines(lines[-1:], 20/60.0*end,30/60.0*end, linewidth=1)

	start, end = gca().get_ylim()
	#text(3050,start+(0.5*(end-start)),"i=%i" % angle, fontsize=16, rotation="vertical")
	xlim(1000,3000)
	#ylim(0,6.9)
	title("Can a BAL outflow produce quasar-like emission lines?")

	annotate('Lyman alpha',xy=(1220,45), xycoords='data',xytext=(1400,50),textcoords='data',arrowprops=dict(arrowstyle="->"))
	annotate('C IV',xy=(1558,25), xycoords='data',xytext=(1630,34),textcoords='data',arrowprops=dict(arrowstyle="->"))
	annotate('Mg II',xy=(2795,8.5), xycoords='data',xytext=(2650,13),textcoords='data',arrowprops=dict(arrowstyle="->"))
	text(2200,40, "YES!", fontsize=60)
	text(2000,35, "(But they're a little on the puny side...)", fontsize=12)
	#if i<4:
		#gca().set_xticklabels([])
	#if i == 2:
	ylabel("Flux", fontsize=20)

	float_legend()

subplots_adjust(hspace=0, wspace=0)
xlabel("Wavelength", fontsize=20)
savefig("xk_uvspec_estfix.png", dpi=600)
#savefig("../../figures/uvspec.png", dpi=500, bbox="tight")
#savefig("../../figures/uvspec.eps", bbox="tight")

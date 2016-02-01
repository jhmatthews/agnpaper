#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import pyfits as py
import sdss
from pretty import *


set_pretty()
big_tick_labels(18)
rcParams["axes.linewidth"]=1.5
NORM = True
almost_black = '#262626'
#nameb = "nextgen_alpha1_long"
suffix = "_a075_pre_he"

PYTHON = os.environ["PYTHON"]


namea = "grid11/nextgen%s" % suffix
namea = "../specs/webgrid/run5_thmin70_rmin50_a0p5_rv1e19_f0p01.spec"
namea = "../../../runs/referee_m16/rerun_adiab/run42_thmin70_rmin50_a0p5_rv1e19_f0p01.spec"
#nameb = "../specs/grid11/nextgen_uv"
namec = "asi_thmin70_rmin50_a1_f0p01_mdotw5"
#namec = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"

name2 = "fiducial"
smoda = r.read_spectrum(namea)
#smodb = r.read_spectrum("../specs/"+nameb)
#smodc = r.read_spectrum("../specs/"+namec)


sdisk = r.read_spectrum("../../generate/specs/disk")
sh13 = r.read_spectrum("../../generate/specs/h13")

lines = [1215, 1240, 1400, 1550, 1640, 1860, 1909, 2800]
textx = [1215, 1400, 1550, 1650, 1860, 2800]
ycoords = [2, 50, 50, 50, 50, 30, 10]
line_labels = [r"Ly $\alpha$", r"N~\textsc{v}", r"Si~\textsc{iv}", r"C~\textsc{iv}", r"He~\textsc{ii}", r"Al~\textsc{iii}", "C~\\textsc{iii}]", r"Mg~\textsc{ii}"]
line_labels = [r"Ly $\alpha$~\&~N~\textsc{v}", r"Si~\textsc{iv}", r"C~\textsc{iv}", r"He~\textsc{ii}", r"Al~\textsc{iii}~\&~C~\textsc{iii}]", r"Mg~\textsc{ii}"]

angles = [20,60,73,76]

colors = get_colors()

figure(figsize=(10,12))

NORM_WAVELENGTH = 1300


def get_hst(f2000):


	files1 = [f for f in os.listdir('0946_fits/') if "x1d" in f]
	z = 1.223
	for i in range(4,len(files1)):

		d = py.getdata('0946_fits/'+files1[i])
		w = d["WAVELENGTH"][0] / (1.0 + z)
		f = d["FLUX"][0]

		f2000_comp = util.get_flux_at_wavelength(w,f,NORM_WAVELENGTH )

		
	ff = f/f2000_comp*f2000
	plot(w[:-10],ff[:-10], label="PG0946+301 STIS", c=colors[0], linewidth=2, alpha=1)

	d = py.getdata('0946_fits/o5ao06020_sx1.fits')
	w = d["WAVELENGTH"][0] / (1.0 + z)
	f = d["FLUX"][0]
	f2000_comp = util.get_flux_at_wavelength(w,f,NORM_WAVELENGTH )
	ff = f/f2000_comp*f2000
	plot(w[:-10],ff[:-10], c=colors[0], linewidth=2, alpha=1)

	# name = "spSpec-53431-1947-408.fit"
	# w, flux = sdss.read_sdss(name)

	# w = w / (1.0 + z)

	# f2000_comp = util.get_flux_at_wavelength(w,flux,2000)
	#plot(w,util.smooth(flux/2/f2000_comp*f2000), c=colors[0], linewidth=2, alpha=1)

	return 0

#run_agnspec

for i in range(len(angles)):

	angle = angles[i]
	subplot(4,1,i+1)
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




	f2000= util.get_flux_at_wavelength(wwa,Fmoda,2000)

	if angles[i]<=70:

		w,f = np.loadtxt("%s/examples/sdssedrqsobal/sdssedrqsobal.nonbal.spec" % PYTHON, unpack=True,usecols=(0,1))
		f2000_comp = util.get_flux_at_wavelength(w,f,2000)
		plot(w,f/f2000_comp*f2000, label="non-BAL composite",c=colors[2], linewidth=1.5, alpha=1)

	elif angle > 70 and angle < 75:
		f1300 = util.get_flux_at_wavelength(wwa,Fmoda,1300)
		get_hst(f1300)


	elif angle > 75:
		name = "spSpec-52026-0593-509.fit"
		name = "spSpec-52057-0626-423.fit"
		redshift = py.getval(name,'z',0)


		w, flux = sdss.read_sdss(name)

		w = w / (1.0 + redshift)

		f2000_comp = util.get_flux_at_wavelength(w,flux,2000)
		plot(w,util.smooth(flux/f2000_comp*f2000), label="SDSS J162901.63+453406.0", c=colors[1], linewidth=2, alpha=1)


	#ylim(0,6.9)
	if i == 2: ylim(0,6.8)
	if i == 3: ylim(0,3.4)
	if i == 1: ylim(0,24.9)
	start, end = gca().get_ylim()
	if i == 0:
		for l in range(len(textx)):
			text(textx[l], ycoords[l], line_labels[l])

		vlines(lines[:2], 8./60.0*end,0.25*end, linewidth=1)
		vlines(lines[2:-1], 35./60.0*end,48/60.0*end, linewidth=1)
		vlines(lines[-1:], 15/60.0*end,25/60.0*end, linewidth=1)


	start, end = gca().get_ylim()
	text(3050,start+(0.5*(end-start)),"$%i^\circ$" % angle, fontsize=20, rotation="vertical")
	xlim(800,3000)

	plot(wwa,util.smooth(Fmoda,window_len=5), label=r"Model", c='k', linewidth=2)
	plot(sdisk["Lambda"],util.smooth(Fdisk, window_len=100), label="Disc continuum", c="k", linestyle="--")


	if i<3:
		gca().set_xticklabels([])
	#if i == 2:
	#	ylabel("$F_{\lambda}$ at 10~Gpc ($10^{-14}$~erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)", fontsize=20)

	float_legend()

text(600, 10.05,"$F_{\lambda}$ at 10~Gpc ($10^{-14}$~erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)", rotation="vertical", fontsize=24)
subplots_adjust(hspace=0, wspace=0, left=0.1, right=0.95, top=0.97, bottom=0.07)
xlabel("Wavelength (\AA)", fontsize=24)
savefig("uvspec_new.png", dpi=100)
#savefig("../../figures/uvspec.png", dpi=500, bbox="tight")
savefig("../../figures/fig3.eps", bbox="tight")

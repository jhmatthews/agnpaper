#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *


figure(figsize=(12,6))

#subplot(2,1,1)

#fnames = ["1e45", "1e45_uv"]
#fnames = ["nextgen_alpha1_long_wide","1e45", "1e45_uv"]
#fnames = ['nextgen_p3_3',"nextgen_3_30", 'nextgen_30_300', 'nextgen_uv', 'nextgen_opt']
#fnames = ['checkvol_a6',"checkvol_a6_xray","h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
#fnames = ["h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
fnames = ["../specs/webgrid/run5_thmin70_rmin50_a0p5_rv1e19_f0p01.spec"]
NN = len(fnames) - 1
NN = len(fnames) - 1

alphas = [0.5,0.5,0.75,0.75]
labels=["1e45, f=0.01", "H13"]
big_tick_labels(18)

lines = np.log10(C/np.array([6563, 4861, 1215, 1240, 1400, 1550, 1909, 2800])/ANGSTROM)
textx = np.log10(C/np.array([6563, 4861, 1215, 1400, 1550, 1909, 2800])/ANGSTROM)
#ycoords = [2, 50, 50, 50, 50, 10, 10]
line_labels = [r"H$\alpha$", r"H$\beta$", r"Ly $\alpha$~\&~N~\textsc{v}", r"Si~\textsc{iv}", r"C~\textsc{iv}", r"C~\textsc{iii}]", r"Mg~\textsc{ii}"]


colors = get_colors()

ccc = [colors[0],colors[0], "k", "k"]
modcols = [colors[0], "k", "k"]

angles = [20, 89]
geolab = ["`Quasar-like'", "`Below Wind'"]

lims = ((30.5,31.7), (27.5,30.1))

nu_norm = C/2000.0/ANGSTROM

for iang in range(len(angles)):

	#subplot(2,1,iang+1)

	set_pretty()

	all_nu = np.array([])
	all_f = np.array([])
	concatenate

	#bounds = [(1e19,1e18),(1e18,1e17),(1e17,1e16),(1e16,1e15),(1e15,1e14)]
	bounds = [(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14)]

	for i in range(len(fnames)):


		s = r.read_spectrum(fnames[NN-i])

		ev = 1.60217657e-19

		angle = angles[iang]

		nu = s["Freq."]
		w = s["Lambda"]

		flambda = s["A%iP0.50" % angle]


		fl_to_fnu = s["Lambda"] / s["Freq."]

		norm = 4.0 * PI * (100.0 * PARSEC) * (100.0 * PARSEC)

		lnu = flambda * fl_to_fnu * norm

		f2000 = util.get_flux_at_wavelength(nu, lnu,nu_norm)

		w,f = np.loadtxt("/Users/jamesmatthews/Documents/Python/examples/sdssedrqsobal/sdss_hst_combined.dat", unpack=True,usecols=(0,1))
		nuc = C / (w*ANGSTROM)
		fl_to_fnu = w / nuc 
		lnuc = f * fl_to_fnu

		f2000_comp = util.get_flux_at_wavelength(nuc, lnuc,nu_norm)
		lnuc = lnuc/f2000_comp

		if iang == 0: 
			plot(np.log10(nuc), np.log10(lnuc), label="non-BAL composite",c=colors[2], linewidth=2)
			text(14.8, 0.5, "non-BAL composite", color=colors[2], fontsize=20)
		plot(np.log10(nu), util.smooth(np.log10(lnu/f2000)-(2+((1-iang)*2)), window_len = 5), c=colors[iang], label="$%i^\circ$" % angles[iang], linewidth=2)
		
		dd = 0.2*(1 - iang) #little offset thing
		text(14.8, 0.5-dd-(2+((1-iang)*2)), "%s, $%i^\circ$" % (geolab[iang], angles[iang]), color=colors[iang], fontsize=20)

		if iang == 0:
			for ii in range(len(textx)):
				text(textx[ii], -0.5, line_labels[ii])
			vlines(lines, -0.8,-0.6, linewidth=1.5)
		#float_legend(loc=3)


		xlim(14.64,15.5)
		ylim(-4.6,1.2)
		#ylim(lims[iang][0], lims[iang][1])
		ylabel(r"$\log [F_\nu / F_{2000}]$ (Arb.)", fontsize=20)

		#if iang == 0:
		#	gca().set_xticklabels([])

		big_tick_labels(18)
		long_ticks()

hlines([0,-2,-4],14.64,np.log10(nu_norm), linestyle="--")
#ylim(1e41,1e45)
subplots_adjust(hspace=0, wspace=0, left=0.08, right=0.97)

xlabel(r"$\log \nu$", fontsize=20)

savefig("sed.png",dpi=300)
savefig("../../figures/fig5.eps")
#sys.exit()




#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *


figure(figsize=(7,10))

#fnames = ["1e45", "1e45_uv"]
#fnames = ["nextgen_alpha1_long_wide","1e45", "1e45_uv"]
#fnames = ['nextgen_p3_3',"nextgen_3_30", 'nextgen_30_300', 'nextgen_uv', 'nextgen_opt']
#fnames = ['checkvol_a6',"checkvol_a6_xray","h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
#fnames=["nextgen_a05_pre_xray", "nextgen_a06_pre_xray", "nextgen_a075_pre_xray", "nextgen_a05_pre_xray_he","nextgen_a06_pre_xray_he"]
fnames = ["nextgen_a05_pre_xray_he", "nextgen_a05_pre_xray"]

#fnames = ["h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
NN = len(fnames) - 1
NN = len(fnames) - 1

alphas = [0.5,0.5,0.75,0.75]
labels=["1e45, f=0.01", "H13"]
big_tick_labels(16)

colors = get_colors()

#ccc = [colors[0],colors[0], "k", "k"]
#


angles = [60,72,78,75]

for iang in range(len(angles)):

	subplot(4,1,iang+1)

	set_pretty()

	all_nu = np.array([])
	all_f = np.array([])
	
	text(2e17,1e44,"$i=%2i^\circ$" % angles[iang], fontsize=20)

	#bounds = [(1e19,1e18),(1e18,1e17),(1e17,1e16),(1e16,1e15),(1e15,1e14)]
	#bounds = [(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14)]

	for i in range(len(fnames)):

		lname = fnames[i][8:11]
		if "he" in fnames[i]:
			lname ="He and H"
		else:
			lname = "H"


		s = r.read_spectrum("../specs/grid11/"+fnames[i])

		ev = 1.60217657e-19

		angle = angles[iang]

		nu = s["Freq."]
		w = s["Lambda"]

		flambda = s["A%iP0.50" % angle]


		fl_to_fnu = s["Lambda"] / s["Freq."]

		norm = 4.0 * PI * (100.0 * PARSEC) * (100.0 * PARSEC)

		lnu = flambda * fl_to_fnu * norm


		#all_nu = np.concatenate( (all_nu, nu[select]) )
		#all_f = np.concatenate( (all_f, nu[select]*lnu[select])) 

		#x = np.arange(14.5,18,0.05)
		#xx = np.interp(x, np.log10(all_nu), np.log10(all_f))
		plot(nu, util.smooth(nu*lnu, window_len = 10), label=lname, linewidth=2)
		loglog()
		float_legend()
		xlim(1e17,2e18)
		ylim(1e41,1e45)
		vlines([2000.0/HEV],1e41,1e45)

	ylabel(r"$\nu L_\nu$ (erg~s$^{-1}$)", fontsize=20)



big_tick_labels(16)
long_ticks()







xlabel(r"$\nu$ (Hz)", fontsize=20)


# ax2 = gca().twiny()
# ev = 10.0**np.arange(-3,2,1)
# new_tick_locations = (ev*1000.0) / HEV

# def tick_function(X):
#     V = 1/(1+X)
#     return ["%.3f" % z for z in V]

# #semilogx()
# ax2.set_xscale("log")
# ax2.set_xlim(1e14,1e19)
# ax2.set_xticks(new_tick_locations)
# print new_tick_locations
# ax2.set_xticklabels(ev)
# ax2.set_xlabel("E (keV)", fontsize=20)


savefig("xray_angles2.png",dpi=200)
#savefig("../../figures/sed_all_balqsos.eps",bbox="tight")
clf()


	# set_pretty()
	# plot(HEV*nu/1000.0, util.smooth(lnu) )

	# semilogy()
	# xlim(0.2,10)

	# xlabel(r"Energy (keV)", fontsize=20)
	# ylabel(r"$L_\nu$ (erg~s$^-1$~Hz$^-1$)", fontsize=20)
	# savefig("xray_spectrum_lnu_ev.png")


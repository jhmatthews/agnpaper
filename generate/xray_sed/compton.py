#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *


figure(figsize=(12,7))

#subplot(2,1,1)

#fnames = ["1e45", "1e45_uv"]
#fnames = ["nextgen_alpha1_long_wide","1e45", "1e45_uv"]
#fnames = ['nextgen_p3_3',"nextgen_3_30", 'nextgen_30_300', 'nextgen_uv', 'nextgen_opt']
#fnames = ['checkvol_a6',"checkvol_a6_xray","h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
#fnames = ["h13","fiducial_xray_fine","fiducial_uv","fiducial_xray"]
fnames = ["grid11/nextgen_a05_pre_xray"]
NN = len(fnames) - 1
NN = len(fnames) - 1

alphas = [0.5,0.5,0.75,0.75]
labels=["1e45, f=0.01", "H13"]
big_tick_labels(16)

colors = get_colors()

ccc = [colors[0],colors[0], "k", "k"]
modcols = [colors[0], "k", "k"]

angles = [85]

for iang in range(len(angles)):

	set_pretty()

	all_nu = np.array([])
	all_f = np.array([])
	concatenate

	#bounds = [(1e19,1e18),(1e18,1e17),(1e17,1e16),(1e16,1e15),(1e15,1e14)]
	bounds = [(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14),(1e19,1e14)]

	for i in range(len(fnames)):

		s = r.read_spectrum("../specs/"+fnames[NN-i])

		ev = 1.60217657e-19

		angle = angles[iang]

		nu = s["Freq."]
		w = s["Lambda"]

		flambda = s["A%iP0.50" % angle]


		fl_to_fnu = s["Lambda"] / s["Freq."]

		norm = 4.0 * PI * (100.0 * PARSEC) * (100.0 * PARSEC)

		lnu = flambda * fl_to_fnu * norm

		select = (nu < bounds[NN-i][0]) * (nu > bounds[NN-i][1]) 

		all_nu = np.concatenate( (all_nu, nu[select]) )
		all_f = np.concatenate( (all_f, nu[select]*lnu[select])) 

		x = np.arange(14.5,18,0.05)
		xx = np.interp(x, np.log10(all_nu), np.log10(all_f))
		plot(nu*HEV*1e-3, util.smooth(lnu, window_len = 20), c=ccc[i], label="Model", linewidth=2)



big_tick_labels(16)
long_ticks()


f_2kev = 2000.0/HEV
f_2500 = C/(2500.0 / 1e8)
vlines([f_2kev, f_2500],0.2e47,1e47,linestyle="-", linewidth=1.2)
text(f_2kev*1.1,5e46,"$2$~keV")
text(f_2500*1.1,5e46,"$2500$~\AA")

#lines = 1e-3 * HEV * C / (lines * ANGSTROM)
#vlines(lines, 1e41,1e45)

float_legend(loc=2)
#semilogy()
loglog()
#xlim(5e14,1e18)
xlim(0.4,10)
ylim(1e23,1e27)

xlabel(r"$E$ (keV)", fontsize=20)
ylabel(r"$\nu L_\nu$ (erg~s$^{-1}$)", fontsize=20)



savefig("sed_balqsos_data.png",dpi=300)
#sys.exit()




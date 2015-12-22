#!/usr/bin/env python
'''
alpha_ox.py makes plots of alpha_ox and L_2kev plotted
against L_2500 
'''
from pylab import *
import py_read_output as r 
import numpy as np 
import os, sys
import py_plot_util as util
from constants import *
from astropy.io import ascii
from pretty import *
#import disk

# some functions for getting mean monochromatic fluxes in a bin
def get_flux_inwbin(lambda_array, flux_array, w, binsize=1, log=False):
	wstart = w - (0.5*binsize)
	wstop = w + (0.5*binsize)
	select = (lambda_array > wstart) * (lambda_array < wstop)
	return np.mean(flux_array[select])

def get_flux_infbin(f_array, flux_array, f, binsize=0.1, log=False):

	binsize *= f

	fstart = f - (0.5 * binsize)
	fstop = f + (0.5 * binsize)
	select = (f_array > fstart) * (f_array < fstop)
	return np.mean(flux_array[select])

def get_flux_at_f(f_array, flux_array, f):

	i = np.abs(f_array - f).argmin()

	return flux_array[i]



r.setpars()
rcParams["legend.frameon"] = "False"
rcParams["legend.numpoints"] = 1
colors = ["r", "g"]
colors_xx = get_colors()

# fnames to use
# fnames = [("fiducial_xray", "fiducial_uv"), ("../specs/nextgen_a5_xray", "../specs/nextgen_a5")]
# fnames = [("fiducial_xray", "fiducial_uv"), ("../specs/nextgen_3_30", "../specs/nextgen_uv")]
# fnames = [("../specs/fiducial_xray_fine", "../specs/h13.spec"), ("../specs/checkvol_a6_xray", "../specs/checkvol_a6")]

aaa = "06"
vers = "pre"
suffix1 = "_a%s_%s" % (aaa,vers)
suffix2 = "_a%s_%s_xray" % (aaa,vers)
fnames = [("../specs/fiducial_xray_fine", "../specs/h13.spec"), ("restart_xrays", "/Users/jmatthews/Documents/runs/qso_he_test/with_estfix/nextgen_a05_pre_he")]
#fnames = 
# 9 before
label=fnames[1][0][16:]


# load the Steffen et al. data
dat = np.loadtxt("steffen.dat", usecols=(10,12,13), dtype="string", unpack=True)
dat2 = np.loadtxt("steffen_bqs.dat", usecols=(10,12,13), dtype="float", unpack=True)

islim = np.empty(len(dat[1]), dtype=bool)

for i in range(len(dat[1])):
	islim[i] = ("<" in dat[1][i])

l2500 = dat[0][~islim].astype(np.float)
lx = dat[1][~islim].astype(np.float)
alphaox = dat[2][~islim].astype(np.float)



L_2500s= []
alpha_oxs = []
angles = []
markers = []
Lxs = []

mode = 'both'

ang = [20,60,72,75,78]

xray_do = [(1,0)]



def get_lx(s, name):

	from scipy.integrate import trapz

	L_05_2 = trapz(s["Lambda"], s[name])

	print L_05_2

	return 0


big_tick_labels(20)
long_ticks()
for ff in fnames:

	s1 = r.read_spectrum(ff[0])
	s2 = r.read_spectrum(ff[1])

	ncomponents = 9
	nspecs = len(s1.dtype.names) - ncomponents

	lambda_2kev = 6.2

	#bin1 = 

	#figure(figsize=(8,12))


	norm = 4.0 * PI * (100.0 * PARSEC) * (100.0 * PARSEC)

	fl_to_fnu1 = s1["Lambda"] / s1["Freq."]
	fl_to_fnu2 = s2["Lambda"] / s2["Freq."]

	L_2500s_temp = []
	alpha_oxs_temp = []
	angles_temp = []
	markers_temp = []
	Lxs_temp = []

	for i in range(len(ang)):

		#subplot(4,2,i+1)

		name = s1.colnames[ncomponents + i]
		an = ang[i]
		if "checkvol" in ff[0]:
			#print "YAS"
			if an == 70:
				an = 72
			if an == 75:
				an = 75
			if an == 80:
				an = 78
		name = "A%iP0.50" % an

		fnorm1 = s1[name] * fl_to_fnu1
		fnorm2 = s2[name] * fl_to_fnu2




	 	#f_2kev = util.get_flux_at_wavelength(s1["Lambda"],fnorm1, 6.2) * norm
	 	f_2kev = get_flux_infbin(s1["Freq."], fnorm1, 2000.0/HEV) * norm
	 	#f_2kev = get_flux_at_f(s1["Freq."], fnorm1, 2000.0/HEV) * norm

	 	#f_2kev = get_lx (s1, name)

	 	#f_2500 = util.get_flux_at_wavelength(s2["Lambda"],fnorm2, 2500.0) * norm
	 	nu_2500 = C / (2500.0 * ANGSTROM)
	 	#f_2500 = get_flux_inwbin(s2["Lambda"], fnorm2, 2500.0) * norm
	 	f_2500 = get_flux_infbin(s2["Freq."], fnorm2, nu_2500) * norm
	 	#f_2500 = get_flux_at_f(s2["Freq."], fnorm2, nu_2500) * norm

	 	alpha_ox = 0.3838 * np.log10(f_2kev/f_2500)


	 	L_2500s_temp.append(f_2500)
	 	Lxs_temp.append(f_2kev)
	 	alpha_oxs_temp.append(alpha_ox)

	 	angle = int(name[1:3])
	 	angles_temp.append(angle)

	 	if angle > 80.0:
	 		marker="*"
	 	elif angle < 70:
	 		marker = "o"
	 	else:
	 		marker = "x"

	 	markers_temp.append(marker)

	 	#ylim(0.001,100)

	L_2500s.append(L_2500s_temp)
	alpha_oxs.append(alpha_oxs_temp)
	angles.append(angles_temp)
	markers.append(markers_temp )
	Lxs.append(Lxs_temp)

bal_l2500, bal_l2, bal_alphaox, bal_alphaoxerror = np.loadtxt('balqsos.dat', usecols=(1,3,7,8), unpack=True)



figure(figsize=(12,8))
#suptitle(titles[icolor])

#axes.color_cycle    : 332288, 88CCEE, 44AA99, 117733, 999933, DDCC77, CC6677, 882255, AA4499
colors = ["#332288", "#88CCEE", "#44AA99", "#117733"]
colors = [colors_xx[0], colors_xx[2]]


bbox_props = [dict(boxstyle="square,pad=0.2", fc="None", ec=colors[0], lw=2), dict(boxstyle="square,pad=0.2", fc="None", ec=colors[1], lw=2)]


#clf()

scatter(l2500, lx, c="k", s=80, label="Steffen+ 2006, COMBI-7")
scatter(dat2[0].astype(np.float), dat2[1].astype(np.float), edgecolors="k", s=80, facecolors="None", label="Steffen+ 2006, BQS")
scatter(bal_l2500, bal_l2, marker="^", edgecolors=colors_xx[1], facecolors="None", s=180, lw=2.5, label="Saez+ 2012, BALQSO")

for jj in range(2):
	for i in range(len(angles[jj])):
		#text(L_2500s[jj][i], Lxs[jj][i], angles[jj][i], color=colors[jj],fontsize=16)

		xy = np.array([L_2500s[jj][i], Lxs[jj][i]])
		xytext = np.array([L_2500s[jj][i], Lxs[jj][i]])
		xytext[1] += 2e27
		xy[1] += 0.5e27

		#annotate(angles[jj][i], xy, xytext, color=colors[1],fontsize=16, arrowprops=dict(arrowstyle="-"))
		#if jj == 0 or angles[jj][i]!=60:
		text(np.log10(L_2500s[jj][i]-(L_2500s[jj][i]*0.13)), np.log10(Lxs[jj][i]-(Lxs[jj][i]*0.15)), angles[jj][i], color=colors[jj],fontsize=16, bbox=bbox_props[jj])
		#else:
		#	scatter(L_2500s[jj][i], Lxs[jj][i], edgecolors=colors[jj], facecolors="None", marker="s", s=400, lw=2)
		#	annotate(angles[jj][i], xy, xytext, color=colors[jj],fontsize=16, arrowprops=dict(arrowstyle="-", edgecolor="k", connectionstyle="arc3"))

logl2500_line = np.arange(27.5,33.5,0.5)
loglx_line = (0.721 * logl2500_line) + 4.531

plot(logl2500_line, loglx_line, linestyle="--", c="k")

#scatter(L_2500s[0],Lxs[0], edgecolors=colors[0], facecolors="None",marker="s", s=400, label="H13", lw=2)
#scatter(L_2500s[1],Lxs[1], edgecolors=colors[1], facecolors="None",marker="s",  s=400, label="This work", lw=2)

scatter([0,0],[0,0], edgecolors=colors_xx[0], facecolors="None",marker="s", s=100, label="H13", lw=2)
scatter([0,0],[0,0], edgecolors=colors[1], facecolors="None",marker="s",  s=100, label="This work", lw=2)

#loglog()
ylim(22,28.4)
xlim(27.5,32.5)


xlabel(r"$\log[L_{2500}$~(erg~s$^{-1}$~Hz$^{-1}$)]", fontsize=20)
ylabel(r"$\log[L_{2kev}$~(erg~s$^{-1}$~Hz$^{-1}$)]", fontsize=20)

subplots_adjust(hspace=0)
long_ticks()
legend(loc=2, scatterpoints = 1)

savefig("lx%s.png" % suffix1)
#savefig("lx.png", bbox="tight")
#
savefig("../../figures/lx.png", bbox="tight")
#savefig("../../figures/alpha_ox_%s.eps" % label, bbox="tight")
clf()

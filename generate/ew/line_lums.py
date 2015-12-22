#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import pretty


def get_continuum(wave, spec, lines = [1215,1400,1550,1640,2800,1900], lmin=1000, lmax = 2000, deg = 3):

	i_jump = (wave > lmin) * (wave < lmax)

	i_lines = np.zeros(len(wave), dtype=bool)

	for i in range(len(wave)):

		comp = lines - wave[i]
		comp = np.absolute(comp)

		i_line = (comp < 100)

		if np.sum(i_line) == 0:
			i_lines[i] = True
		else:
			i_lines[i] = False		

	i_select = i_jump * i_lines

	continuum_points = wave[i_select]

	#print continuum_points, i_select

	coefficients = np.polyfit(continuum_points, spec[i_select], deg)

	return coefficients

def make_f(a, x):

 	deg = len(a)
 	y = 0

 	for i in range(deg):
 		y += a[i] * (x**(deg - i - 1))

 	return y

def find_line(w,f, wline = 6563, window_size=100, rw = False):

	fabs_w = np.fabs(w - wline)

	if rw:
		w_bool = (fabs_w < window_size) * (w > wline)
	else:
		w_bool = (fabs_w < window_size) 

	# create windows to look in
	w_window = w
	f_window = w_bool * f

	# return arrays with length over specified 2*window_size
	return w_window[f_window >0], f_window[f_window >0]


def normalize(w, f, lmin=1300, lmax=3000):

	a = get_continuum(w, f, [1215,1550,1400,2800,6563], lmin=lmin, lmax=lmax, deg=2)

	#fnorm = f/make_f(a, w)

	#fnorm = f / util.get_flux_at_wavelength(w,f,1450)

	return make_f(a, w)



l_3000,l_1350,l_ha, ews_c4, l_c4, ews_mg2, l_mg2, bal_flag, r_flag, mbh, edd = np.loadtxt("/Users/jmatthews/Documents/J_ApJS_194_45.tar/catalog.dat", unpack=True, usecols=(20,22,24, 107,103,87,83, 13, 14, 136,140))

fname = "/Users/jmatthews/Documents/runs/qso_he_test/with_estfix/nextgen_a05_pre_he"

s = r.read_spectrum(fname)


lines = [1550,2800, 6563]
ll = [l_c4, l_mg2, l_ha]
labels=["c4", "Mg", "Ha"]
pretty.set_pretty()
colors=pretty.get_colors()
figure(figsize=(8,10))
pretty.long_ticks()
pretty.big_tick_labels(18)

for il in range(len(ll)):
	w = s["Lambda"]
	nu = s["Freq."]

	select = (mbh > 8.5) * (mbh < 9.5) * (edd < 0) * (edd > -1.5)
	select_qsos = select * (ll[il] > 0.0)

	line_lums = []
	incs = []

	angles = ["A20P0.50", "A73P0.50"]

	for i in range(len(angles)):
		angle = angles[i]
		f = s[angle]

		if labels[il] != "Ha":
			f_continuum = normalize(w,f)
			ww,ff=find_line(w,f-f_continuum,wline=lines[il],rw=True)

		else:
			f_continuum = normalize(w,f,lmin=3000, lmax=7000)
			ww,ff=find_line(w,f-f_continuum,wline=lines[il],rw=True)


		l_line = 2.0*np.fabs(np.trapz(ff,x=ww)) * 4.0 * PI * (100.0 * PARSEC) **2

		flambda_to_fnu = s["Lambda"] / s["Freq."]

		l3000 = C/(3000.0*ANGSTROM)*util.get_flux_at_wavelength(w,flambda_to_fnu*f,3000 ) * 4.0 * PI * (100.0 * PARSEC) **2
		l1350 = C/(1350.0*ANGSTROM)*util.get_flux_at_wavelength(w,flambda_to_fnu*f,1350) * 4.0 * PI * (100.0 * PARSEC) **2 

		print angle, labels[il]
		print np.log10(l_line), np.mean(ll[il][select_qsos]), np.std(ll[il][select_qsos])
		print np.log10(l3000), np.mean(l_3000[select*(l_3000 > 0)]), np.std(l_3000[select*(l_3000 > 0)])
		print np.log10(l1350), np.mean(l_1350[select*(l_1350 > 0)]), np.std(l_1350[select*(l_1350 > 0)])

		line_lums.append(l_line)
		incs.append(int(angle[1:3]))

	subplot(211)
	hist(ll[il][select_qsos], orientation="horizontal", bins=np.arange(42,46,0.1), facecolor=colors[il],alpha=0.5)
	ylim(42,46)
	xlabel("Counts", fontsize=20)
	ylabel(r"$\log [L_{line} ($erg~s$^{-1}$)]", fontsize=20)

	subplot(212)
	plot(incs, np.log10(line_lums), c=colors[il], linewidth=2, label=labels[il])
	ylim(42,46)
	xlabel(r"Inclination ($^\circ$)", fontsize=20)
	ylabel(r"$\log [L_{line} ($erg~s$^{-1}$)]", fontsize=20)

	#print l_line, 

pretty.float_legend()
suptitle("Line Luminosity Comparison: Model, v SDSS Quasar catalog\n($8.5<\log M_{BH}<9.5$, $-1.5<\log (L/L_{edd}) < 0$)", fontsize=20)
savefig("line_lums.png",dpi=200)
# w = s["Lambda"]
# for i in range(len(s.colnames[9:])):
# 	angle_string = s.colnames[9 + i]
# 	f = s[angle_string]

# 	f_continuum = normalize(w,f)



#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import pretty
from pretty import *

D = 100.0 * PARSEC

def Ledd (m):

	'''
	calculates eddington luminosity for a solar mass m
	
	Args:
		m	mass in solar masses
	
	returns:
		Eddington luminosity in ergs s^-1, float
	'''
	
	m *= MSOL
	
	consts = (4.0 * PI * G * C * MPROT ) / THOMPSON
	
	L = consts * m
	
	return L

def run_agnspec(fname, mass=1e9, mdot=5, angm=0, alpha=0.01, inc=80):

	mu = np.cos(inc / 180.0 * np.pi)

	print inc, mu

	param_string = "mass=%8.4e,mdot=%8.4e,angm=%.6f,alpha=%.6f,mu=%8.4e,savename=\"%s\"" % (mass, mdot, angm, alpha, mu, fname)

	print param_string
	os.system("idl -e 'agnspec,%s' > idl.out" % param_string)

	return 0

def read_agnspec(fname):

	nu, lnu = np.loadtxt(fname, unpack=True, usecols=(0,1))

	wave = (C / nu) / ANGSTROM

	llambda = nu * lnu / wave

	flambda = llambda / 4.0 / PI / (D*D)

	return wave, flambda

def mdot_from_edd ( edd, m , eta = 0.1):

	''' 
	calculates an accretion rate from an eddington fraction
	Args:
		edd_frac		eddington fraction
		m			mass of central object in solar masses
		eta			accretion efficiency, set to 1 (0.1 more realistic)
	
	returns:
		mdot in solar masses / yr
	'''
	
	L = Ledd (m)		# eddington luminosity
	
	mdot = edd * L / ( (C ** 2) )
	
	mdot *= 1.0 / eta
	
	mdot = mdot * ( YR ) / MSOL	# normalise units
	
	return mdot

def run_grid(masses, edd_frac, spin, incs):

	for ispin in range(len(spin)):
		for im in range(len(masses)):
			for iedd in range(len(edd_frac)):

				mdot = mdot_from_edd(edd_frac[iedd], masses[im])

				# figure(figsize=(8,8))

				# pretty.set_pretty()

				for ii in range(len(incs)):
					fname = "sp%i_inc%2i_mbh%i_edd%itminus2" % (int(spin[ispin]+0.5),incs[ii], np.log10(masses[im]), edd_frac[iedd]/0.01)

					run_agnspec(fname, mass=masses[im], mdot=mdot, inc=incs[ii], angm=spin[ispin])

def set_set2():
	import brewer2mpl
	color = "Set1"
	set2 = brewer2mpl.get_map(color, 'qualitative', 8).mpl_colors
	g = set2[2]
	b = set2[1]
	r = set2[0]
	set2[0] = b
	set2[1] = g
	set2[2] = r

	return set2


set_pretty()
edd_fracs = [0.2]
spins = [0,0.998]
incs = np.arange(10,90,10)
lstyles=["-", "--"]
set2 = set_set2()
figure(figsize=(8,8))
big_tick_labels(14)

lls = [500,1000,2000]

for il in range(len(lls)):
	subplot(3,1,il+1)
	flambda_look = lls[il]
	for iedd in range(len(edd_fracs)):
		for ispin in range(len(spins)):

			fl = []
			for ii in range(len(incs)):

				fname = "sp%i_inc%i_mbh%i_edd%itminus2" % (int(spins[ispin]+0.5),incs[ii], np.log10(1e9), edd_fracs[iedd]/0.01)


				w,f=read_agnspec(fname)
				ff = util.get_flux_at_wavelength(w,f,flambda_look)

				fl.append(ff)

			fl= np.array(fl)
			plot(incs, fl/fl[0], linewidth=2, linestyle=lstyles[ispin], c=set2[ispin], label="$\epsilon=%.2f, a_* = %.3f$" % (edd_fracs[iedd], spins[ispin]) )



	s = r.read_spectrum('../specs/disk_wide')
	fl = []

	for ii in range(len(incs)):
		f = s["A%iP0.50" % incs[ii]]

		fl.append(util.get_flux_at_wavelength(s["Lambda"],f,flambda_look))

	plot(incs, fl/fl[0], linewidth=2, c=set2[2], label="Classical AD")

	ylim(0,1.5)

	if il == 1: ylabel("$F_{\lambda}$ / $F_{\lambda} (i=0^\circ)$", fontsize=20)
	text(15,0.4,"$\lambda = %i$\AA" % flambda_look, fontsize=20)

	if il<2:
		gca().set_xticklabels([])
	#long_ticks()



subplots_adjust(hspace = 0)
xlabel("$i~(^\circ$)", fontsize=20)
#ylabel("$F_{%i}$ / $F_{%i} (i=0^\circ)$" % (flambda_look,flambda_look), fontsize=20)
float_legend()


savefig("agnspec.png", dpi=300)
savefig("../../figures/agnspec.png", dpi=300)







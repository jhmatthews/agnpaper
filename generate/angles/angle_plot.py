#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from collections import OrderedDict
from plot_norm import make_f, get_continuum
import plot_norm as _plot_norm
from pretty import *
try:
	from thread import get_ident as _get_ident
except ImportError:
	from dummy_thread import get_ident as _get_ident


def add_after(pf, key_after, newkey, newval):

	before_dict = OrderedDict()
	after_dict = OrderedDict()
	new_dict = OrderedDict([(newkey, newval)])
	after = False 

	for key,val in pf.viewitems():

		if after:
			after_dict[key] = val
		else:
			before_dict[key] = val

		if key == key_after:
			after = True


	all_dicts = OrderedDict(before_dict.items() + new_dict.items()  + after_dict.items())

	return all_dicts

def vel_nu (nu, nu_0):
	'''
	Returns the velocity of a frequency nu with respect to 
	line centre nu_0 in units of cm/s
	'''
	delta_nu = nu - nu_0
	
	vel = C * delta_nu / nu_0
	
	return vel

def nu_vel(v, nu_0):
	'''
	Returns the frequency of a velocity v with respect to 
	line centre nu_0 in units of cm/s
	'''
	
	nu =  (nu_0 * v / C) + nu_0
	
	return nu


def normalize(w, f):

	a = get_continuum(w, f, [1215,1550,1400,2800], lmin=1300, lmax=3000, deg=2)

	#fnorm = f/make_f(a, w)

	#fnorm = f / util.get_flux_at_wavelength(w,f,1450)

	return make_f(a, w)



def make_c4_plot(s, name, angmin=70, angmax=82):

	nangles = len(s.colnames) - 9

	nu_c4 = C / (1550.0 * ANGSTROM)
	big_tick_labels(18)

	figure(figsize=(12,9))
	rcParams["axes.linewidth"]=1.5

	velocity_c4 = -C * ( (1550.0 - s["Lambda"]) / 1550.0)
	velocity_mg = -C * ( (2800.0 - s["Lambda"]) / 2800.0)
	velocity_al = -C * ( (1857.0 - s["Lambda"]) / 1857.0)

	lambda_fe = 2382
	velocity_fe = -C * ( (lambda_fe - s["Lambda"]) / lambda_fe)

	set_pretty()

	start = 3
	end = 3
	ntodo = 15

	big_tick_labels(14)

	angles = np.arange(72,80)
	iplot = 0
	colors = get_colors()

	subplots_adjust(hspace=0,wspace=0, left=0.08, right=0.97, top=0.97, bottom=0.08)


	for i in range(len(angles)):

		iplot += 1
		subplot(2,4,iplot)

		#title(s.colnames[i])
		angle_name = "A%iP0.50" % angles[i]

		fnorm = normalize(s["Lambda"], s[angle_name])
		#fnorm = s[s.colnames[i]]/fnorm

		#plot(s["Lambda"], fnorm)
		offset = 2
		plot(-velocity_c4/1e8, (3*offset)+s[angle_name]/fnorm, linewidth=3, label=r"C~\textsc{iv}~$1550$\AA")
		plot(-velocity_mg/1e8, (2*offset)+s[angle_name]/fnorm, linewidth=3, label=r"Mg~\textsc{ii}~$2800$\AA")
		plot(-velocity_al/1e8, (1*offset)+s[angle_name]/fnorm, linewidth=3, label=r"Al~\textsc{iii}~$1857$\AA")
		plot(-velocity_fe/1e8, 0+s[angle_name]/fnorm, linewidth=3, label=r"Fe~\textsc{ii}~$2382$\AA")


		bi = BALnicity(nu_c4, s["Freq."], s[angle_name])

		#xlim(1452,1638)
		xlim(30,-29)
		ylim(0,10.5)

		text(-18,8.3,r"$%i^\circ$" % int(angle_name[1:3]),fontsize=20)
		#text(-10,4.1,r"$BI=%i$" % bi,fontsize=16)



		print angle_name

		if iplot != 1 and iplot != 5 and iplot != 11:
			gca().set_yticklabels([])

		if iplot < 5:
			gca().set_xticklabels([])

		if iplot == 5:
			#xlabel("Velocity (1000 km/s)", fontsize=20)
			text(-60,-1.5,"Velocity (1000 km/s)", fontsize=24)

		if iplot == 5:
			#ylabel("Flux / Continuum", fontsize=20)
			text(43,13,"Flux / Continuum", fontsize=24, rotation="vertical")

		if iplot == 1: 
			handles, labels = gca().get_legend_handles_labels()
			for iii in range(4):
				text(29,7.3 - (iii*offset), labels[iii], color=colors[iii], fontsize=16)

	#subplots_adjust(hspace=0,wspace=0, left=0.08, right=0.97)
	#subplot(2,4,5)
	#text(-43,12,"Flux / Continuum", fontsize=24, rotation="vertical")
	#text(60,-1.3,"Velocity (1000 km/s)", fontsize=24)

	#subplot()

	savefig("c4_angles.png",dpi=300)
	savefig("../../figures/c4_angles.png",dpi=300)
	savefig("../../figures/c4_angles.eps",bbox="tight")


def BALnicity ( nu_line, nu_array, spec_array):

	'''
	calculate the BALnicity index of a line
	
	:INPUT
		nu_line		float
					The frequency of the line at line centre
		
		nu_array		float
					array of frequencies
		
		spec_array	float
					array of intensities/fluxes
		
	:OUTPUT
		BALnicity index		float
		
	'''
	
	# we integrate over positive velocities as looking for blueshifted BALs
	
	vlow = 3000.0*1.0e5		# need to give velocity in CGS units
	nu_low = nu_vel(vlow, nu_line)
	
	vhigh = 25000.0*1.0e5
	nu_high = nu_vel(vhigh, nu_line)
	
	ilow = 0; ihigh = -999
	
	for i_nu in range(len(nu_array)):
		 if nu_array[i_nu]>= nu_low and ilow==0:
		 	ilow = i_nu
		 	ihigh = 0
		 
		 if nu_array[i_nu]>= nu_high and ihigh==0:
		 	ihigh = i_nu
		 

	ww = C / nu_array / ANGSTROM

	#spec_array = normalise(spec_array, 0.05, normalise_point=682)
	#spec_array = spec_array / util.get_flux_at_wavelength(1450, ww, spec_array)
	
	constants = np.zeros(len(spec_array))
	
	for i in range(ilow, ihigh):
		if spec_array[i] <= 0.9: 
			constants[i] = 1
			
	#pylab.plot(nu_array[ilow:ihigh], spec_array[ilow:ihigh])
	
	#pylab.show()
	
	BI = 0.0	
	for i in range(ilow, ihigh):
		dv = ( vel_nu(nu_array[i], nu_line) - vel_nu(nu_array[i-1], nu_line) )
		contribution = constants[i] * ((1.0 - (spec_array[i] / 0.9)) * dv)

		BI += contribution
	
	BI = BI / 1.0e5
	
	return BI

# def make_bi_plot(s, angmin=70, angmax=82):

# 	nangles = len(s.colnames) - 9

# 	for i in range(nangles):


#name = "../specs/nextgen_uv"

#name = "../specs/checkvol_a6"
name = "../specs/grid10/nextgen_a05"
name = "../specs/webgrid/run6_thmin70_rmin50_a0p6_rv1e19_f0p01.spec"
name = "../specs/webgrid/run5_thmin70_rmin50_a0p5_rv1e19_f0p01.spec"
s = r.read_spectrum(name)
make_c4_plot(s, name)

# clf()

# figure(figsize=(18,6))
# d = np.loadtxt("data/atomic78/lines_linked_ver_2.py", unpack=True, usecols=(1,2,3,4))
# select = (d[0] == 26) * (d[1] == 2)
# fe2 = d[2][select]

# plot(s["Lambda"], util.smooth(s["A78P0.50"]), linewidth=2)
# vlines(fe2, 0.05,0.2 )
# xlabel("Wavelength \AA", fontsize=16)
# ylabel("Flux", fontsize=16)

# title("78 deg with Fe II lines marked")
# xlim(1000,3000)
# savefig("fe.png")

# clf()

# plot(s["Lambda"], util.smooth(s["A72P0.50"]), linewidth=2, label="72")
# plot(s["Lambda"], util.smooth(s["A75P0.50"]), linewidth=2, label="75")
# plot(s["Lambda"], util.smooth(s["A78P0.50"]), linewidth=2, label="78")

# float_legend()
# xlabel("Wavelength \AA", fontsize=16)
# ylabel("Flux", fontsize=16)
# semilogy()
# savefig("broadband.png")








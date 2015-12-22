#from plot_norm import *
import numpy as numpy
from pylab import *
import os, sys
from cobra_sub import smooth
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plot_norm import *
from constants import *
from pretty import *

#set_pretty()
xkcd()

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


def get_line_EW(s, angle_string = "A20P0.50", wline=6563, delta = 1.2):



	line_w, line_f = find_line(s, angle_string, w=wline)

	#print "Finding EW of line at %i A" % wline
	
	cont_t = line_f[1:10]
	cont = np.mean(cont_t)  

	delta = (np.mean(cont_t) + np.std(cont_t) ) / np.mean(cont_t)

	wlast = line_w[0]
	integral = 0.0
	sum_i = 0.0

	for i in range(1,len(line_f)):

		if line_f[i] > (delta * cont):

			dw = line_w[i] - wlast
			integral += dw * (1.0 - (line_f[i] / cont) )
			sum_i +=  line_f[i]
			wlast = line_w[i]

	#print "EW is %8.4f" % integral
	#print "I_cont times b is %8.4e" % (integral * cont) 
	#print "This should be similar to line flux %8.4e" % sum_i

	return integral

def get_line_RWEW(w,f, wline=1500, delta = 1.2):
	'''
	get the "red wing equivalent width" of an emission line with
	wavelenght wline 
	'''
	# first get the continuum fit
	f_continuum = normalize(w,f)

	# get an array of values around the line, rw =true means red wing only
	line_w, line_f = find_line(w, f/f_continuum, wline=wline, rw=True)

	# get the wavelength bins for the integration
	dw = np.fabs(line_w[1:] - line_w[:-1])
	dw = np.append(dw, dw[-1])

	# do the integral
	integrand = (1.0 - line_f) * dw
	integral = np.sum(integrand)

	return integral



def plot_ew():

	path = "/Users/jmatthews/Documents/runs/QSO_clumped/grid7_posterrorfixes/"
	#specs = [f for f in os.listdir(path) if f.endswith('.spec') and 'run' in f]
	specs = []
	specs.append("nextgen_a05_pre_he")

	print specs
	figure(figsize=(14,6))
	for ifile in range(len(specs)):

		if ifile != len(specs) - 1:
			fname = path + specs[ifile]
		else:
			fname = specs[ifile]
		s = r.read_spectrum(fname)

		ws = [1550,2800]
		labels=[r"$W_{\lambda, RW}$~(C\textsc{iv}, Model)",r"$W_{\lambda, RW}$~(Mg\textsc{ii}, Model)"]
		labels2 = [r"$1/2$~$W_{\lambda}$ (C\textsc{iv}, non-BAL, D12)", r"$1/2$~$W_{\lambda}$ (Mg\textsc{ii}, non-BAL, D12)", r"$1/2$~$W_{\lambda}$~(Mg\textsc{ii}, BAL, D12)"]
		sigmas = np.array([28.6, 30.9, 28.1])/2.0
		mus = np.array([45.1,47.2, 41.2])/2.0
		ymax = [150,100]
		w = s["Lambda"]
		
		big_tick_labels(18)
		long_ticks()


		for j in range(len(ws)):

			subplot(1,2,j+1)

			ew = np.zeros(20)
			angle = np.zeros(20)

			for i in range(len(s.colnames[9:])):
				angle_string = s.colnames[9 + i]
				f = s[angle_string]
				ew[i] = np.fabs(get_line_RWEW(w,f, wline=ws[j], delta = 1.2))
				angle[i] = angle_string[1:3]

			if ifile == len(specs) - 1:
				c = get_colors()
				alpha = 1
				color_plot = c[2]
				width_plot = 2
				label_data = labels[j] 
			else:
				color_plot = 'k'
				alpha = 0.2
				width_plot = 1
				label_data = None

			if ifile == len(specs) - 1:

				mu = mus[j]
				sigma = sigmas[j]
				y = np.zeros(len(np.arange(20,71))) + mu
				plot(np.arange(20,71), y, label=labels2[j], linewidth=2, c=c[1])
				fill_between(np.arange(20,71), y-sigma, y+sigma, facecolor=c[1], alpha=0.5, edgecolor="None")


			# plot the actual data
			plot(angle, ew, label=label_data, linewidth=width_plot, c=color_plot, alpha=alpha)


			if ifile == len(specs) - 1:
				if j == 1:
					mu = mus[j+1]
					sigma = sigmas[j+1]
					ybal = np.zeros(len(np.arange(70,83))) + mu
					plot(np.arange(70,83), ybal, label=labels2[j+1], linewidth=2, c=c[0])
					fill_between(np.arange(70,83), ybal-sigma, ybal+sigma, facecolor=c[0], alpha=0.5, edgecolor="None")
				

				theta0 = 20.0/180.0*np.pi
				theta = np.arange(20,90,0.1)/180.0*np.pi
				ew0 = ew[0]*np.cos(theta0)*(1.0 + (2.0/3.0)*np.cos(theta0))
				denom = np.cos(theta)*(1.0 + (2.0/3.0)*np.cos(theta))
				#ew0 = ew[0]*np.cos(theta0)
				#denom = np.cos(theta)
				plot(np.arange(20,90,0.1), (ew0/denom), c=c[3], label=r"$\cos i (1 + 2/3 \cos i)$ extrapolation", linewidth=2,linestyle="--")

				float_legend(loc=2)
				xlim(20,85)
				#semilogy()

				xlabel(r"Inclination ($^\circ$)", fontsize=20)

				if j == 0: ylabel(r"$W_{\lambda,RW}$ (\AA)", fontsize=20)

				#gca().set_xscale('log', basex=2)
				#gca().set_yscale('log', basey=2)
				vlines([70,82],2,ymax[j]*2,color="k",linewidth=2,linestyle="--")
				ylim(2,ymax[j])


		# tickl= [str(int(i+0.01)) for i in 2.0**np.arange(1,int(np.log2(ymax[j]))+1)]
		# gca().set_yticklabels([2,4,8,16,32,64])

		# labels = [int(item.get_text() for item in gca().get_yticklabels()]
		# # #labels[1] = 'Testing'
		# yticks(ew, labels)

		# gca().set_yticklabels(labels)

	subplots_adjust(wspace = 0.1)

	savefig("ew_xkcd.png", dpi=300, bbox_inches='tight')
	clf()

def normalize(w, f):

	a = get_continuum(w, f, [1215,1550,1400,2800], lmin=1300, lmax=3000, deg=2)

	#fnorm = f/make_f(a, w)

	#fnorm = f / util.get_flux_at_wavelength(w,f,1450)

	return make_f(a, w)

def plot_fit():

	s = r.read_spectrum("nextgen_a05_pre")

	w = s["Lambda"]

	figure(figsize=(8,6))

	for i in range(len(s.colnames[9:])):

		f = s[s.colnames[9+i]]

		a = get_continuum(w,f, lines = [1215,1400,1550,1640,2800,1900], lmin=1100, lmax = 2000, deg = 1)
	
		ff = normalize(w,f)
		f2000 = util.get_flux_at_wavelength(w,ff,2000)
		plot(w,f/ff,label="fit")
		#f2000 = util.get_flux_at_wavelength(w,f,2000)
		#plot(w,f, label="model")

		xlim(1500,1600)
		ylim(0,7)

		savefig("fits_%s.png" % s.colnames[9+i][1:3], dpi=300)
		clf()

plot_ew()


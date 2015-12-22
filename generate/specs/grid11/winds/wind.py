#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from pretty import *

set_pretty()

from matplotlib.colors import LinearSegmentedColormap


#cm_use = get_viridis()
cm_use=get_optb()
# f = open("nextgen_alpha1_long.complete", "r")

# data_old = []
# count = 0
# for line in f:
# 	data = line.split()

# 	if len(data) != len(data_old):
# 		print count, len(data), len(data_old)

# 	count+=1
# 	data_old = data

# #sys.exit()

RG = 8.85e14
x1000 = np.logspace(14.0, np.log10(RG*1000-10000000),num=1000)


y1000 = np.sqrt((RG*1000)**2 - x1000**2)

x1000 = np.append(x1000,RG*1000)
y1000 = np.append(y1000,1e-99)
#n
xuse= np.append(x1000,1e19)


print x1000, y1000


root = "asi_thmin70_rmin50_a0p75_f0p01_mdotw5"
d = r.read_pywind_summary(root)


strings = ["te", "ne", "ionC4", "ionAl3", "ionMg2", "lum_lines", "ly", "IP"]
levs = [np.arange(3.4,4.8,0.1), np.arange(7,11,0.1), np.arange(-4,0,0.1),np.arange(-4,0,0.1),np.arange(-4,0,0.1),np.arange(40,45,0.2), np.arange(37,42.5,0.2), np.arange(-3,2,0.2)]
labels = [r'$\rm{log (electron~temperature/K)}$', r'$\rm{log (electron~density/cm^{-3})}$',r'$\rm{log (C\textsc{iv}~ion~fraction)}$)', r'$\rm{log (Al\textsc{iii}~ion~fraction)}$)', r'$\rm{log (Mg\textsc{ii}~ion~fraction)}$)', r'$\rm{log (total~line~luminosity/erg~s^{-1})}$)',
r'$\rm{log (Ly}~\rm{\alpha~line~luminosity)/erg~s^{-1}}$)', r'$\log U$']
log = [1,1,1,1,1,1,1,1]


fig = figure(figsize=[12,15])

for i in range(len(strings)):

	subplot(4,2,i+1)

	if strings[i] == "ionMg2" or strings[i] == "ionAl3":
		x,z,var = r.read_pywind(root+".%s.dat" % strings[i])

	elif strings[i] != "ly":
		x,z,var = util.wind_to_masked(d, strings[i])
	else:
		x,z,var = r.read_pywind(root+".lev2_emiss.dat")
		#x,z,var = util.wind_to_masked(d, "var")

	if log[i]:
		var = np.log10(var)

	contourf(x,z,var, extend='both', levels=levs[i], cmap=cm_use)
	grid()
	loglog()
	xlim(3e16,1e19)
	ylim(1e15,1e19)
	ylabel(r"$z$ (cm)")
	xlabel(r"$x$ (cm)")
	ax = gca()

	cbar = colorbar(orientation='horizontal')
	cbar.ax.set_xlabel(labels[i])

	plot(x1000, y1000,c="k")

	angles=[75,80]
	


	for angle in angles:
		zcoord=[]
		for xcoord in xuse:
			zcoord.append(xcoord*np.tan(np.radians(90.-angle)))
		plt.plot(xuse,zcoord,color='k',linewidth=0.75, linestyle="--")
		plt.text(10.0**16.55,(10**16.6)*np.tan(np.radians(90.-angle))-((angle - 76.2)*1e15),angle,fontsize=12)



subplots_adjust(bottom=0.05, top=0.95)

savefig("wind.png", dpi=300)
savefig("../../figures/wind.eps", bbox="tight")








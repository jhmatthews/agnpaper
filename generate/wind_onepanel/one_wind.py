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


cm_use = get_viridis()
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


d = r.read_pywind_summary("nextgen_a05")


strings = ["te"]
levs = [np.arange(-3,2,0.2)]
levs = [np.arange(3,5,0.1)]
labels = [r'$\log U$']
log = [1]


fig = figure(figsize=[1,1])

x,z,var = util.wind_to_masked(d, strings[0])

if log[0]:
	var = np.log10(var)



# select = (var.mask == False)*(z<1e19)

# z = z[select]
# x = x[select]
# # e = e[(e.mask == False)]

# import matplotlib.tri as tri
# T = tri.Triangulation(x.flatten(), z.flatten())


# var = var[select]
# var = var.flatten()


#tricontourf(T,np.log10(var), mask = var.mask, levels = levs[0], extend="both", cmap=cm_use)

contourf(x,z,var, extend='both', levels=levs[0], cmap=cm_use)
#grid()
#loglog()
pretty_axes()
xlim(-8.85e14,1e18)
ylim(-8.85e14,1e18)

# text(5e17,6e17,"`Quasar-like', $20-70^\circ$", fontsize=16)

#text(5e17,8e17,"BLR", fontsize=16)

# text(5e17,4e17,"`LoBALQSO'", fontsize=16)

# #text(5e17,2e17,"`Below Wind', $82^\circ+$", fontsize=16)

# text(5e17,2e17,"`Below Wind', $82^\circ+$", fontsize=16)
#ylabel(r"$z$ (cm)")
#xlabel(r"$x$ (cm)")
ax = gca()

#cbar = colorbar(orientation='horizontal')
#cbar.ax.set_xlabel(labels[0])

axis('off')
#plot(x1000, y1000,c="k")

angles=[75,80]

gca().get_xaxis().set_visible(False)
gca().get_yaxis().set_visible(False)

for angle in angles:
	zcoord=[]
	for xcoord in xuse:
		zcoord.append(xcoord*np.tan(np.radians(90.-angle)))
	#plt.plot(xuse,zcoord,color='k',linewidth=0.75, linestyle="--")
	#plt.text(10.0**16.55,(10**16.6)*np.tan(np.radians(90.-angle))-((angle - 76.2)*1e15),angle,fontsize=12)



# xdisk = np.arange(0,1e19,1e17)
# ydisk = np.zeros(len(xdisk))
# #plot(xdisk,ydisk,linewidth=4, c="k")



# thetamin=72.0
# thetamax=80.0
# diskmax=25.0
# xmin=10.0
# xmax=18.0
# rstar=2.0

# xwindmax=40.0
# zwindmax=(xwindmax-xmin)*np.tan(np.radians(90.0-thetamin))


# posx=np.arange(xmin,xwindmax,(xmax-xmin)/10.0)
# ymin1=[]
# ymax1=[]

# x=np.arange(-8.85e16,8.85e14,1e12)
# y= np.sqrt((8.85e14**2) - (x**2))

# plot(x,y)





# plot(x,y1,x,y2,color='black')
# plot(-1.0*x,y1,-1.0*x,y2,color='black')
# plot(-1.0*x,-1.0*y1,-1.0*x,-1.0*y2,color='black')
# plot(x,-1.0*y1,x,-1.0*y2,color='black')
# #set_xlim([-1.0*xwindmax,xwindmax])
# fill_between(x,y1, y2, facecolor='0.9', interpolate=True)
# fill_between(-1.0*x,y1, y2, facecolor='0.9', interpolate=True)
# fill_between(x,-1.0*y1, -1.0*y2, facecolor='0.9', interpolate=True)
# fill_between(-1.0*x,-1.0*y1, -1.0*y2, facecolor='0.9', interpolate=True)
# plot([-1.0*diskmax,diskmax],[0.0,0.0],color='k',linewidth=2)
# plot(0.0,0.0,'o',color='w',ms=20.0)


savefig("onepanel_wind.png", dpi=500, transparent=True)
#savefig("onepanel_wind.png", dpi=500)
#savefig("../../figures/onepanel_wind.eps", bbox="tight")








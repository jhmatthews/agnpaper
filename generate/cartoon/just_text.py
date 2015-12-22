'''
SV_cartoon.py 

This script creates the cartoon figure for 
the 2014 CV paper. The data were scraped from Dexter.
'''

import matplotlib.mlab as mlab
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#plt.xkcd()


#mpl.rcParams['pdf.use14corefonts'] = True
#mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#mpl.rc('text', usetex=True)

thetamin=72.0
thetamax=80.0
diskmax=25.0
xmin=10.0
xmax=18.0
rstar=2.0

xwindmax=40.0
zwindmax=(xwindmax-xmin)*np.tan(np.radians(90.0-thetamin))


posx=np.arange(xmin,xwindmax,(xmax-xmin)/10.0)
ymin1=[]
ymax1=[]

x=np.array([xmin,xmax,xwindmax])
y1=np.array([0.0,(xmax-xmin)*np.tan(np.radians(90.0-thetamin)),(xwindmax-xmin)*np.tan(np.radians(90.0-thetamin))])
y2=np.array([0.0,0.0,(xwindmax-xmax)*np.tan(np.radians(90.0-thetamax))])

lsize=14
tsize=16 

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.plot(x,y1,x,y2,color='black', linewidth=0)
ax.plot(-1.0*x,y1,-1.0*x,y2,color='black', linewidth=0)
ax.plot(-1.0*x,-1.0*y1,-1.0*x,-1.0*y2,color='black', linewidth=0)
ax.plot(x,-1.0*y1,x,-1.0*y2,color='black', linewidth=0)
ax.set_xlim([-1.0*xwindmax,xwindmax])

# plt.fill_between(x,y1, y2, facecolor='0.9', interpolate=True)
# plt.fill_between(-1.0*x,y1, y2, facecolor='0.9', interpolate=True)
# plt.fill_between(x,-1.0*y1, -1.0*y2, facecolor='0.9', interpolate=True)
# plt.fill_between(-1.0*x,-1.0*y1, -1.0*y2, facecolor='0.9', interpolate=True)
# plt.plot(0.0,0.0,'o',color='w',ms=20.0)
# plt.plot(0.0,0.0,'o',color='k',ms=5.0)

#ax.text(-20,9,r'$\rm{Model~elements~producing~photons}$',fontsize=tsize)
#ax.text(-22,-9,r'$\rm{Geometric~parameters~defining~model}$',fontsize=tsize)


''' This is for the radiating sources'''
# ax.annotate(r'$\rm{Central~Source}$', xy=(-0.5,0.6),xytext=(-17,8),fontsize=lsize,arrowprops=dict(arrowstyle="->"),)
# ax.annotate(r'$\rm{Accretion~Disc}$', xy=(-5.0,0.1),xytext=(-26,5.5),fontsize=lsize,arrowprops=dict(arrowstyle="->"),)
# plt.text(-35,5,r'$\rm{Biconical~wind}$',fontsize=lsize,rotation=-37.0)

'''This is for the wind parameters'''
#ax.plot([x[0],x[0]],[0,-4],'--',color='k')
#ax.plot([x[1],x[1]],[0,-5],'--',color='k')
ax.text(4.9,-4,r'$\rm{r_{min}}$',fontsize=lsize)
ax.text(12.8,-5,r'$\rm{r_{max}}$',fontsize=lsize)
#ax.annotate('',xy=(14.3,-1.4), xycoords='data',xytext=(10,-2),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3 ,rad=0.3"))
ax.text(11.1,-2.8,r'$\rm{\theta_{min}}$',fontsize=lsize)
#ax.annotate('',xy=(24.3,-1.0), xycoords='data',xytext=(18,-2),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3 ,rad=0.3"))
ax.text(21,-2.7,r'$\rm{\theta_{max}}$',fontsize=lsize)
plt.text(27,-3,r'$\rm{\dot{M}_{wind}}$',fontsize=lsize,rotation=-37.0)
#ax.annotate('',xy=(35,-5.0), xycoords='data',xytext=(27,-3),textcoords='data',arrowprops=dict(arrowstyle="->"))

'''This is for the disk parameters'''
#ax.plot([-1.9,-1.9],[0,-4],'--',color='k')
#ax.plot([-25,-25],[0,-6],'--',color='k')
ax.text(-6.3,-4,r'$\rm{r_{X}}$',fontsize=lsize)
ax.text(-24,-6,r'$\rm{r_{disc}(max)}$',fontsize=lsize)
#ax.annotate('',xy=(-26,0), xycoords='data',xytext=(-35,0),textcoords='data',arrowprops=dict(arrowstyle="->"))
plt.text(-33,-1,r'$\rm{\dot{M}_{acc}}$',fontsize=lsize)

ax.set_xlim(-40,40)
plt.axis('off')
plt.savefig('text.png',bbox_inches='tight', dpi=300, transparent=True)
#plt.savefig('../figures/fig2_cartoon.eps',bbox_inches='tight')




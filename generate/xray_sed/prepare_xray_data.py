#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
fname = "grid11/nextgen_a05_pre_xray"
#fname = "fiducial_xray_fine"
angle = 75

s = r.read_spectrum("../specs/"+fname)

norm = (100.0 * PARSEC) * (100.0 * PARSEC) / (1e10 * PARSEC) / (1e10 * PARSEC)

fl_to_fnu = s["Lambda"] / s["Freq."]

flambda = s["A%iP0.50" % angle]

fnu = flambda * fl_to_fnu * norm

data_to_write = np.ndarray((2, len(fnu)))

data_to_write[0] = s["Freq."]
data_to_write[1] = fnu

data_to_write = np.transpose(data_to_write )

header_string = "Frequency(Hz)  |   F_nu (erg/s/cm^2/Hz)\nPower law has spectral index 0.9"
np.savetxt('qsomod_a05_xray_75deg.spec', data_to_write, header=header_string)
#np.savetxt('h13mod_xray_75deg.spec', data_to_write, header=header_string)

plot(s["Freq."], fnu)
savefig("testspec.png")
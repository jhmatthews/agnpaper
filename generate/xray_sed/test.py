#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np


sed = np.loadtxt("pg1001_sed_rest.dat",unpack=True)
yerr = [sed[2],sed[3]]

nu = sed[0]
nulnu = sed[1]

yerr[0] = -10.0 ** (nulnu + yerr[0]) + (10.0**nulnu)
yerr[1] = 10.0 ** (nulnu + yerr[1]) - (10.0**nulnu)

errorbar(10.0**sed[0], 10.0**sed[1], yerr=yerr, fmt=".")

print yerr

loglog()
savefig("test.png")
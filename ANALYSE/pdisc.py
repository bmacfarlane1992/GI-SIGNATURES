#
# pdisc.py
#
# Python program to read in pdisc (radially varied, temporally static) data structures.
# Uses information of file name to generate plots analogous to Stamatellos,
# Whitworth & Hubber (2012) work.
#
# Author: Benjamin MacFarlane
# Date: 29/07/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#

#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
from scipy.optimize import curve_fit
from scipy import interpolate
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def read(arch_dir, dat_dir, snaparr):
#
    print("pdisc files being read")
#
    # Define number of files to be read dependent on snaparr dimensions
#
    if (snaparr.ndim == 1):
        file_n = len(snaparr)
    else:
        file_n = len(snaparr)*len(snaparr[0])
#
    # Arrays to fill
#
    r = [[] for i in range(file_n)] ; vkep = [[] for i in range(file_n)]
#
    # Create string array as pointer to pdisc files based on snaparr
#
    file_list = []
    fcount = 0
    if (snaparr.ndim ==1):
        snaparr_tmp = np.array([0]*len(snaparr))
        for a in range(0, len(snaparr)):
                if (snaparr[a] < (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                       str(snaparr[a]+69)+'.pdisc.1')
                elif (snaparr[a] > (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                    str(snaparr[a]+69)+'.pdisc.1')
                    snaparr_tmp[a] = fcount
                    fcount = fcount + 1
    else:
        snaparr_tmp = np.array([[0]*len(snaparr[0])]*len(snaparr))
        for a in range(0, len(snaparr)):
            for b in range(0, len(snaparr[0])):
                if (snaparr[a][b] < (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                       str(snaparr[a][b]+69)+'.pdisc.1')
                elif (snaparr[a][b] > (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                    str(snaparr[a][b]+69)+'.pdisc.1')
                    snaparr_tmp[a][b] = fcount
                    fcount = fcount + 1
#
    # Loop over files in pdisc directory to read all temporally evolved values
#
    for i in range(0, file_n):
#
    # Now loop over lines in each file and extract variables of interest
    # then convert to numpy format for future data extraction for plotting
#
        f = open(file_list[i], 'r')
#
    # Firstly, read and ignore file header
#
        header = f.readline()
#
        for line in f:
            line = line.strip() ; columns = line.split()
            r[i].append(float(columns[1])) ; vkep[i].append(float(columns[13]))
        f.close()
    r = np.array(r) ; vkep = np.array(vkep)
#
#
    return r, vkep

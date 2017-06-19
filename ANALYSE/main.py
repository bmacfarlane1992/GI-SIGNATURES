#
# main.py
#
# Python program to execute read and analysis of simulation files.
# See script headers for more details on code I/O.
#
# Author: Benjamin MacFarlane
# Date: 08/11/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
arch_dir = "/Users/bmacfarlane1992/Documents/ACADEMIA_PHD/PROJECTS/GI-SIGNATURES/"    # Location of project directory architecture
dat_dir = "/Volumes/SD_bm1992/DATA/"            # Location of data from analyse_disc.f90
ic_dir = "/Volumes/SD_bm1992/ICs/"              # Location of ICs to read in sink data
plot_form = "png"								# Format of plots (either "eps" or "png")
col_arr = ["b", "g", "r", "c", "m", "k"]        # Array of color pointers for plots throughout
#
v_K = ["90"]            # Keplerian velocity percentage restriction on disc mass/radius
inclin = ["60","90"]    # Inclination of disc being analysed
#
ea_run = [0,1,3,4,5]          # Select EA runs to process
part_mass = 5.e-6            # Mass of SPH particles
#
r_limit = 150           # Limit of radial plots in pdisc analyses
#
pv = "TRUE"             # Choose whether ("TRUE") or not ("FALSE") to generate PV diagram
mcomp_tmp = "TRUE"      # Choose whether ("TRUE") or not ("FALSE") to generate mass comparison
#                       # of simulation vs. PV diagram analysis system masses
#
EA_timeref = ["BEFORE","POST-ONSET","PRE-OFFSET","AFTER"]       # Define names of EA snapshots for EA length reference
EA_lenref = ["SHORT","MEDIUM","LONG"]                           # and time reference to EA outburst event
#
d = {}
snap0 = 'snap0' ; snap0_IA = 'snap0_IA'
snap1 = 'snap1' ; snap1_IA = 'snap1_IA'
snap3 = 'snap3' ; snap3_IA = 'snap3_IA'     # snaparr formatted as [before,during,after]  for [short,medium,long]
snap4 = 'snap4' ; snap4_IA = 'snap4_IA'     # accretion events. Uses dictionary to ensure that snaparr can be
snap5 = 'snap5'    ; snap5_IA = 'snap5_IA'  # manipulated into snaparr{i,j} array. DE05 incdices also listed in ###
#
d[snap0] = [410, 610, 800] ; d[snap0_IA] = []
    ### [479, 529, 579, 629, 679, 729, 779, 829, 849, 869] ###
d[snap1] = [300, 1050, 1800] ; d[snap1_IA] = []
    ### [219, 369, 519, 669, 819, 969, 1119, 1269, 1419, 1569, 1719, 1869] ###
d[snap3] = [[100,150,200,250],[400,450,550,600],[1050,1150,1350,1450]]
    ### [[00169,00219,00269,00319],[00469,00519,00619,00669],[01119,1219,01419,01519]] ###
d[snap4] = [[200,209,212,220],[860,875,880,900]]
    ### [[00269,00278,00281,00289],[00929,00944,00949,00969]] ###
d[snap5] = [[150,183,191,224],[719,755,765,801]]
    ### [[00219,00252,00260,00293],[00788,00824,00834,00870]] ###
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
pcAU = 206265.        # Conversion from parsec -> AU
pcm = 3.086e16        # Conversion from parsec -> m
AUm = 1.496e11        # Conversion from AU -> m
cm2m2 = 1e-4          # Conversion for 1/m**2 to 1/cm**2
Jerg = 1e7              # Conversion from Joules to erg
cgsSIopac = 1e-1       # Conversion of cgs opacity to SI units
G = 6.67e-11        # Gravitational constant (SI)
Msol_kg = 1.998e30    # Conversion from Solar masses -> kg
kb = 1.38e-23         # Bolztmann constant (SI)
c = 2.998e8           # Spped of light (SI)
h = 6.626e-34         # Planck constant (SI)
mH = 1.67e-27         # Mass of Hydrogen atom (SI)
mu = 2.3             # Mean molecular weight for PPD
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import os            # Standard Python modules
import numpy as np
#
import r1            # Local modules
import sink
import pdisc
import pv_diag
import mass_comp
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    # Loop over run number
#
for i in range(0, len(ea_run)):
#
# Dummy arrays to which PV masses are entered. Prevents overwriting of data in
# pv_diag.py script
#
    pv_MASK60 = [] ; pv_MASK90 = []
#
    # Loop over ea inclination values
#
    for j in range(0, len(inclin)):
#
    # Loop over Keplerian restrictions
#
        for k in range(0, len(v_K)):
#
            print "\n Run "+str(ea_run[i])+ " from " \
               "vK"+str(v_K[k])+"_"+str(inclin[j])+"i/ directory"+ \
               " is now being analysed \n"
#
            plotdir = arch_dir+"PLOTS/"+str(ea_run[i])
            os.system('mkdir '+str(plotdir))
            plotdir = plotdir+"/"+"vK"+v_K[k]+"_"+inclin[j]+"i/"
            os.system('mkdir '+str(plotdir))
#
            dat_dir_tmp = dat_dir+str(ea_run[i])+"/"+"vK"+v_K[k]+"_"+inclin[j]+"i/"
#
            snaparr = d['snap'+str(ea_run[i])] ; snaparr = np.array(snaparr)
#
    # Read in rdisc.1 file
#
            hasharr_app, n_accr, r_d_kep, r_d_sigALMA, m_s, m_mri_d, m_d_kep, \
            m_d_piv, m_d_sigALMA, = r1.read(dat_dir_tmp, plotdir, ea_run[i], \
            snaparr, v_K[k], inclin[j], plot_form)
#
#
    # DE05.sink file(s) read (Planet sink data vs. radius vs. time)
#
            pmass, pradius, timearr = sink.read(arch_dir, ic_dir, plotdir, \
               ea_run[i], snaparr, pcAU)
#
    # pdisc read (Disc parameters vs. radius for defined time)
#
            r, vkep = pdisc.read(arch_dir, dat_dir_tmp, snaparr)
#
    # Position-Velocity diagram plot
#
            if ((pv == "TRUE")):
                pv_mass60, pv_mass90, kep_fit, raw_fit = pv_diag.pv(dat_dir_tmp, \
                plotdir, ea_run[i], snaparr, v_K[k], inclin[j], r, vkep, pmass, \
                EA_lenref, EA_timeref, pcAU, AUm, G, Msol_kg, pv_MASK60, \
                pv_MASK90, pcm, mH, plot_form, kb, c, h, cgsSIopac, mu, Jerg, \
                cm2m2, part_mass )
#
                if (inclin == "60"):
                    pv_MASK60 = pv_mass60
#
    # Mass comparison of simulation to PV data (both 60 and 90 deg.)
#
    if ( (raw_fit == "FALSE") and (kep_fit == "TRUE") and (mcomp_tmp == "TRUE")):
        mass_comp.comp(dat_dir_tmp, plotdir, ea_run[i], snaparr, \
        timearr, v_K, inclin, m_s, m_mri_d, pmass, m_d_kep, \
        m_d_piv, m_d_sigALMA, pv_mass60, pv_mass90, EA_lenref, plot_form)
#
    # Convert all images to .eps format, suitable for publication
#
#    os.system('for f in '+plotdir+'*.pdf; do pdftops -eps $f; done')
#    os.system('rm -r '+plotdir+'*.pdf')
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Remove garbage from analysis directory and end program
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
os.system("rm -r *.pyc")
os.system("rm -r *~")
exit()

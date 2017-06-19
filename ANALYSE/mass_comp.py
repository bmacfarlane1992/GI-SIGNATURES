#
# mass_comp.py
#
# Programme to plot comparison of masses in simulation and PV diagram analyses
#
# Author: Benjamin MacFarlane
# Date: 04/11/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
print_term1 = "FALSE"    # Choose whether ("TRUE") or not ("FALSE") to output mass comparisons to terminal
print_term2 = "TRUE"
EWASS = "FALSE"            # Choose whether ("TRUE") or not ("FALSE") to compare inclination dependant PV masses
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def comp(dat_dir, plotdir, ea_run, snaparr, timearr, v_K, inclin, m_s, m_mri_d, \
   pmass, m_d_kep, m_d_piv, m_d_sigALMA, pv_mass60, pv_mass90, EA_lenref, \
   plot_form):
#
    print "Mass comparisons of simulation and PV data now being plotted"
#
    if ea_run == 0:
        paper_tag = 'A'
    elif ea_run == 1:
        paper_tag = 'B'
    elif ea_run == 3:
        paper_tag = 'C'
    elif ea_run == 4:
        paper_tag = 'D'
    elif ea_run == 5:
        paper_tag = 'E'
#
#
    # For EA [2, 3, 4, 5, 6] runs
#
    if (snaparr.ndim == 1):
        file_n = len(snaparr)
    else:
        file_n = len(snaparr)*len(snaparr[0])
    snaparr_tmp = [0]*file_n ; timearr_tmp = [0]*file_n ; pmass_tmp = [0]*file_n
#
    fcount = 0
    if (snaparr.ndim == 1):
        for i in range(0, len(snaparr)):
            snaparr_tmp[fcount] = snaparr[i]
            timearr_tmp[fcount] = timearr[i]
            fcount = fcount + 1
    else:
        for i in range(0, len(snaparr)):
            for j in range(0, len(snaparr[0])):
                snaparr_tmp[fcount] = snaparr[i][j]
                timearr_tmp[fcount] = timearr[i][j]
                fcount = fcount + 1
#
#                for a in range(1, len(pmass)-1):
#                    pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i][j]
#
#
    # Read in accretion parameters for runs with accretion
#
    t_s = [] ; t_e = []
    f = open(dat_dir+'../acc_params.dat','r')
    for line in f:
        line = line.strip() ; columns = line.split()
        t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
    f.close()
    t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Loop over snapshots, and compute masses for comparison
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    m_sys_kep = [0] * file_n ; m_sys_piv = [0] * file_n ; m_sys_sigALMA = [0] * file_n
    m_proto = [0] * file_n
#
    for i in range(0, file_n):
        m_sys_kep[i] = m_s[i] +  m_mri_d[i]  + pmass_tmp[i] + m_d_kep[i]
        m_sys_piv[i] = m_s[i] +  m_mri_d[i] + pmass_tmp[i] + m_d_piv[i]
        m_sys_sigALMA[i] = m_s[i] + m_mri_d[i] + pmass_tmp[i] + m_d_sigALMA[i]
        m_proto[i] = m_s[i] + m_mri_d[i]
#
        if (print_term1 == "TRUE"):
            print "for snaparr value: "+str(snaparr_tmp[i])
            print "\tMass of central protostar: "+str(m_s[i]+m_mri_d[i])
            print "\tMass of planets in disc: "+str(pmass_tmp[i])
            print "\tMass of star+disc+planets (Keplerian velocity criterion): "+str(m_sys_kep[i])
            print "\tMass of star+disc+planets (Oribit infall criterion): "+str(m_sys_piv[i])
            print "\tMass of star+disc+planets (ALMA density criterion): "+str(m_sys_sigALMA[i])
            print "\tMass of system (Keplerian fitted P-V diagram): "+str(pv_mass[i])
#
    if (print_term2 == "TRUE"):
        print "Mass % discrepancy between Sigma criterion and PV mass:"
        for i in range(0, file_n):
            print round(pv_mass90[i]/m_sys_sigALMA[i], 3)
#        print "Mass % discrepancy between 90 deg. and 60 deg. inclination PV  mass fits:"
#        for i in range(0, file_n):
#            print pv_mass60[i]/pv_mass90[i]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Plotting of mass comparisons, for different accretion events
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    f = plt.figure(1)
    f.set_rasterized(True)
#
    ax1 = plt.subplot(111)
    ax1.set_rasterized(True)

#
    mass_stack60 = pv_mass60 + m_sys_kep + m_sys_piv + m_sys_sigALMA
    mass_stack90 = pv_mass90 + m_sys_kep + m_sys_piv + m_sys_sigALMA
#
    for i in range(0,len(t_e)):
        plt.fill_between([t_s[i],t_e[i]], 0, 1.0, color='k', alpha = 0.5)
    plt.scatter(timearr_tmp, pv_mass90, \
       label = 'PV: 90 deg.', s=80, facecolors='none', edgecolors='b')
    plt.scatter(timearr_tmp, pv_mass60, \
       label = 'PV: 60 deg.', s=80, facecolors='b', edgecolors='b')
    plt.scatter(timearr_tmp, m_sys_kep, \
       label = 'Keplerian fraction cutoff', marker = 's', s=80, facecolors='none', edgecolors='g')
    plt.scatter(timearr_tmp, m_sys_piv, \
       label = 'Radial infall', marker = '+', s=80, facecolors='none', edgecolors='r')
    plt.scatter(timearr_tmp, m_sys_sigALMA, \
       label = 'Surface Density cutoff', marker = '^', s=80, facecolors='none', edgecolors='k')
    plt.ylim(0, ax1.get_ylim()[1])
    plt.scatter(timearr_tmp, m_proto, label = 'Protostar mass', marker = '*', \
    facecolors='k', edgecolors='k', s = 80)
    if (ea_run == 0):
        legend = plt.legend(loc = 'upper left', fontsize=16, scatterpoints=1)
    if ((snaparr.ndim == 2) and (ea_run == 3)):
        plt.xlim(timearr_tmp[4]-0.25, timearr_tmp[11]+0.25)
    elif ((snaparr.ndim == 2) and (ea_run > 3)):
        plt.xlim(timearr_tmp[4]-0.25, timearr_tmp[7]+0.25)
    elif ((snaparr.ndim == 1)):
        plt.xlim(timearr_tmp[0]-0.25, timearr_tmp[2]+0.25)
    plt.ylim(0, 1.0)
    plt.xlabel('Time (kyr)',fontsize=16)
    plt.ylabel('Mass '+(r'(M$_{\odot}$)'),fontsize=16)
#
    plt.savefig(plotdir+'mass_Time_'+paper_tag+'.'+str(plot_form), \
       transparant=True, rasterized=True, format=str(plot_form), dpi=150)
    plt.clf()

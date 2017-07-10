#
# mass_comp.py
#
# Programme to plot comparison of masses in simulation and PV diagram analyses
#
# Author: Benjamin MacFarlane
# Date: 21/06/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
print_term = "FALSE"	# Choose whether ("TRUE") or not ("FALSE") to output mass comparisons to terminal
EWASS = "FALSE"			# Choose whether ("TRUE") or not ("FALSE") to compare inclination dependant PV masses
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import numpy as np
import matplotlib.pyplot as plt
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def comp(dat_dir, plotdir, ea_run, snaparr, snapcore, timearr, v_K, inclin, m_s, m_mri_d, \
   pmass, m_d_kep, m_d_piv, m_d_sigALMA, pv_mass, EA_lenref):
#
	print "Mass comparisons of simulation and PV data now being plotted"
#
	# Read in masses independent of snaparr shape
#
	# For EA [0, 1] runs
#
	if (snaparr.ndim == 1):
		file_n = len(snapcore)
		snaparr_tmp = [0]*file_n ; timearr_tmp = [0]*file_n ; pmass_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snapcore)):
				if (snaparr[i] == snapcore[j]):
					snaparr_tmp[fcount] = snaparr[i] ; timearr_tmp[fcount] = timearr[i]
					for a in range(1, len(pmass)):
						pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i]
					fcount = fcount + 1
		t_s = [] ; t_e = [] ; n_accr = []
#
	# For EA [2, 3, 4, 5, 6] runs
#
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
		snaparr_tmp = [0]*file_n ; timearr_tmp = [0]*file_n ; pmass_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snaparr[0])):
				snaparr_tmp[fcount] = snaparr[i][j] ; timearr_tmp[fcount] = timearr[i][j]
#
#				for a in range(1, len(pmass)-1):
#					pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i][j]
#				
				fcount = fcount + 1
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
#
	for i in range(0, file_n):
#
		m_sys_kep[i] = m_s[i] + m_mri_d[i] + pmass_tmp[i] + m_d_kep[i]
		m_sys_piv[i] = m_s[i] + m_mri_d[i] + pmass_tmp[i] + m_d_piv[i]
		m_sys_sigALMA[i] = m_s[i] + m_mri_d[i] + pmass_tmp[i] + m_d_sigALMA[i]
#
		if (print_term == "TRUE"):
			print "for snaparr value: "+str(snaparr_tmp[i])
			print "\tMass of central protostar: "+str(m_s[i]+m_mri_d[i])
			print "\tMass of planets in disc: "+str(pmass_tmp[i])
			print "\tMass of star+disc+planets (Keplerian velocity criterion): "+str(m_sys_kep[i])
			print "\tMass of star+disc+planets (Oribit infall criterion): "+str(m_sys_piv[i])
			print "\tMass of star+disc+planets (ALMA density criterion): "+str(m_sys_sigALMA[i])
			print "\tMass of system (Keplerian fitted P-V diagram): "+str(pv_mass[i])
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of mass comparisons, for different accretion events
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	plt.figure(1)
#
	ax1 = plt.subplot(111)
#
	mass_stack = pv_mass + m_sys_kep + m_sys_piv + m_sys_sigALMA ; max_mass = max(mass_stack)+0.1
#
	for i in range(0,len(t_e)):
		plt.fill_between([t_s[i],t_e[i]], 0, max_mass, color='k', alpha = 0.5)
	plt.scatter(timearr_tmp, pv_mass, \
	   label = 'PV', s=80, facecolors='none', edgecolors='b')
	plt.scatter(timearr_tmp, m_sys_kep, \
	   label = 'Keplerian fraction cutoff', marker = 's', s=80, facecolors='none', edgecolors='g')
	plt.scatter(timearr_tmp, m_sys_piv, \
	   label = 'Radial infall', marker = '+', s=80, facecolors='none', edgecolors='r')
	plt.scatter(timearr_tmp, m_sys_sigALMA, \
	   label = 'Surface Density cutoff', marker = '^', s=80, facecolors='none', edgecolors='k')
	plt.ylim(0, ax1.get_ylim()[1])
	legend = plt.legend(loc = 'upper left', fontsize=8, scatterpoints=1)
	if ((snaparr.ndim == 2) and (ea_run == 3)):
		plt.xlim(timearr_tmp[4]-0.25, timearr_tmp[11]+0.25)
	elif ((snaparr.ndim == 2) and (ea_run > 3)):
		plt.xlim(timearr_tmp[4]-0.25, timearr_tmp[7]+0.25)
	plt.ylim(0, max_mass)
	plt.xlabel('Time (kyr)' ) ; plt.ylabel('Mass '+(r'(M$_{\odot}$)') )
#
	plt.savefig(plotdir+'mass_comp.pdf') ; plt.clf()
#
	if (EWASS == "TRUE"):
#
		pv_mass60 = [0.1488, 0.1715, 0.1693, 0.2212, 0.3719, 0.4063, 0.3695, 0.3425]
		
		plt.figure(1)
		ax1 = plt.subplot(111)
#
		mass_stack = pv_mass + m_sys_kep + m_sys_piv + m_sys_sigALMA ; max_mass = max(mass_stack)+0.1
#
		for i in range(0,len(t_e)):
			plt.fill_between([t_s[i],t_e[i]], 0, max_mass, color='k', alpha = 0.5)
		plt.scatter(timearr_tmp, pv_mass, label = 'PV: Edge-on', s=80, facecolors='none', edgecolors='b')
		plt.scatter(timearr_tmp, pv_mass60, label = 'PV: 60 deg.', s=80, facecolors='b', edgecolors='b')
		plt.scatter(timearr_tmp, m_sys_kep, label = 'Keplerian fraction cutoff', marker = 's', \
		   s=80, facecolors='none', edgecolors='g')
		plt.scatter(timearr_tmp, m_sys_sigALMA, label = 'Surface Density cutoff', marker = '^', \
		   s=80, facecolors='none', edgecolors='k')
		plt.ylim(0, ax1.get_ylim()[1])
		legend = plt.legend(loc = 'upper right', fontsize=8, scatterpoints=1)
		plt.ylim(0.2, max_mass) ; plt.xlim(85.5, 87.)
		plt.xlabel('Time (kyr)' ) ; plt.ylabel((r'M$_{sys}$')+' '+(r'(M$_{\odot}$)') )
#
		plt.savefig(plotdir+'EWASS_S9.pdf') ; plt.clf()
#

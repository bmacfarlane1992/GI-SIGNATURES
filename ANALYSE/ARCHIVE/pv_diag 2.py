#
# pv_diag.py
#
# Programme designed to read in simulation data and generate synthetic PV diagram 
#
# Author: Benjamin MacFarlane
# Date: 08/06/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
sig_mult = 0		# Threshold sigma value over which edge of Keplerian profile is fitted to
#
kep_fit = "FALSE"	# Choose whether ("TRUE") or not ("FALSE") to fit Keplerian profiles to PV data
#
fit_start = 0		# Index of radial bins that fitting starts at (if kep_fit == "TRUE")
fit_stop = 15		# As above, but final index of fitting
#
RL_fitcheck = "FALSE"	# Choose whether ("TRUE") or not ("FALSE") to compare masses computed to specfifc quadrants
#
raw_fit = "FALSE"	# Choose whether ("TRUE") or not ("FALSE") to plot raw points to which Keplerian fit runs 
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
import scipy
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def pv(dat_dir, plotdir, ea_run, snaparr, snapcore, v_K, inclin, r, vkep, \
   pmass, EA_lenref, EA_timeref, pcAU, AUm, G, Msol_kg):
#
	print "Position-Velocity diagram now being plotted"
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Set filenames to be read in for PV plotting/analysis, dependent on EA run being analysed
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# For EA [0, 1] runs
#
	if (snaparr.ndim == 1):
		file_n = len(snaparr)
		pv_file = [""]*file_n
		Lfit_file = [""]*file_n ; Rfit_file = [""]*file_n
#
		snaparr_tmp = [0]*file_n ; snapcore_tmp = [0]*len(snapcore)
		pmass_tmp = [0.]*len(snapcore)
		
		plot_title = [""]*len(snapcore) ; plot_filename = [""]*len(snapcore)
		fcount = 0
		for i in range(0, len(snaparr)):
#
			if (snaparr[i] < (1000-70)):
				pv_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.hist_pv.1'
				Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.fit_Lpv.1'
				Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.fit_Rpv.1'
			elif (snaparr[i] > (1000-70)):
				pv_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.hist_pv.1'
				Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.fit_Lpv.1'
				Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.fit_Rpv.1'
			for j in range(0, len(snapcore)):
				if (snaparr[i] == snapcore[j]):
					snapcore_tmp[j] = fcount
					plot_title[fcount] = 'P-V Diagram: '+str(ea_run)+'   ('+str(snaparr[i])+')'
					plot_filename[fcount] = plotdir+'PV_diag_'+str(snaparr[i])+'.pdf'
#
					for a in range(1, len(pmass)):
						pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i]
					fcount = fcount + 1
#
	# For EA [2, 3, 4, 5, 6] runs
#
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
		pv_file = [""]*file_n
		Lfit_file = [""]*file_n
		Rfit_file = [""]*file_n
#
		snaparr_tmp = [0]*file_n ; pmass_tmp = [0.]*file_n
		plot_title = [""]*file_n ; plot_filename = [""]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snaparr[0])):
#
				if (snaparr[i][j] < (1000-70)):
					pv_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.hist_pv.1'
					Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
					Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
				elif (snaparr[i][j] > (1000-70)):
					pv_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.hist_pv.1'
					Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
					Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
#
				plot_title[fcount] = 'P-V Diagram: ('+str(EA_lenref[i])+', '+ \
				   str(EA_timeref[j])+')'
				plot_filename[fcount] = plotdir+'PV_diag_'+str(EA_lenref[i])+'_' \
				   +str(EA_timeref[j])+'.pdf'
#
				snaparr_tmp[fcount] = snaparr[i][j]
#
#				for a in range(1, len(pmass)-1):
#					pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i][j]
#
				fcount = fcount + 1
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Loop over list of snapshots and read in data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
	# First, define arrays where PV diagram mass and associated error are stored
#
	pv_mass = [] ; pv_mass_err = []
#
	for i in range(0, fcount):
#
		r_bin = [] ; vz_bin = [] ; count_bin = []
#
	# Read data files, for PV diagram data, and data from which Keplerian fits are made
#	
		f = open(pv_file[i], 'r')
		r_range = float(f.readline()) ; v_range = float(f.readline())
		delt_r = float(f.readline()) ; delt_v = float(f.readline())
		y_min = float(f.readline()) ; y_max = float(f.readline())
		part_mass = float(f.readline())
		for line in f:
			line = line.strip() ; columns = line.split()
			r_bin.append(float(columns[0]))
			vz_bin.append(float(columns[1]))
			count_bin.append(float(columns[2]))
		f.close()
#
		r_bin = np.array(r_bin) ; vz_bin = np.array(vz_bin)
		count_bin = np.array(count_bin)
#
	# Select whether or not Keplerian fits are made to PV data
#
		if (kep_fit == "TRUE"):
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Calculate noise level based off in-situ PV diagram envelope emission
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			noise_count = 0
			count_sum = 0
#
			for j in range(0, len(count_bin)):
				if ((count_bin[j] != 0) and (r_bin[j] < -50.) \
				   and (vz_bin[j] > 0.5)):
					noise_count = noise_count + 1
					count_sum = count_sum + count_bin[j]
				if ((count_bin[j] != 0) and (r_bin[j] > 50.) \
				   and (vz_bin[j] < -0.5)):
					noise_count = noise_count + 1
					count_sum = count_sum + count_bin[j]
#
			sig_noise = sig_mult * (count_sum / noise_count)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# LHS PV diagram analysis
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			Lfit_arr = [[],[],[]]
#
			f = open(Lfit_file[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				Lfit_arr[0].append(float(columns[0]))
				Lfit_arr[1].append(float(columns[1]))
				Lfit_arr[2].append(float(columns[2]))
			f.close()
#	
			r_iter = int(math.ceil(r_range / delt_r))
			v_iter = int(math.ceil(v_range / delt_v))
#	
	# Define radial points of analysis, to begin evaluating y-positions of fits
#
			r_arr = [0.]*r_iter	
			rcount = 0
			for a in range(0, r_iter):
				rcount = rcount + 1
#
	# Determine the velocity value of Keplerian edge for fitting procedure
	# Ensure radial bins not at R = 0, as fitting routine crashes
#
				r_arr[a] = ((rcount-1)*delt_r) + 0.5*delt_r - r_range
				Lfit_y = []
				Lfit_y0 = []
				Lfit_r = []
#
				for a in range(0, r_iter):
					Lfit_ydum = []
					Lfit_ndum = []
					Lfit_ystart = []
					Lfit_sig = []
					for b in range(0, len(Lfit_arr[0])):
						if (Lfit_arr[0][b] == r_arr[a]):
							Lfit_ydum.append(Lfit_arr[1][b])
							Lfit_ndum.append(Lfit_arr[2][b])
					for b in range(0, len(Lfit_ydum)):
						if (Lfit_ndum[b] > sig_noise ):
							Lfit_y.append(Lfit_ydum[b])
							Lfit_r.append(r_arr[a])
							break
					for b in range(0, len(Lfit_ydum)):
						if (Lfit_ndum[b] > 0):
							Lfit_y0.append(Lfit_ydum[b])
							break
#
	# Compute sigma value for errors on curve fitting procedure, by defining absolute edge of
	# PV diagram and calculating rough sigma from difference between sigma_noise bin
#	
			Lfit_sig = []
			for a in range(0, len(Lfit_y)):
				if (Lfit_y[a] == Lfit_y0[a]):
					Lfit_sig.append(0.5*delt_v)
				else:
					Lfit_sig.append( 0.5*( abs(Lfit_y[a] - Lfit_y0[a] ) ) )	
#
	# Change quadrant of LHS fit arrays for correct power index fitting to curve
#
			Lfit_r = [-k for k in Lfit_r] ; Lfit_y = [-k for k in Lfit_y]
			Lfit_r = Lfit_r[::-1] ; Lfit_y = Lfit_y[::-1]
#
	# Define function of Keplerian fit
#	
			def func(x, a):
			  return a* np.power(x, -0.5)
#
	# Fitting routine to PV-diagram, and determination of system mass using Keplerian edge fit
#
			coeffs1, fitcov1 = scipy.optimize.curve_fit(func, Lfit_r[fit_start:fit_stop], \
			   Lfit_y[fit_start:fit_stop], [1.], \
			   sigma=Lfit_sig[fit_start:fit_stop])
#
			Lfit_r1 = np.array(Lfit_r)*AUm ; Lfit_y1 = np.array(Lfit_y)*1000.
			Lfit_sig1 = np.array(Lfit_sig)*1000.
			coeffs11, fitcov11 = scipy.optimize.curve_fit(func, Lfit_r1[fit_start:fit_stop], \
			   Lfit_y1[fit_start:fit_stop], [1e3], \
			   sigma=Lfit_sig1[fit_start:fit_stop])
			y_fitted1 = func(Lfit_r[fit_start:fit_stop], coeffs1[0])
#
			perr1 = np.sqrt(np.diag(fitcov11))
#
			M_sys1 = (coeffs11[0]**2.0)/(G*Msol_kg)
			M_err1 = (perr1[0]**2.0)/(G*Msol_kg) 
#
	# Re-convert arrays to original quadrant, and change
#
			Lfit_r = [-k for k in Lfit_r] ; Lfit_y = [-k for k in Lfit_y]
			Lfit_r = Lfit_r[::-1] ; Lfit_y = Lfit_y[::-1]
#
			y_fitted1 = [-k for k in y_fitted1]
			y_fitted1 = y_fitted1[::-1]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# RHS PV diagram analysis
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			Rfit_arr = [[],[],[]]
#
			f = open(Rfit_file[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				Rfit_arr[0].append(float(columns[0]))
				Rfit_arr[1].append(float(columns[1]))
				Rfit_arr[2].append(float(columns[2]))
			f.close()
#
			r_iter = int(math.ceil(r_range / delt_r))
			v_iter = int(math.ceil(v_range / delt_v))
#
			r_arr = [0.]*r_iter	
			rcount = 0
			for a in range(0, r_iter):
				rcount = rcount + 1
#
				r_arr[a] = ((rcount-1)*delt_r) + 0.5*delt_r
				Rfit_y = []
				Rfit_y0 = []
				Rfit_r = []
#
				for a in range(0, r_iter):
					Rfit_ydum = []
					Rfit_ndum = []
					for b in range(0, len(Rfit_arr[0])):
						if (Rfit_arr[0][b] == r_arr[a]):
							Rfit_ydum.append(Rfit_arr[1][b])
							Rfit_ndum.append(Rfit_arr[2][b])
					ind = len(Rfit_ndum)-1
					for b in range(0, len(Rfit_ydum)):
						if (Rfit_ndum[ind] > sig_noise):
							Rfit_y.append(Rfit_ydum[ind])
							Rfit_r.append(r_arr[a])
							break
						ind = ind - 1
					ind = len(Rfit_ndum)-1
					for b in range(0, len(Rfit_ydum)):
						if (Rfit_ndum[ind] > 0):
							Rfit_y0.append(Rfit_ydum[ind])
							break
						ind = ind - 1
#
	# Compute sigma value for errors on curve fitting procedure, by defining absolute edge of
	# PV diagram and calculating rough sigma from difference between sigma_noise bin
#	
			Rfit_sig = []
			for a in range(0, len(Rfit_y)):
				if (Rfit_y[a] == Rfit_y0[a]):
					Rfit_sig.append(0.5*delt_v)
				else:
					Rfit_sig.append( 0.5*( abs(Rfit_y[a] - Rfit_y0[a] ) ) )
#
	# Fitting routine to quadrant data
#
			coeffs2, fitcov2 = scipy.optimize.curve_fit(func, Rfit_r[fit_start:fit_stop], \
			   Rfit_y[fit_start:fit_stop], [1.], \
			   sigma=Rfit_sig[fit_start:fit_stop])
			y_fitted2 = func(Rfit_r[fit_start:fit_stop], coeffs2[0])
#
			Rfit_r2 = np.array(Rfit_r)*AUm ; Rfit_y2 = np.array(Rfit_y)*1000.
			Rfit_sig2 = np.array(Rfit_sig)*1000.
			coeffs22, fitcov22 = scipy.optimize.curve_fit(func, Rfit_r2[fit_start:fit_stop], \
			   Rfit_y2[fit_start:fit_stop], [1e3], \
			   sigma=Rfit_sig2[fit_start:fit_stop])
#
			perr2 = np.sqrt(np.diag(fitcov22))
#
			M_sys2 = (coeffs22[0]**2.0)/(G*Msol_kg)
			M_err2 = (perr2[0]**2.0)/(G*Msol_kg)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Outputs of PV diagram mass analysis
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Mass average
#
			pv_mass.append( (M_sys1 + M_sys2) / 2. )
			pv_mass_err.append( (M_err1 + M_err2) / 2. )
#
	# Quadrant mass/error checks
#
			if (RL_fitcheck == "TRUE"):
				print "Mass of the system determined from LHS PV fit is: ", \
				   round(M_sys1,4), " Solar masses"
				print "Mass of the system determined from RHS PV fit is: ", \
				   round(M_sys2,4), " Solar masses"
				print "Averaged mass from PV diagram fits is: ", \
				   round(pv_mass[i],4), "Solar masses"
				print "Error on mass determined from PV diagram is: ", \
				   round(pv_mass_err[i],5), "Solar masses"
#
	# Estimation of companion mass vs. sink mass	
#
			if (pmass_tmp[i] > 0):
				companion_mass = abs(M_sys1 - M_sys2)
				discrep = ( abs(companion_mass - pmass_tmp[i]) / pmass_tmp[i] ) * 100.			
				print "\n Estimate of companion mass from kinematic analysis is: ", \
				   round(companion_mass,4), " Solar masses"
				print "Actual sink particle mass is: ", round(pmass_tmp[i],4), "\n"
				print "% discrepancy between kinematic/sink companion mass: ", \
				   round(discrep,4), " Solar masses"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of PV diagram - EA RUNS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Set up grid of interpolation points, interpolate on arbitrary grid, and plot
#
		if (snaparr.ndim == 2):
#
			xi, yi = np.linspace(r_bin.min(), r_bin.max(), 500), \
			   np.linspace(vz_bin.min(), vz_bin.max(), 500)
			xi, yi = np.meshgrid(xi, yi)
#
	# Before defining the2D histogram density, convert number density to (logarithmic) 
	# surface density by computing cgs conversions for area and mass, then convert number 
	# histogram to surface density in PV diagram radius-velocity bins
#
			sd_norm = []
			bin_area = (delt_r*(AUm*100.)) * ( (y_max-y_min) * (AUm*100.))	# area in cm^2
			part_mass = part_mass * (Msol_kg * 1000.)			# mass in g
			for j in range(0, len(count_bin)):
				if (count_bin[j] != 0):
					sd_norm.append( (count_bin[j]*part_mass) / bin_area)
				if (count_bin[j] == 0):
					sd_norm.append(1e-10)  
			sd_norm = np.log10(np.array(sd_norm))
#
			zi = scipy.interpolate.griddata((r_bin, vz_bin), sd_norm, (xi, yi), \
			   method='nearest')
#
			plt.figure(1)
			ax1 = plt.subplot(111)
			plt.imshow(zi, vmin=-1., vmax=2.5, origin='lower', \
			   extent=[r_bin.min(), r_bin.max(), vz_bin.min(), vz_bin.max()], aspect = 'auto')
#
	# Overplot Keplerian fits to data, plotting points fitted to if chosen
#
#
			if (kep_fit == "TRUE"):
#
	# PV Keplerian distributions
#
				if (raw_fit == "TRUE"):
					Lraw = plt.plot(Lfit_r, Lfit_y, linewidth = 4, \
					   linestyle = 'solid', color = 'g')		
					Rraw = plt.plot(Rfit_r, Rfit_y, linewidth = 4, \
					   linestyle = 'solid', color = 'g')		
#
				if (raw_fit == "FALSE"):
					Lfit_start = len(Lfit_r)-fit_stop
					Lfit_stop = len(Lfit_r)-fit_start
					Lfit_start = len(Lfit_r)-fit_stop + \
					   (len(Lfit_r[Lfit_start:Lfit_stop]) - len(y_fitted1))
#
					Lfit = plt.plot(Lfit_r[Lfit_start:Lfit_stop], y_fitted1, \
					   linewidth = 4, linestyle = 'solid', color = 'r')
					Rfit = plt.plot(Rfit_r[fit_start:fit_stop], y_fitted2, \
					   linewidth = 4, linestyle = 'solid', color = 'r')
#
	# Simulation Keplerian distributions
#
					plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
					   linestyle = 'dashed', color = 'k')
					r[i] = [-k for k in r[i]] ; vkep[i] = [-k for k in vkep[i]]
					plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
					   linestyle = 'dashed', color = 'k')
#
			cbar = plt.colorbar() ; cbar.ax.set_ylabel('Log. Surface density,'+' (g '+(r'cm$^{-2}$')+')')
			plt.xlabel('Radius (AU)') ; plt.ylabel('Line of sight velocity'+' (km'+(r's$^{-1}$')+')')
			plt.xlim(-100,100) ; plt.ylim(-10, 10)
#
			plt.savefig(plot_filename[i]) ; plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of PV diagram - NON EA RUNS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Set up grid of interpolation points, interpolate on arbitrary grid, and plot
#
		if (snaparr.ndim == 1):
#
				for jdash in range(0, len(snapcore_tmp)):
					if (snaparr_tmp[i] == snapcore_tmp[jdash]):
#
						xi, yi = np.linspace(r_bin.min(), r_bin.max(), 500), \
						   np.linspace(vz_bin.min(), vz_bin.max(), 500)
						xi, yi = np.meshgrid(xi, yi)
#
						sd_norm = []
						bin_area = (delt_r*(AUm*100.)) * ( (y_max-y_min) * (AUm*100.))
						part_mass = part_mass * (Msol_kg * 1000.)
						for j in range(0, len(count_bin)):
							if (count_bin[j] != 0):
								sd_norm.append( (count_bin[j]*part_mass) / bin_area)
							if (count_bin[j] == 0):
								sd_norm.append(1e-10)  
						sd_norm = np.log10(np.array(sd_norm))
#
						zi = scipy.interpolate.griddata((r_bin, vz_bin), sd_norm, (xi, yi), \
						   method='nearest')
#
						plt.figure(1)
						ax1 = plt.subplot(111)
						plt.imshow(zi, vmin=-1., vmax=2.5, origin='lower', \
						   extent=[r_bin.min(), r_bin.max(), vz_bin.min(), vz_bin.max()], aspect = 'auto')
#
						if (kep_fit == "TRUE"):
#
							if (raw_fit == "TRUE"):
								Lraw = plt.plot(Lfit_r, Lfit_y, linewidth = 4, \
								   linestyle = 'solid', color = 'g')		
								Rraw = plt.plot(Rfit_r, Rfit_y, linewidth = 4, \
								   linestyle = 'solid', color = 'g')		
#
							if (raw_fit == "FALSE"):
								Lfit_start = len(Lfit_r)-fit_stop
								Lfit_stop = len(Lfit_r)-fit_start
								Lfit_start = len(Lfit_r)-fit_stop + \
								   (len(Lfit_r[Lfit_start:Lfit_stop]) - len(y_fitted1))
#
								Lfit = plt.plot(Lfit_r[Lfit_start:Lfit_stop], y_fitted1, \
								   linewidth = 4, linestyle = 'solid', color = 'r')
								Rfit = plt.plot(Rfit_r[fit_start:fit_stop], y_fitted2, \
								   linewidth = 4, linestyle = 'solid', color = 'r')		
#
								plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
								   linestyle = 'dashed', color = 'k')
								r[i] = [-k for k in r[i]] ; vkep[i] = [-k for k in vkep[i]]
								plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
								   linestyle = 'dashed', color = 'k')
#
						cbar = plt.colorbar() ; cbar.ax.set_ylabel('Log. Surface density,'+' (g '+(r'cm$^{-2}$')+')')
						plt.xlabel('Radius (AU)') ; plt.ylabel('Line of sight velocity'+' (km'+(r's$^{-1}$')+')')
						plt.xlim(-100,100) ; plt.ylim(-10, 10)
#
						plt.savefig(plot_filename[jdash]) ; plt.clf()
#	
#
	return pv_mass, kep_fit, raw_fit

#
# pv_diag.py
#
# Programme designed to read in simulation data and generate synthetic PV diagram
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
sig_mult = 0        # Threshold sigma value over which edge of Keplerian profile is fitted to
#
kep_fit = True    # Choose whether to ***fit*** Keplerian profiles to PV data, to compute kinematic mass
kep_plot = False     # Choose whether to ***plot*** Keplerian profiles (raw or fitted) onto PV diagram
ana_kep = False   # Choose whether to over plot analytic Keplerian profile on PV plot
#
fit_start = 0        # Index of radial bins that fitting starts at (if kep_fit == "TRUE")
fit_stop = [[10,16,18],[6,20,20],\
   [20,20,20,20,12,16,12,16,16,16,20,16],[20,20,20,20,17,15,20,20],[20,20,20,20,20,13,17,11]]
#
RL_fitcheck = "FALSE"    # Choose whether ("TRUE") or not ("FALSE") to compare masses computed to specfifc quadrants
#
raw_fit = "FALSE"    # Choose whether ("TRUE") or not ("FALSE") to plot raw points to which Keplerian fit runs
#
arrows = True        # Decide whether keplerian fits are plotted, or arrows in ea5 (for publication Fig. 11)
#
intensity = True    # Decide whether or not flux scale is used
isotherm = False      # If intensity is True, choose whether to use actual, or isothermal assumption in flux calculation
#
Jynorm = True          # If called for, show plot in terms of Jy/beam
Jy = 1.e26           # milli-Jansky value
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
import scipy
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def pv(dat_dir, plotdir, ea_run, snaparr, v_K, inclin, r, vkep, pmass, EA_lenref, \
   EA_timeref, pcAU, AUm, G, Msol_kg, pv_MASK60, pv_MASK90, pcm, mH, plot_form, \
   kb, c, h, cgsSIopac, mu, Jerg, cm2m2, part_mass):
#
    print "Position-Velocity diagram now being plotted"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Set filenames to be read in for PV plotting/analysis, dependent on EA run being analysed
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    if (ea_run <= 1):
        ea_run_IND = ea_run
    elif (ea_run > 1):
        ea_run_IND = ea_run -1
#
    # For EA [2, 3, 4, 5, 6] runs
#
    if (ea_run <= 1):
        file_n = len(snaparr)
    else:
        file_n = len(snaparr)*len(snaparr[0])
#
    pv_file = [""]*file_n
    Lfit_file = [""]*file_n
    Rfit_file = [""]*file_n
#
    snaparr_tmp = [0]*file_n ; pmass_tmp = [0.]*file_n
    plot_title = [""]*file_n ; plot_filename = [""]*file_n
    fcount = 0

    if (ea_run <= 1):
        for i in range(0, len(snaparr)):
            if (snaparr[i] < (1000-70)):
                pv_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i]+69)+'.du.hist_pv.1'
                Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i]+69)+'.du.fit_Lpv.1'
                Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i]+69)+'.du.fit_Rpv.1'
            elif (snaparr[i] > (1000-70)):
                pv_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i]+69)+'.du.hist_pv.1'
                Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i]+69)+'.du.fit_Lpv.1'
                Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i]+69)+'.du.fit_Rpv.1'
#
            plot_title[fcount] = 'P-V Diagram: Snapshot '+str(snaparr[i])
            plot_filename[fcount] = plotdir+'PV_diag_'+str(snaparr[i])+'.'+str(plot_form)
            snaparr_tmp[fcount] = snaparr[i]
            fcount = fcount + 1

    else:
        for i in range(0, len(snaparr)):
            for j in range(0, len(snaparr[0])):
#
                if (snaparr[i][j] < (1000-70)):
                    pv_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i][j]+69)+'.du.hist_pv.1'
                    Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
                    Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.00'+str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
                elif (snaparr[i][j] > (1000-70)):
                    pv_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i][j]+69)+'.du.hist_pv.1'
                    Lfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
                    Rfit_file[fcount] = dat_dir+'pv_diag/DE05.du.0'+str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
#
                plot_title[fcount] = 'P-V Diagram: ('+str(EA_lenref[i])+', '+str(EA_timeref[j])+')'
                plot_filename[fcount] = plotdir+'PV_diag_'+str(EA_lenref[i])+'_'+str(EA_timeref[j])+'.'+str(plot_form)
#
                snaparr_tmp[fcount] = snaparr[i][j]
                fcount = fcount + 1
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Loop over list of snapshots and read in data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    for i in range(0, fcount):
#
        r_bin = [] ; vz_bin = [] ; temp_bin = [] ; count_bin = []
#
    # Read data files, for PV diagram data, and data from which Keplerian fits are made
#
        f = open(pv_file[i], 'r')
        r_range = float(f.readline()) ; v_range = float(f.readline())
        delt_r = float(f.readline()) ; delt_v = float(f.readline())
        y_min = float(f.readline()) ; y_max = float(f.readline())
        trash = float(f.readline())
        for line in f:
            line = line.strip() ; columns = line.split()
            r_bin.append(float(columns[0]))
            vz_bin.append(float(columns[1]))
            if (float(columns[2]) == 0.):
                temp_bin.append(10.)
            else:
                temp_bin.append(float(columns[2]))
            count_bin.append(float(columns[3]))
        f.close()
#
        r_bin = np.array(r_bin) ; vz_bin = np.array(vz_bin)
        temp_bin = np.array(temp_bin) ; count_bin = np.array(count_bin)
#
    # Select whether or not Keplerian fits are made to PV data
#
        if kep_fit is True:
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
            coeffs1, fitcov1 = scipy.optimize.curve_fit(func, Lfit_r[fit_start:fit_stop[ea_run_IND][i]], \
               Lfit_y[fit_start:fit_stop[ea_run_IND][i]], [1.], \
               sigma=Lfit_sig[fit_start:fit_stop[ea_run_IND][i]])
#
            Lfit_r1 = np.array(Lfit_r)*AUm ; Lfit_y1 = np.array(Lfit_y)*1000.
            Lfit_sig1 = np.array(Lfit_sig)*1000.
            coeffs11, fitcov11 = scipy.optimize.curve_fit(func, Lfit_r1[fit_start:fit_stop[ea_run_IND][i]], \
               Lfit_y1[fit_start:fit_stop[ea_run_IND][i]], [1e3], \
               sigma=Lfit_sig1[fit_start:fit_stop[ea_run_IND][i]])
            y_fitted1 = func(Lfit_r[fit_start:fit_stop[ea_run_IND][i]], coeffs1[0])
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
            coeffs2, fitcov2 = scipy.optimize.curve_fit(func, Rfit_r[fit_start:fit_stop[ea_run_IND][i]], \
               Rfit_y[fit_start:fit_stop[ea_run_IND][i]], [1.], \
               sigma=Rfit_sig[fit_start:fit_stop[ea_run_IND][i]])
            y_fitted2 = func(Rfit_r[fit_start:fit_stop[ea_run_IND][i]], coeffs2[0])
#
            Rfit_r2 = np.array(Rfit_r)*AUm ; Rfit_y2 = np.array(Rfit_y)*1000.
            Rfit_sig2 = np.array(Rfit_sig)*1000.
            coeffs22, fitcov22 = scipy.optimize.curve_fit(func, Rfit_r2[fit_start:fit_stop[ea_run_IND][i]], \
               Rfit_y2[fit_start:fit_stop[ea_run_IND][i]], [1e3], \
               sigma=Rfit_sig2[fit_start:fit_stop[ea_run_IND][i]])
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
            pv_mass_err = [] ; pv_mass60 = pv_MASK60 ; pv_mass90 = pv_MASK90
            if (inclin == "60"):
                pv_mass60.append( (M_sys1 + M_sys2) / 2. )
                pv_mass90 = []
            if (inclin == "90"):
                pv_mass90.append( (M_sys1 + M_sys2) / 2. )
            pv_mass_err.append( (M_err1 + M_err2) / 2. )
#
    # Quadrant mass/error checks
#
            if RL_fitcheck is True:
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
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Generation of z-axis values
    #
    # Note that due to difference in SI and cgs treatments, part_mass and bin_area
    # computation are slightly different in flux and surface density based PV diagrams
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    # Common parameters, either for flux or surface density PV diagrams
#
        dist = 140. * pcm        # Distance to source (in pc, converted to m)
        d_to_g = 0.01           # Assumed dust-to-gas ratio
        lam = 1.3e-3            # Assumed lambda value, equivalent to ALMA band 6
        nu = 230.e9             # As above, for frequency
        kappa_lam = 2.3 * cgsSIopac      # Assumed opacity for disc, in cgs, converted to SI
        bin_area = (delt_r*(AUm)) * ( (y_max-y_min) * (AUm))    # area in m^2
        delt_omega = bin_area / dist**2.     # Solid angle for resolved area (SI)
        gas_pmass = part_mass * Msol_kg            # gas particle mass in kg
        B_nu = []              # Blank array for spectral irradiance values to be filled in flux PV diagrams
#
    # For isothermal assummed observation, computing flux
#
#
        if (intensity is True) and (isotherm is True):
#
            if (ea_run == 0):               # Dependent on either RF and outburst state
                t_iso = 30.                 # for ERF models, set assumed temperature value
            elif (ea_run == 1):             # - motivated by average azimuthal values in-disc
                t_iso = 100.
            else:
                for ii in range(0, len(snaparr)):
                    for jj in range(0, len(snaparr[0])):
                        if (jj == 0) or (jj == 3):
                            t_iso = 30.
                        else:
                            t_iso = 100.
#
            for ii in range(0, len(count_bin)):
                B_nu.append( ( (2. * h * nu**3. ) / c**2. ) * \
                   ( 1. / (math.exp( (h*nu) / (kb*t_iso) ) - 1.) ) )  # Planck function for t_iso
#
    # For non-isothermal case (T computed as average in projected v_los and r bin)
#
        if (intensity is True) and (isotherm is False):
#
            print np.mean(temp_bin), min(temp_bin), max(temp_bin)
#
            for ii in range(0, len(count_bin)):
                B_nu.append(( (2. * h * nu**3. ) / c**2. ) * \
                   ( 1. / (math.exp( (h*nu) / (kb * temp_bin[ii]) ) - 1.) ) )
#
    # For surface density case
#
        elif intensity is False:
#
            gas_pmass = part_mass * (Msol_kg * 1000.)            # mass in g
            bin_area = (delt_r*(AUm*100.)) * ( (y_max-y_min) * (AUm*100.))    # area in cm^2
#
#
    # Loop over files, computing the z-axis value based on whether set by Flux, or surface density
#
#
        zaxis = []
#
        for j in range(0, len(count_bin)):
#
            if intensity is True:
#
                if (count_bin[j] != 0):
#
                    sd_dust = ( (count_bin[j]*gas_pmass) / bin_area ) * d_to_g
#
                    if Jynorm is True:
                        zaxis.append(  B_nu[j] * kappa_lam * sd_dust * Jy * delt_omega )
                    else:
                        zaxis.append( B_nu[j] * Jerg * cm2m2 * kappa_lam * sd_dust * delt_omega)
#
                if (count_bin[j] == 0):
                    zaxis.append(1.e-99)

#
            if intensity is False:
                if (count_bin[j] != 0):
                    zaxis.append( (count_bin[j]*gas_pmass) / bin_area)
                if (count_bin[j] == 0):
                    zaxis.append(1e-99)
#
        zaxis = np.log10(np.array(zaxis))
#        zaxis = np.array(zaxis)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Plotting of PV diagram
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    # Set up grid of interpolation points, interpolate on arbitrary grid, and plot
#
        xi, yi = np.linspace(r_bin.min(), r_bin.max(), 500), \
           np.linspace(vz_bin.min(), vz_bin.max(), 500)
        xi, yi = np.meshgrid(xi, yi)
#
        zi = scipy.interpolate.griddata((r_bin, vz_bin), zaxis, (xi, yi), \
           method='nearest')
#
        zmax = zaxis.max() ; zmin = zaxis.max() - 3.
#
        plt.figure(1)
        ax1 = plt.subplot(111)
        plt.imshow(zi, vmin=zmin, vmax=zmax, origin='lower', \
           extent=[r_bin.min(), r_bin.max(), vz_bin.min(), vz_bin.max()], \
           aspect = 'auto')
#
    # Overplot Keplerian fits to data, plotting points fitted to if chosen
#
        if kep_plot is True:
#
    # PV Keplerian distributions
#
            if raw_fit is True:
                Lraw = plt.plot(Lfit_r, Lfit_y, linewidth = 4, \
                   linestyle = 'solid', color = 'g')
                Rraw = plt.plot(Rfit_r, Rfit_y, linewidth = 4, \
                   linestyle = 'solid', color = 'g')
#
            else:
                Lfit_start = len(Lfit_r)-fit_stop[ea_run_IND][i]
                Lfit_stop = len(Lfit_r)-fit_start
                Lfit_start = len(Lfit_r)-fit_stop[ea_run_IND][i] + \
                   (len(Lfit_r[Lfit_start:Lfit_stop]) - len(y_fitted1))
#
                Lfit = plt.plot(Lfit_r[Lfit_start:Lfit_stop], y_fitted1, \
                   linewidth = 4, linestyle = 'solid', color = 'r')
                Rfit = plt.plot(Rfit_r[fit_start:fit_stop[ea_run_IND][i]], y_fitted2, \
                   linewidth = 4, linestyle = 'solid', color = 'r')
#
    # Simulation Keplerian distributions
#
                if ana_kep is True:
                    plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
                    linestyle = 'dashed', color = 'k')
                    r[i] = [-k for k in r[i]] ; vkep[i] = [-k for k in vkep[i]]
                    plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 4, \
                    linestyle = 'dashed', color = 'k')
#
        if arrows is True:
            if (ea_run == 5) and (i == 4) and isothermal is True:
                ax1.arrow(-50., -7.5, 0., 5., head_width=5., head_length=1., fc='g', ec='g', width=2.)
                ax1.arrow(40., 7.5, 0., -4.5, head_width=5., head_length=1., fc='g', ec='g', width=2.)
            if (ea_run == 5) and (i == 5):
                ax1.arrow(-50., -7.5, 0., 4.5, head_width=5., head_length=1., fc='g', ec='g', width=2.)
                ax1.arrow(-60., -7.5, 0., 3., head_width=5., head_length=1., fc='b', ec='b', width=2.)
                ax1.arrow(40., 7.5, 0., -4.5, head_width=5., head_length=1., fc='g', ec='g', width=2.)
            if (ea_run == 5) and (i == 6):
                ax1.arrow(-50., -7.5, 0., 5., head_width=5., head_length=1., fc='g', ec='g', width=2.)
                ax1.arrow(-80., -7.5, 0., 3., head_width=5., head_length=1., fc='b', ec='b', width=2.)
                ax1.arrow(50., 7.5, 0., -5., head_width=5., head_length=1., fc='g', ec='g', width=2.)
                ax1.arrow(60., 7.5, 0., -3., head_width=5., head_length=1., fc='b', ec='b', width=2.)
#
        cbar = plt.colorbar()
        if Jynorm is True:
            cbar.ax.set_ylabel('Log. F'+(r'$_\nu$')+' (Jy)')
        else:
            if intensity is True:
                cbar.ax.set_ylabel('Log. Flux,'+' (erg '+(r'cm$^{-2}$')+')')
            if intensity is False:
                cbar.ax.set_ylabel('Log. Surface density,'+' (g '+(r'cm$^{-2}$')+')')
        plt.xlabel('Radius (AU)') ; plt.ylabel('Line of sight velocity'+' (km'+(r's$^{-1}$')+')')
        plt.xlim(-100,100) ; plt.ylim(-10, 10)
#
        plt.savefig(plot_filename[i], format=str(plot_form), dpi=150) ; plt.clf()
#
#
    return pv_mass60, pv_mass90, kep_fit, raw_fit

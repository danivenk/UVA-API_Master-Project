import numpy as np
import scipy
import scipy.stats as sps
import scipy.fftpack as spfft
import scipy.signal as spsig
import scipy.interpolate as spinterp
import math
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

def find_nearest(array,value):
    """ Returns nearest value in an array and its index. """
    # Courtesy of Stack Overflow: 
    # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    idx = np.abs(array-value).argmin()
    return array[idx], idx

def calc_dispfrac(rad,rad_area,rin,rcor,seedff_norm,seedff_ind,heatff_norm,heatff_ind):
    """ Calculates viscous dissipation for each disk annulus normalised by the total summed over all radii. 
        Also calculates the fraction of dissipated power from radii inside the corona which go into heating 
        and seed photons. The sum of heating and seed luminosity fraction is not automatically kept equal
        to 1, to allow for unconstrained heating/cooling, e.g. from energy transport to the flow from outside 
        the coronal radius, but unity total fraction can be forced by specifying seedff_norm < 0, although
        care must be taken to ensure that heat_frac_flow remains < 1 at all radii.
        For seedff_norm = 0 no seed luminosity will be produced inside the corona."""
        
    dispfrac = rad_area*(rad**(-3))*(1.-np.sqrt(rin/rad))   # Assumes zero-torque inner boundary
    dispfrac = dispfrac/np.sum(dispfrac)
    seed_frac_flow = np.zeros(len(rad))
    heat_frac_flow = np.zeros(len(rad))
    
    heat_frac_flow[rad <= rcor] = heatff_norm*(rad[rad <= rcor]/rin)**heatff_ind
    
    if (seedff_norm >= 0):
        seed_frac_flow[rad <= rcor] = seedff_norm*(rad[rad <= rcor]/rin)**heatff_ind
    else:
        seed_frac_flow[rad <= rcor] = 1.0-heat_frac_flow[rad <= rcor]
        
    return dispfrac, seed_frac_flow, heat_frac_flow
    
    

def calc_illumination_fracs(rad,rad_area,geomod,params):
    """ Front-end to function which outputs the coronal solid angle and 'illumination fractions':
    fraction of disk emission intercepted by the corona and fraction of coronal emission
    intercepted by the disk, for an array of input radii and the areas of the associated annuli.
    The calculation depends on the input model name (geomod) and its parameters, rcor (coronal radius)
    and auxiliary parameters in params. """
    
    rcor = params[0]
    omega_cor = np.zeros(len(rad))
    frad_disktocor = np.zeros(len(rad))
    frad_cortodisk = np.zeros(len(rad))
    for i in range(0, len(rad)):
        if (rad[i] > rcor):
            omega_cor[i], frad_disktocor[i], frad_cortodisk[i] = geomod(rad[i],params)
    frad_cortodisk = np.multiply(frad_cortodisk,rad_area)
    if (geomod == bknpow_emiss):
        total_frac = np.sum(frad_cortodisk)
        frad_cortodisk = frad_cortodisk*params[6]/total_frac
        frad_disktocor = frad_disktocor*params[6]/total_frac
    if (geomod == piecewise_emiss):
        total_frac = np.sum(frad_cortodisk)
        frad_cortodisk = frad_cortodisk*params[2]/total_frac
        
    return omega_cor, frad_disktocor, frad_cortodisk


def bknpow_emiss(r,params):
    """ Broken power-law empirical function for illumination fractions.
    Solid angle set to zero (not used anyway). """
    [r_cor,r_bk1,r_bk2,em_ind1,em_ind2,em_ind3,total_cortodisk] = params
    if (r < r_bk1):
        f_cortodisk = (r/r_bk1)**em_ind1
    elif (r >= r_bk1 and r < r_bk2):
        f_cortodisk = (r/r_bk1)**em_ind2
    elif (r >= r_bk2):
        f_cortodisk = ((r_bk2/r_bk1)**em_ind2)*(r/r_bk2)**em_ind3
    f_disktocor = f_cortodisk*2.*np.pi*np.square(r_cor)
    omega = 0.0
    return omega, f_disktocor, f_cortodisk

def piecewise_emiss(r,params):
    """ Piece-wise power-law empirical function for illumination fractions.
    Solid angle set to zero (not used anyway). """
    r_cor = params[0]
    r_outer = params[1]
    i = 3
    while i+2 < len(params):
        r_inner = params[i]
        em_ind = params[i+1]
        norm = params[i+2]
        if ((r < r_outer) & (r >= r_inner)):
            f_disktocor = norm*(r/r_inner)**em_ind
        r_outer = r_inner
        i = i+3
    f_cortodisk = f_disktocor/(2.*np.pi*r_cor**2)
    omega = 0.0
    return omega, f_disktocor, f_cortodisk


def an_sphere(r,params):
    """ Quick analytical formulae for illumination fractions for spherical 'solid' corona of radius r_cor. """
    r_cor = params[0]
    f_disktocor = (np.arcsin(r_cor/r)-np.sqrt(np.square(r_cor/r)-np.power(r_cor/r,4.)))/np.pi # exact derivation
    f_cortodisk = f_disktocor/(2.*np.pi*np.square(r_cor))
    omega = 0.0 # Dummy output for solid angle, which we do not calculate here.
    return omega, f_disktocor, f_cortodisk

def sphere(r,params):
    """ Calculates (numerically) solid angle, illumination fractions for a spherical 'solid' corona with single
    parameter, the coronal radius r_cor. The model finds solid angle and illumination fractions by numerical
    integration using the approach given in Appendix in the UM20 paper. Number of cell spacings on latitudinal
    and azimuthal axes given by nphi and ntheta respectively (params[1] and params[2])."""
    r_cor = params[0]
    nphi = int(params[1])
    ntheta = int(params[2])
    dtheta = np.pi/(2.*ntheta)
    dphi = 2.*np.pi/nphi
    theta = np.linspace(dtheta/2,(np.pi-dtheta)/2.,ntheta)
    phi = np.linspace(dphi/2,2*np.pi-(dphi/2),nphi)
    phi_arr, theta_arr = np.meshgrid(phi, theta)
    x = np.multiply(np.multiply(np.sin(theta_arr), np.cos(phi_arr)),r_cor)
    y = np.multiply(np.multiply(np.sin(theta_arr), np.sin(phi_arr)),r_cor)
    z = np.multiply(np.cos(theta_arr), r_cor)
    da = np.multiply(np.square(r_cor)*dtheta*dphi, np.sin(theta_arr))
    xpos = r-x
    ypos = -1.*y
    zpos = -1.*z
    rpos = np.sqrt(np.square(xpos)+np.square(ypos)+np.square(zpos))
    xpos = xpos/rpos
    ypos = ypos/rpos
    zpos = zpos/rpos
    daomega = np.multiply(da,np.multiply(x,xpos)+np.multiply(y,ypos)+np.multiply(z,zpos))
    daomega = np.divide(daomega,(np.square(rpos)*r_cor))
    daomegawt = np.multiply(daomega,np.abs(zpos))
    omega = daomega[daomega > 0].sum()
    omegawt = daomegawt[daomegawt > 0].sum()
    f_disktocor = omegawt/np.pi  # Fraction of disk photons from this radius which hit the corona
    f_cortodisk = f_disktocor/(2.*np.pi*np.square(r_cor))  # Fraction of coronal photons
    # hitting a unit area of disk surface at this radius
    return omega, f_disktocor, f_cortodisk

def cylinder(r,params):
    """ Calculates (numerically) solid angle, illumination fractions for a cylindrical 'solid'
    corona with radius r_cor, height h_cor. The model finds solid angle and illumination fractions
    by numerical integration using the approach given in Appendix in the UM20 paper.  Number of
    cell spacings on vertical and azimuthal axes given by nz and nphi respectively
    (params[2] and params[3])."""
    r_cor = params[0]
    h_cor = params[1]
    nz = params[2]
    nphi = params[3]
    dz = h_cor/nz
    dphi = 2.*np.pi/nphi
    z = np.arange(0,h_cor,dz)
    phi = np.arange(0,2*np.pi,dphi)
    phi_arr, z_arr = np.meshgrid(phi, z)
    x = np.multiply(np.cos(phi_arr),r_cor)
    y = np.multiply(np.sin(phi_arr),r_cor)
    da = dz*dphi*r_cor
    xpos = r-x
    ypos = -1.*y
    zpos = -1.*z_arr
    rpos = np.sqrt(np.square(xpos)+np.square(ypos)+np.square(zpos))
    xpos = xpos/rpos
    ypos = ypos/rpos
    zpos = zpos/rpos
    daomega = np.multiply(da,np.multiply(x,xpos)+np.multiply(y,ypos))
    daomega = np.divide(daomega,(np.square(rpos)*r_cor))
    daomegawt = np.multiply(daomega,np.abs(zpos))
    omega = daomega[daomega > 0].sum()
    omegawt = daomegawt[daomegawt > 0].sum()
    f_disktocor = omegawt/np.pi  # Fraction of disk photons from this radius which hit the corona
    f_cortodisk = f_disktocor/((2.*np.pi*r_cor*h_cor)+(np.pi*np.square(r_cor))) # Fraction of coronal photons
    # hitting a unit area of disk surface at this radius
    return omega, f_disktocor, f_cortodisk

def inv_cone(r,params):
    """ Calculates (numerically) solid angle, illumination fractions for a (frustrum of) an inverted conical
    'solid' corona with base radius r_cor, height h_cor, radius at the top r_top. The model finds solid angle
    and illumination fractions by numerical integration using the approach given in Appendix in
    the UM20 paper. Number of cell spacings on vertical and azimuthal axes given by nphi and nz respectively
    (params[3] and params[4])"""
    r_cor = params[0]
    h_cor = params[1]
    r_top = params[2]
    nphi = int(params[3])
    nz = int(params[4])
    dz = h_cor/nz
    dphi = 2.*np.pi/nphi
    z = np.linspace(dz/2,h_cor-dz/2,nz)
    phi = np.linspace(dphi/2,2*np.pi-(dphi/2),nphi)
    phi_arr, z_arr = np.meshgrid(phi, z)
    cone_angle = np.arctan((r_top-r_cor)/h_cor)
    r_cone = r_cor + z_arr*(r_top-r_cor)/h_cor
    da = dphi*r_cone*dz/np.cos(cone_angle)
    x1 = np.multiply(np.cos(phi_arr),r_cone)
    y1 = np.multiply(np.sin(phi_arr),r_cone)
    z1 = -1.*r_cone*(r_top-r_cor)/h_cor
    xpos = r-x1
    ypos = -1.*y1
    zpos = -1.*z_arr
    rpos = np.sqrt(np.square(xpos)+np.square(ypos)+np.square(zpos))
    xpos = xpos/rpos
    ypos = ypos/rpos
    zpos = zpos/rpos
    daomega = np.multiply(da,np.multiply(x1,xpos)+np.multiply(y1,ypos)+np.multiply(z1,zpos))
    rcorpos = np.sqrt(np.square(r_cone)+np.square(z1))
    daomega = np.divide(daomega,(np.square(rpos)*rcorpos))
    daomegawt = np.multiply(daomega,np.abs(zpos))
    omega = daomega[daomega > 0].sum()
    omegawt = daomegawt[daomegawt > 0].sum()
    f_disktocor = omegawt/np.pi  # Fraction of disk photons from this radius which hit the corona
    f_cortodisk = f_disktocor/(np.pi*(r_top**2+(r_cor+r_top)*np.sqrt((r_top-r_cor)**2+h_cor**2))) # Fraction of
    # coronal photons hitting a unit area of disk surface at this radius
    return omega, f_disktocor, f_cortodisk
    


def calc_timing_params(rad,i_rsigmax,rcor,i_rcor,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par):
    """ Function to calculate variability time-scale for each radial bin, and from that assign a signal 
    power spectrum to each bin (defining the Lorentzian frequency and assigning the appropriate q and rms).
    There are 3 models for assigning Lorentzians to radii: continuous (Lorentzian signals at all radii 
    within the maximum signal radius), multi-frequency (multiple Lorentzian parameters already given, and 
    assigned to radii according to the radial dependence of tau), multi-radius (radii pre-specified 
    instead of frequencies and time-scales and hence frequency then set according to the radial 
    tau dependence). 

    Models and parameters examples:

    lor_model = 'continuous' # Driving signals are split across all radial bins to rsigmax (equivalent to basic Lyubarskii 
                         # model with mdot fluctuations at all radii)
    lor_par = [1.0,0.3,0.0] # For continuous model: [q,total_rms_in_disk,total_rms_in_corona], if the last parameter is 
                        # < 0, the total_rms_in_disk is split evenly across corona and disk radii
                        # If the 2nd parameter is < 0, it corresponds to the rms per radial bin of the disk, and same for
                        # corona if 2nd parameter also negative

    lor_model = 'multi_frequency'  # Driving signals at discrete (observed units) frequencies, mapped to radii using 
                                # the tau parameters
    lor_par = [1e-5, 0.5, 0.2, 2e-6, 0.5, 0.1] # [peak frequency_1, q_1, rms_1, peak_frequency_2, q_2,...] can be any 
                                            # 3*n (n is an integer) in length

    lor_model = 'multi_radius' # Driving signals in discrete radial bins, with frequency worked out according to the tau
                            # parameters.
    lor_par = [25.0,0.5,0.1]  # [radius_1, q_1, rms_1, radius_2, q_2,...] can be any 3*n (n is an integer) in length

    lor_model = 'multi_radius_frequency'  # The frequencies for chosen radii are selected, so there is no link to the
    model for variability time-scales.
    lor_par = [20.0,0.3,0.5,0.2,15.0,2.6,0.5,0.15,11.0,8.4,0.5,0.15] # [radius_1, peak frequency_1, q_1, rms_1, radius_2, 
    peak_frequency_2, q_2,...] can be any 4*n (n is an integer) in length

    """
    
    tau = np.zeros(len(rad))
    lfreq = np.zeros(len(rad))
    q = np.zeros(len(rad))
    rms = np.zeros(len(rad))
    
    # Assign time-scales to radii depending on whether the radius is in the disk or corona
    # Normalisation is defined as value at 10 r_g (whether or not disk or corona reach 10 r_g)
    for i in range(0, len(rad)):
        if (rad[i] > rcor):
            tau[i] = (disk_tau_par[0]*(rad[i]/rcor)**disk_tau_par[1])*rad[i]**1.5
        else:
            tau[i] = (cor_tau_par[0]*(rad[i]/rcor)**cor_tau_par[1])*rad[i]**1.5

    # The continuous model assigns a Lorentzian to every radial bin inside (and including) the maximum signal radius.
    # The rms is calculated by splitting the rms^2 equally between all bins, so produce the required total rms.
    # Parameter lor_par[2] sets the rms of the coronal part of the flow - if it is negative the coronal radii are set
    # to the same rms as the disk (the total rms^2 is split equally between radial bins including the disk and 
    # coronal bins).  A static corona (or disk) can be set by making the rms equal zero.
    if (lor_model == 'continuous'):
        for i in range(0, len(rad)):
            lfreq[i] = 1./(tau[i]*t_scale)
            q[i] = lor_par[0]
            if (i > i_rsigmax):
                rms[i] = 0.
            else:
                if (lor_par[1] >= 0. and lor_par[2] >= 0.):
                    if (rad[i] > rcor):
                        rms[i] = np.sqrt(lor_par[1]**2/((i_rsigmax+1)-i_rcor))
                    else:
                        rms[i] = np.sqrt(lor_par[2]**2/i_rcor)                

                elif (lor_par[1] >= 0. and lor_par[2] < 0.):
                    rms[i] = np.sqrt(lor_par[1]**2/(i_rsigmax+1))

                elif (lor_par[1] < 0. and lor_par[2] >= 0.):
                    if (rad[i] > rcor):
                        rms[i] = -1.*lor_par[1]
                    else:
                        rms[i] = np.sqrt(lor_par[2]**2/i_rcor)

                else:
                    if (rad[i] > rcor):
                        rms[i] =  -1.*lor_par[1]
                    else:
                        rms[i] =  -1.*lor_par[2]
            
    elif (lor_model == 'multi_frequency'):
        for j in range(0, len(lor_par), 3):
            ltau = (1./lor_par[j])/t_scale
            if (ltau <= tau[len(rad)-1] and ltau >= tau[0]):
                ltau2, i = find_nearest(tau,ltau)
                lfreq[i] = 1./(tau[i]*t_scale)
                q[i] = lor_par[j+1]
                if (i > i_rsigmax):
                    rms[i] = 0.
                else:                
                    rms[i] = lor_par[j+2]
                
    elif (lor_model == 'multi_radius'):
        for j in range(0, len(lor_par), 3):
            if (lor_par[j] <= rad[len(rad)-1] and lor_par[j] >= rad[0]):
                lrad, i = find_nearest(rad,lor_par[j])
                lfreq[i] = 1./(tau[i]*t_scale)
                q[i] = lor_par[j+1]
                if (i > i_rsigmax):
                    rms[i] = 0.
                else:
                    rms[i] = lor_par[j+2]        
                    
    elif (lor_model == 'multi_radius_frequency'):
        for j in range(0, len(lor_par), 4):
            if (lor_par[j] <= rad[len(rad)-1] and lor_par[j] >= rad[0]):
                lrad, i = find_nearest(rad,lor_par[j])
                if (rms[i] > 0): # In case this radius already has a Lorentzian assigned, go to next available radius:
                    for i_next in range(1,i,1):
                        if (rms[i-i_next] == 0. and rms[i] > 0):
                            i = i-i_next
                lfreq[i] = lor_par[j+1]
                q[i] = lor_par[j+2]
                if (i > i_rsigmax):
                    rms[i] = 0.
                else:
                    rms[i] = lor_par[j+3]        
                print("Signal radius for ",lor_par[j+1]," Hz is ",rad[i])
                print("with rms ",rms[i])

    return tau, lfreq, q, rms


def calc_propagation_params(rad,rad_edge,rcor,disk_prop_par,cor_prop_par):
    """ Calculates the propagation delay across each bin, using either disk or coronal parameters depending on
    whether the radius is inside or outside the coronal radius rcor. """
    deltau = np.zeros(len(rad))
    delrad = np.diff(rad_edge)
    deltau[rad <= rcor] = cor_prop_par[0]*delrad[rad <= rcor]*np.sqrt(rad[rad <= rcor])*(rad[rad <= rcor]/rcor)**cor_prop_par[1]
    deltau[rad > rcor] = disk_prop_par[0]*delrad[rad > rcor]*np.sqrt(rad[rad > rcor])*(rad[rad > rcor]/rcor)**disk_prop_par[1]
    return deltau            
    
    
def calc_radial_time_response(rad,rcor,disp_frac,disktocor_frac,cortodisk_frac,thermal_frac,seed_frac_flow,heat_frac_flow):
    """ Calculates the direct disk, seed and heating impulse responses due to the dissipation generated in the 
    radial bins.  Treats the direct disk and seed components due to dissipation and reverberation separately.
    Note that the output impulse responses are integrated values over each radial delay range - for the final
    IRFs used for calculation these will be divided by the time-delay ranges for each radial bin to give
    the response per unit time-delay. """

    # First set up the arrays for the responses:
    ldisk_disp = np.zeros(len(rad))

    # Disk emission is simply the viscous dissipation.
    # To conserve luminosity we should subtract the fraction entering the corona from the disk luminosity, the resulting
    # disk luminosity is therefore the angle-averaged value seen by the observer (not the intrinsic value).
    ldisk_disp[rad > rcor] = disp_frac[rad > rcor]*(1-disktocor_frac[rad > rcor])
    
    # Within rcor the dissipated energy associated with the corona is split between internal coronal seed photons
    # and coronal heating. We need add the seed component to anything which may be present from the inner disk.
    # We must also correct for the fraction which is reprocessed and returned by the disk.
    f_rev = np.sum(cortodisk_frac) # The fraction of coronal luminosity that is reprocessed by the disk
    f_return = np.sum(cortodisk_frac*disktocor_frac) # the fraction of coronal luminosity that is reprocessed by the disk
    # and returned to the corona
    lseed_disp = disp_frac*(disktocor_frac+seed_frac_flow)
    lheat = disp_frac*heat_frac_flow
    # The reverberation seed is just the luminosity returning to the corona
    lseed_rev = thermal_frac*f_return*(lseed_disp+lheat)/(1-f_return)
    # Finally the reverberation which escapes the system without returning to the corona
    ldisk_rev = thermal_frac*(f_rev-f_return)*(lseed_disp+lheat)/(1-f_return)
    
    return ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev

def calc_irfs_mono(gamma_par,e_seed,energy,ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev):
    """ Calculates the IRFs of the energy dependent photon flux (for single, i.e. 'mono' energies), photon index, 
    the total disk emission (dissipation + reverberation) and the total seed emission, using the disk, seed and 
    heating impulse responses as input. 
    Disk and seed emission IRFs are normalised to the total disk and seed emission.
    Again, note that the IRF values used here are integrated over the delay range of each radial bin - they can
    be divided by the delay ranges later, to give correct IRF 'per unit time-delay'.  This also works for the
    flux and photon index IRFs, since the relevant equations (7 & 8) are the same whether or not we multiply
    both sides by a delay increment delta-tau. """
        
    lheat_sum = np.sum(lheat)
    lseed_sum = np.sum(lseed_disp+lseed_rev)
    gamma_mean = gamma_par[0]*(lseed_sum/lheat_sum)**gamma_par[1]
    gamma_irf = gamma_mean*gamma_par[1]*((lseed_disp + lseed_rev)/lseed_sum - lheat/lheat_sum)    
    u = gamma_mean*gamma_par[1]
    v = (gamma_mean-1.)
    w = np.reshape(np.log(energy/e_seed),(len(energy),1))
    n = lseed_disp.size
    flux_irf = (1-u*w+(u/v))*np.reshape(lseed_disp+lseed_rev,(1,n))/lseed_sum + (u*w-(u/v))*np.reshape(lheat,(1,n))/lheat_sum
    disk_irf = ldisk_disp + ldisk_rev
    seed_irf = lseed_disp + lseed_rev

    return gamma_mean, gamma_irf, flux_irf, disk_irf, seed_irf
    

def calc_disk_band(disp_frac,cortodisk_frac,disktocor_frac,lseed_disp,lheat,thermal_frac,rcor,rad,rad_area,
                   eband,kTvar):
    # Given the fractional of total flux due to
    # dissipation, and the area
    # of each radial annulus of the disk, calculate the fraction of the flux at each radius which is emitted in a fixed
    # energy band, [E_min,E_max].  The energy band is defined not in absolute energy but relative to the maximum
    # blackbody temperature in the disk.  E.g. we observe an inner disk temperature of 0.3 keV and assume to probe
    # this to correspond to the maximum.  If we then want to probe an energy range 0.6-0.9 keV, we should specify
    # eband = [2.0,3.0]
    
    f_rev = np.sum(cortodisk_frac) # The fraction of coronal luminosity that is reprocessed by the disk
    f_return = np.sum(cortodisk_frac*disktocor_frac) # the fraction of coronal luminosity that is reprocessed by the disk
    radldisk_rev = cortodisk_frac*thermal_frac*np.sum((lseed_disp+lheat)/(1-f_return))
    
    kT_rad = np.zeros(len(disp_frac))
    band_frac = np.zeros(len(disp_frac))
    kT_rad[rad > rcor] = ((disp_frac[rad > rcor]+
                                 radldisk_rev[rad > rcor])/rad_area[rad > rcor])**0.25
    kT_max = np.amax(kT_rad)
    kT_rad = kT_rad/kT_max
    
    E_min = eband[0]
    E_max = eband[1]
    # Total ph flux of bb_phlux for 1 keV = 2*Apery's constant *15/pi^4 for non-variable bb,
    # 1.5*Apery's constant *15/pi^4 for variable bb
    if kTvar == False:
        bb_phintnorm = 2*1.2020569031595942*15/np.pi**4
    else:
        bb_phintnorm = 1.5*1.2020569031595942*15/np.pi**4
    for i, radval in enumerate(rad):
        if radval > rcor:
            band_frac[i] = bb_phflux(E_min, E_max, kT_rad[i], 1000, kTvar)/(bb_phintnorm/kT_rad[i])
        
    ldiskband_disp = disp_frac*(1-disktocor_frac)*band_frac
    lseedband_disp = disp_frac*disktocor_frac*band_frac
    
    ldiskband_rev_arr = thermal_frac*band_frac*(cortodisk_frac*(1-
                    disktocor_frac))*(lseed_disp.reshape(len(rad),1)+lheat.reshape(len(rad),1))/(1-f_return)
    
    lseedband_rev_arr = thermal_frac*band_frac*(cortodisk_frac*disktocor_frac)*\
        (lseed_disp.reshape(len(rad),1)+lheat.reshape(len(rad),1))/(1-f_return)
    ldiskband_rev = np.sum(ldiskband_rev_arr,axis=1)
    lseedband_rev = np.sum(lseedband_rev_arr,axis=1)
    ldisk_band = ldiskband_disp + ldiskband_rev
    lseed_band = lseedband_disp + lseedband_rev
    
    return band_frac, kT_rad, ldisk_band, lseed_band, ldiskband_disp, ldiskband_rev

        
def bb_phflux(E_min, E_max, kT, nE, kTvar): # Calculates blackbody photon flux for blackbody with temperature kT and
    # integrated over an energy range E_min to E_max using nE discrete bins
    # The model is normalised to maintain a constant total energy flux regardless of kT (based on XSPEC bbody model)

    E_edge = np.logspace(np.log10(E_min), np.log10(E_max), nE+1)
    E = np.sqrt(np.multiply(E_edge[1:],E_edge[:-1]))
    dE = E_edge[1:]-E_edge[:-1]
    if kTvar == False:
        bb_flux = np.sum(np.nan_to_num(dE*(E**2/kT**4)/(np.exp(E/kT)-1.)))
    else:
        bb_flux = np.sum(np.nan_to_num(dE*0.25*(E**3/kT**5)*np.exp(E/kT)*(np.exp(E/kT)-1)**-2))

    return bb_flux/((np.pi**4)/15)

def bb_enflux(E_min, E_max, kT, nE, kTvar): # Calculates blackbody energy flux for blackbody with temperature kT and
    # integrated over an energy range E_min to E_max using nE discrete bins
    # The model is normalised to maintain a constant total energy flux regardless of kT (based on XSPEC bbody model)

    E_edge = np.logspace(np.log10(E_min), np.log10(E_max), nE+1)
    E = np.sqrt(np.multiply(E_edge[1:],E_edge[:-1]))
    dE = E_edge[1:]-E_edge[:-1]
    if kTvar == False:
        bb_flux = np.sum(np.nan_to_num(dE*(E**3/kT**4)/(np.exp(E/kT)-1.)))
    else:
        bb_flux = np.sum(np.nan_to_num(dE*(0.25*(E**4/kT**5)*np.exp(E/kT)*(np.exp(E/kT)-1)**-2)))

    return bb_flux/((np.pi**4)/15)


def linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,input_irf,deltau_scale,nirf):
    """ Rebins the given input irf to linear rebinning where the uniform binsize is dt.
    The input IRFs, which are integrated over delay bins, are converted to be per unit time-delay.
    Values are interpolated at the edges of the radial delay bins which do not fit exactly into an 
    integer multiple of dt.

    Note that the ouput IRFs only extend to the maximum signal radius - this is to save memory in case the number
    of radii considered is large but the driving signals start from relatively small radii (because at this point
    the IRF arrays become linear and hence large for a broad range of delays...).

    Note that in this function we also apply an efficiency trick in the case where the innermost coronal 
    propagation delays are set to zero, which is a good approximation when the coronal delays are small compared
    to the delays in the inner disk (as would be expected for a geometrically thick corona and geometrically thin 
    disk).  The function bins all the IRF bins for which deltau_scale = 0 into one single bin. """
    
    rebinned_irf = np.zeros(nirf)
    # We create a new array for the inputs which are converted to be per unit time-delay:
    input_irf2 = np.zeros(len(input_irf))
    for i in range(i_rsigmax,-1,-1):
        if (deltau_scale[i] > 0):
            input_irf2[i] = input_irf[i]/deltau_scale[i]

    for i in range(i_rsigmax,-1,-1):
        if (deltau_scale[i] > 0):
            irfbin_start = int(np.round(np.sum(irf_nbins[i:i_rsigmax+1]) - irf_nbins[i]))
            irfbin_stop = int(np.round((irfbin_start + irf_nbins[i])-1))
#            rebinned_irf[irfbin_start:irfbin_stop+1] = input_irf[i]/(dt*irf_nbins[i])
            rebinned_irf[irfbin_start:irfbin_stop+1] = input_irf2[i]
            if (i < i_rsigmax and irf_binedgefrac[i] > 0):
                rebinned_irf[irfbin_start] = input_irf2[i]+(irf_binedgefrac[i]*(input_irf2[i+1]-input_irf2[i]))
            if (i > 0 and irf_binedgefrac[i-1] <= 0):
                rebinned_irf[irfbin_stop] = input_irf2[i]+(irf_binedgefrac[i-1]*(input_irf2[i]-input_irf2[i-1]))
        else:
            rebinned_irf[irfbin_stop+1] = rebinned_irf[irfbin_stop+1] + input_irf[i]

    rebinned_irf[irfbin_stop+1] = rebinned_irf[irfbin_stop+1]/dt # Will be zero if deltau_scale[0] > 0., otherwise
    # this just outputs the combined IRF of coronal deltau_scale = 0 into a single bin (normalised assuming the 
    # width of the bin, which is dt)
    
    # Now calculate flux component outside rsigmax, needed for constant flux contribution
    if (i_rsigmax < len(input_irf)-1):
        flux_outer = np.sum(input_irf[i_rsigmax+1:])
    else:
        flux_outer = 0.
    
    return rebinned_irf, flux_outer

def calc_cross_psd(freq,dt,ci_irf,ref_irf,irf_nbins,deltau_scale,i_rsigmax,f_pk,q,rms):
    """ Produces cross/power-spectra. Takes as input a pair of linear-binned IRFs, the radial binning 
    information and the Lorentzian signal parameters per radial bin, then produces the combined weighted 
    power spectra and cross-spectra by summing for each radial bin the calculated power spectra and 
    cross-spectra, weighted by the driving signal PSD for that bin. """
    
    nirf = len(ref_irf)
    wt_cross_spec = np.zeros(nirf//2)
    wt_pow_spec_ref = np.zeros(nirf//2)
    wt_pow_spec_ci = np.zeros(nirf//2)
    mod_sig_psd = np.zeros(nirf//2)
    ref_irf_dum = np.multiply(ref_irf,dt)
    ci_irf_dum = np.multiply(ci_irf,dt)
    for i in range(i_rsigmax,-1,-1):
        if (i == i_rsigmax or deltau_scale[i+1] > 0.):
            if (i < i_rsigmax):
                irfbin_start = int(np.round(np.sum(irf_nbins[i+1:i_rsigmax+1]) - irf_nbins[i+1]))
                irfbin_stop = int(np.round((irfbin_start + irf_nbins[i+1])-1))
                ref_irf_dum[irfbin_start:irfbin_stop+1] = 0.
                ci_irf_dum[irfbin_start:irfbin_stop+1] = 0.
            if (rms[i] > 0.):
                lor_psd = lorentz_q(freq, f_pk[i], q[i], rms[i])
                ref_fft = spfft.fft(ref_irf_dum)[1:nirf//2+1]
                ci_fft = spfft.fft(ci_irf_dum)[1:nirf//2+1]
                cross_spec = np.multiply(np.conj(ci_fft),ref_fft)
                wt_cross_spec = wt_cross_spec + np.multiply(lor_psd[0,:],cross_spec)
                pow_spec_ref = np.square(np.abs(ref_fft))
                pow_spec_ci = np.square(np.abs(ci_fft))
                wt_pow_spec_ref = wt_pow_spec_ref + np.multiply(lor_psd[0,:],pow_spec_ref)
                wt_pow_spec_ci = wt_pow_spec_ci + np.multiply(lor_psd[0,:],pow_spec_ci)
                mod_sig_psd = mod_sig_psd + lor_psd                
    return wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd

def lorentz_qold(f, f_pk, q, rms):  
    """ Lorentzian function defined in terms of peak frequency and quality factor q
    e.g. see Pottschmidt et al. 2003, A&A, 407, 1039 for more info """
    
    f_res=f_pk/np.sqrt(1.0+(1.0/(4.0*q**2)))
    r=rms/np.sqrt(0.5-np.arctan(-2.0*q)/np.pi)
    lorentz=((1/np.pi)*2*r**2*q*f_res)/ \
    (f_res**2+(4*q**2*np.square(f-f_res)))
    return lorentz

def lorentz_q(f, f_pk, q, rms):  
    """ Returns array of Lorentzian functions, each defined in terms of peak frequency and quality factor q
    e.g. see Pottschmidt et al. 2003, A&A, 407, 1039 for more info """
    
    f_res=np.divide(f_pk,np.sqrt(1.0+(1.0/(4.0*np.square(q)))))
    r=np.divide(rms,np.sqrt(0.5-np.arctan(-2.0*q)/np.pi))
    lorentz1=(1/np.pi)*2*np.multiply(np.power(r,2),np.multiply(q,f_res))
    lorentz2 = 4*np.multiply(np.square(q),np.square(np.subtract(np.reshape(f,[len(f),1]),f_res)))
    lorentz = np.divide(lorentz1,np.square(f_res)+lorentz2)
    return np.transpose(lorentz)

def calculate_stprod_mono(nirf_mult,energy,encomb,flux_irf,disk_irf,gamma_irf,deltau,min_deltau_frac,i_rsigmax,\
                        lfreq,q,rms,t_scale):
    """ Calculates lag-frequency and power spectra for disk, seed and the power-law emission from different 
    bands, based on simplified model and not full-energy-dependent version (which uses calculate_stprod_endep).

    This function takes as input the energy-dependent photon flux, disk and photon index IRFs obtained for the
    original geometrically spaced radial delay bins, as well as the Lorentzian parameters per radius
    and the delay bin width per radius deltau.  It then does the following:
    
    1. Rescales the deltau to the observational time units (probably s) using t_scale and sets the time bin size dt 
    to be the linear rebinned IRFs to be some fraction min_deltau_frac of the minimum delay bin width, combines 
    this with the total delay time from the maximum signal radius and a multiplication factor to get a total 
    length for the linear-rebinned IRFs.
    
    2. Using the bin width dt, assigns the time-delays of the geometrically-spaced radial bins to integer numbers
    of bins of width dt, for use in calculating the linear-rebinned IRFs, also works out overlaps at the bin edges
    for interpolation.
     
    3. Proceeds to linearly-rebin the IRFs and then use these with calculated Fourier frequencies and the radial
    Lorentzian parameters to work out the energy/band-dependent power-spectra and cross-spectra and hence lags 
    for the combinations of disk/PL-energies determined by the encomb array: e.g. [[1,0],[2,1],[3,1]] will compare:
    PL-energy[0] vs. disk, PL-energy[1] vs. PL-energy[0] and PL-energy[2] vs PL-energy[0] """


    print("#######################################")
#    print("Calculating mono-energetic spectral-timing products")
    
    deltau_scale = deltau*t_scale  # Rescale to actual time units used
    dt = min_deltau_frac*np.amin(deltau_scale[deltau_scale > 0])  # The time bin size of the irfs is a fraction 
    # min_deltau_frac of the smallest non-zero value of the propagation time across a bin
    
    nirf = int(nirf_mult*np.power(2,np.math.ceil(np.log2(np.sum(deltau_scale[:i_rsigmax+1])/dt)))) # Next highest 
    # integer powers of 2 that encompasses the full IRF, for the IRF used for spectral-timing with FFT.
    print("Time bin size dt is: ",dt)
    print("The maximum propagation delay is: ",np.sum(deltau_scale[:i_rsigmax+1])," and there are ",nirf," irf bins.")
    irf_nbins = np.zeros(len(deltau_scale))    
    irf_binedgefrac = np.zeros(len(deltau_scale))
    deltau_sum_max = np.sum(deltau_scale[:i_rsigmax+1])
    for i in range(0,i_rsigmax+1):
        i_irf_max = int(np.round((deltau_sum_max - np.sum(deltau_scale[:i]))/dt))-1
        i_irf_min = int(np.round((deltau_sum_max - np.sum(deltau_scale[:i+1]))/dt))
        irf_nbins[i] = (i_irf_max-i_irf_min)+1
        irf_binedgefrac[i] = (dt*np.sum(irf_nbins[:i+1])-np.sum(deltau_scale[:i+1]))/dt

    ci_irf = np.zeros((len(energy)+1,nirf))
    ci_outer = np.zeros(len(energy)+1)
    ci_mean = np.zeros(len(energy)+1) 
    for i in range(len(energy)+1):
        if i == 0:  # For 0 use the disk irf
            ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,disk_irf,deltau_scale,nirf)
            ci_mean[i] = dt*np.sum(ci_irf[i,:])+ci_outer[i]
        else:            
            ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,flux_irf[(i-1),:],
                                                deltau_scale,nirf)
            ci_mean[i] = 1.0 # The PL IRF is already normalised by the mean PL flux at each energy
        
    minfreq = 1./(dt*nirf)
    freq = np.linspace(minfreq, nirf*minfreq/2, nirf//2)

    phlag = np.zeros((len(encomb),len(freq)))
    tlag = np.zeros((len(encomb),len(freq)))
    psd_ci = np.zeros((len(encomb),len(freq)))
    psd_ref = np.zeros((len(encomb),len(freq)))

    for i in range(len(encomb)): # note that encomb values j,k = 0 correspond to the disk, values j, k > 0 correspond to
      #  an index (j,k) -1 in the 'energy' array of monochromatic energies for PL fluxes
        j = encomb[i,0]
        k = encomb[i,1]
#        print("CI mean, ref mean, CI outer, ref outer : ", ci_mean[j], ci_mean[k], ci_outer[j], ci_outer[k])
              
        wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd = calc_cross_psd(freq,dt,ci_irf[j,:],ci_irf[k,:],\
                    irf_nbins,deltau_scale,i_rsigmax,lfreq,q,rms)

        phlag[i,:] = np.angle(wt_cross_spec)
        tlag[i,:] = np.divide(phlag[i,:],2*np.pi*freq)
        psd_ci[i,:] = wt_pow_spec_ci/np.square(ci_mean[j])
        psd_ref[i,:] = wt_pow_spec_ref/np.square(ci_mean[k])

#        print("Calculated for energies ",encomb[i,:])
        
#    print("#######################################")
    
    return freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf,\
        ci_irf, ci_mean, ci_outer
        
        
def calculate_stprod_mono_frange(nirf_mult,energy,encomb,flux_irf,disk_irf,gamma_irf,deltau,min_deltau_frac,i_rsigmax,\
                        lfreq,q,rms,t_scale,freq_range):
    """ Calculates lag-frequency and power spectra for disk, seed and the power-law emission from different
    bands, based on simplified model and not full-energy-dependent version (which uses calculate_stprod_endep).

    This function takes as input the energy-dependent photon flux, disk and photon index IRFs obtained for the
    original geometrically spaced radial delay bins, as well as the Lorentzian parameters per radius
    and the delay bin width per radius deltau.  It then does the following:
    
    1. Rescales the deltau to the observational time units (probably s) using t_scale and sets the time bin size dt
    to be the linear rebinned IRFs to be some fraction min_deltau_frac of the minimum delay bin width, combines
    this with the total delay time from the maximum signal radius and a multiplication factor to get a total
    length for the linear-rebinned IRFs.
    
    2. Using the bin width dt, assigns the time-delays of the geometrically-spaced radial bins to integer numbers
    of bins of width dt, for use in calculating the linear-rebinned IRFs, also works out overlaps at the bin edges
    for interpolation.
     
    3. Proceeds to linearly-rebin the IRFs and then use these with calculated Fourier frequencies and the radial
    Lorentzian parameters to work out the energy/band-dependent power-spectra and cross-spectra and hence lags
    for the combinations of disk/PL-energies determined by the encomb array: e.g. [[1,0],[2,1],[3,1]] will compare:
    PL-energy[0] vs. disk, PL-energy[1] vs. PL-energy[0] and PL-energy[2] vs PL-energy[0] """


    print("#######################################")
#    print("Calculating mono-energetic spectral-timing products")
    
    deltau_scale = deltau*t_scale  # Rescale to actual time units used
    
    if freq_range == None:
        dt = min_deltau_frac*np.amin(deltau_scale[deltau_scale > 0])  # The time bin size of the irfs is a fraction
        # min_deltau_frac of the smallest non-zero value of the propagation time across a bin
        nirf = int(nirf_mult*np.power(2,np.math.ceil(np.log2(np.sum(deltau_scale[:i_rsigmax+1])/dt)))) # Next highest
        # integer powers of 2 that encompasses the full IRF, for the IRF used for spectral-timing with FFT.
    else:
        dt = min(min_deltau_frac*np.amin(deltau_scale[deltau_scale > 0]), min_deltau_frac/(2*freq_range[1]))
        nirf1 = int(nirf_mult*np.power(2,np.math.ceil(np.log2(np.sum(deltau_scale[:i_rsigmax+1])/dt))))
        nirf2 = int(nirf_mult*np.power(2,np.math.ceil(np.log2((1/freq_range[0])/dt))))
        nirf = max(nirf1,nirf2)
        
    print("Time bin size dt is: ",dt)
    print("The maximum propagation delay is: ",np.sum(deltau_scale[:i_rsigmax+1])," and there are ",nirf," irf bins.")
    irf_nbins = np.zeros(len(deltau_scale))
    irf_binedgefrac = np.zeros(len(deltau_scale))
    deltau_sum_max = np.sum(deltau_scale[:i_rsigmax+1])
    for i in range(0,i_rsigmax+1):
        i_irf_max = int((deltau_sum_max - np.sum(deltau_scale[:i]))/dt)-1
        i_irf_min = int((deltau_sum_max - np.sum(deltau_scale[:i+1]))/dt)
        irf_nbins[i] = (i_irf_max-i_irf_min)+1
        irf_binedgefrac[i] = (dt*np.sum(irf_nbins[:i+1])-np.sum(deltau_scale[:i+1]))/dt

    ci_irf = np.zeros((len(energy)+1,nirf))
    ci_outer = np.zeros(len(energy)+1)
    ci_mean = np.zeros(len(energy)+1)
    for i in range(len(energy)+1):
        if i == 0:  # For 0 use the disk irf
            ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,disk_irf,deltau_scale,nirf)
            ci_mean[i] = dt*np.sum(ci_irf[i,:])+ci_outer[i]
        else:
            ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,flux_irf[(i-1),:],
                                                deltau_scale,nirf)
            ci_mean[i] = 1.0 # The PL IRF is already normalised by the mean PL flux at each energy
        
    minfreq = 1./(dt*nirf)
    freq = np.linspace(minfreq, nirf*minfreq/2, nirf//2)

    phlag = np.zeros((len(encomb),len(freq)))
    tlag = np.zeros((len(encomb),len(freq)))
    psd_ci = np.zeros((len(encomb),len(freq)))
    psd_ref = np.zeros((len(encomb),len(freq)))

    for i in range(len(encomb)): # note that encomb values j,k = 0 correspond to the disk, values j, k > 0 correspond to
      #  an index (j,k) -1 in the 'energy' array of monochromatic energies for PL fluxes
        j = encomb[i,0]
        k = encomb[i,1]
#        print("CI mean, ref mean, CI outer, ref outer : ", ci_mean[j], ci_mean[k], ci_outer[j], ci_outer[k])
              
        wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd = calc_cross_psd(freq,dt,ci_irf[j,:],ci_irf[k,:],\
                    irf_nbins,deltau_scale,i_rsigmax,lfreq,q,rms)

        phlag[i,:] = np.angle(wt_cross_spec)
        tlag[i,:] = np.divide(phlag[i,:],2*np.pi*freq)
        psd_ci[i,:] = wt_pow_spec_ci/np.square(ci_mean[j])
        psd_ref[i,:] = wt_pow_spec_ref/np.square(ci_mean[k])

#        print("Calculated for energies ",encomb[i,:])
        
#    print("#######################################")
    
    return freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf,\
        ci_irf, ci_mean, ci_outer


def calc_irfs_mono_supp(gamma_par,e_seed,energy,ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev,var_suppress_params,rad):
    """ Calculates the IRFs of the energy dependent photon flux (for single, i.e. 'mono' energies), photon index,
    the total disk emission (dissipation + reverberation) and the total seed emission, using the disk, seed and
    heating impulse responses as input.
    Disk and seed emission IRFs are normalised to the total disk and seed emission.
    Again, note that the IRF values used here are integrated over the delay range of each radial bin - they can
    be divided by the delay ranges later, to give correct IRF 'per unit time-delay'.  This also works for the
    flux and photon index IRFs, since the relevant equations (7 & 8) are the same whether or not we multiply
    both sides by a delay increment delta-tau. """
        
    lheat_sum = np.sum(lheat)
    lseed_sum = np.sum(lseed_disp+lseed_rev)
    gamma_mean = gamma_par[0]*(lseed_sum/lheat_sum)**gamma_par[1]
    irfs_list = []
    if var_suppress_params != None:
        for params in var_suppress_params:
            [rad_supp,fdisk_supp,fseed_supp,fheat_supp] = params
            lheat2 = np.where(rad > rad_supp, fheat_supp*lheat, lheat)
            ldisk_disp2 = np.where(rad > rad_supp, fdisk_supp*ldisk_disp, ldisk_disp)
            lseed_disp2 = np.where(rad > rad_supp, fseed_supp*lseed_disp, lseed_disp)
            ldisk_rev2 = np.where(rad > rad_supp, ldisk_rev*(fseed_supp*lseed_disp+fheat_supp*lheat)/(lheat+lseed_disp), ldisk_rev)
            lseed_rev2 = np.where(rad > rad_supp, lseed_rev*(fseed_supp*lseed_disp+fheat_supp*lheat)/(lheat+lseed_disp), lseed_rev)
            gamma_irf = gamma_mean*gamma_par[1]*((lseed_disp2 + lseed_rev2)/lseed_sum - lheat2/lheat_sum)
            u = gamma_mean*gamma_par[1]
            v = (gamma_mean-1.)
            w = np.reshape(np.log(energy/e_seed),(len(energy),1))
            n = lseed_disp.size
            flux_irf = (1-u*w+(u/v))*np.reshape(lseed_disp2+lseed_rev2,(1,n))/lseed_sum + (u*w-(u/v))*np.reshape(lheat2,(1,n))/lheat_sum
            disk_irf = ldisk_disp2 + ldisk_rev2
            seed_irf = lseed_disp2 + lseed_rev2
            irfs_list.append([disk_irf,flux_irf,gamma_irf,seed_irf,lheat2])
    else:
        gamma_irf = gamma_mean*gamma_par[1]*((lseed_disp + lseed_rev)/lseed_sum - lheat/lheat_sum)
        u = gamma_mean*gamma_par[1]
        v = (gamma_mean-1.)
        w = np.reshape(np.log(energy/e_seed),(len(energy),1))
        n = lseed_disp.size
        flux_irf = (1-u*w+(u/v))*np.reshape(lseed_disp+lseed_rev,(1,n))/lseed_sum + (u*w-(u/v))*np.reshape(lheat,(1,n))/lheat_sum
        disk_irf = ldisk_disp + ldisk_rev
        seed_irf = lseed_disp + lseed_rev
        irfs_list.append([disk_irf,flux_irf,gamma_irf,seed_irf,lheat])
    return gamma_mean, irfs_list


def calculate_stprod_mono_supp(nirf_mult,energy,encomb,irfs_list,var_suppress_params,rad,deltau,min_deltau_frac,i_rsigmax,\
                        lfreq,q,rms,t_scale):
    """ Calculates lag-frequency and power spectra for disk, seed and the power-law emission from different
    bands, based on simplified model and not full-energy-dependent version (which uses calculate_stprod_endep).

    This function takes as input the energy-dependent photon flux, disk and photon index IRFs obtained for the
    original geometrically spaced radial delay bins, as well as the Lorentzian parameters per radius
    and the delay bin width per radius deltau.  It then does the following:
    
    1. Rescales the deltau to the observational time units (probably s) using t_scale and sets the time bin size dt
    to be the linear rebinned IRFs to be some fraction min_deltau_frac of the minimum delay bin width, combines
    this with the total delay time from the maximum signal radius and a multiplication factor to get a total
    length for the linear-rebinned IRFs.
    
    2. Using the bin width dt, assigns the time-delays of the geometrically-spaced radial bins to integer numbers
    of bins of width dt, for use in calculating the linear-rebinned IRFs, also works out overlaps at the bin edges
    for interpolation.
     
    3. Proceeds to linearly-rebin the IRFs and then use these with calculated Fourier frequencies and the radial
    Lorentzian parameters to work out the energy/band-dependent power-spectra and cross-spectra and hence lags
    for the combinations of disk/PL-energies determined by the encomb array: e.g. [[1,0],[2,1],[3,1]] will compare:
    PL-energy[0] vs. disk, PL-energy[1] vs. PL-energy[0] and PL-energy[2] vs PL-energy[0] """


    print("#######################################")
    print("Calculating mono-energetic spectral-timing products")
    
    deltau_scale = deltau*t_scale  # Rescale to actual time units used
    dt = min_deltau_frac*np.amin(deltau_scale[deltau_scale > 0])  # The time bin size of the irfs is a fraction
    # min_deltau_frac of the smallest non-zero value of the propagation time across a bin
    
    nirf = int(nirf_mult*np.power(2,np.math.ceil(np.log2(np.sum(deltau_scale[:i_rsigmax+1])/dt)))) # Next highest
    # integer powers of 2 that encompasses the full IRF, for the IRF used for spectral-timing with FFT.
    print("Time bin size dt is: ",dt)
    print("The maximum propagation delay is: ",np.sum(deltau_scale[:i_rsigmax+1])," and there are ",nirf," irf bins.")
    irf_nbins = np.zeros(len(deltau_scale))
    irf_binedgefrac = np.zeros(len(deltau_scale))
    deltau_sum_max = np.sum(deltau_scale[:i_rsigmax+1])
    for i in range(0,i_rsigmax+1):
        i_irf_max = int((deltau_sum_max - np.sum(deltau_scale[:i]))/dt)-1
        i_irf_min = int((deltau_sum_max - np.sum(deltau_scale[:i+1]))/dt)
        irf_nbins[i] = (i_irf_max-i_irf_min)+1
        irf_binedgefrac[i] = (dt*np.sum(irf_nbins[:i+1])-np.sum(deltau_scale[:i+1]))/dt

    ci_irfs_list = []
    for j, irfs in enumerate(irfs_list):
        ci_irf = np.zeros((len(energy)+1,nirf))
        ci_outer = np.zeros(len(energy)+1)
        ci_mean = np.zeros(len(energy)+1)
        [disk_irf,flux_irf,gamma_irf,seed_irf,lheat] = irfs
        for i in range(len(energy)+1):
            if i == 0:  # For 0 use the disk irf
                ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,disk_irf,deltau_scale,nirf)
                ci_mean[i] = dt*np.sum(ci_irf[i,:])+ci_outer[i]
            else:
                ci_irf[i,:], ci_outer[i] = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,flux_irf[(i-1),:],
                                                deltau_scale,nirf)
                ci_mean[i] = 1.0 # The PL IRF is already normalised by the mean PL flux at each energy
        ci_irfs_list.append([ci_mean,ci_outer,ci_irf])
        if j == len(irfs_list)-1:
            ci_mean0 = ci_mean
            
                
    minfreq = 1./(dt*nirf)
    freq = np.linspace(minfreq, nirf*minfreq/2, nirf//2)

    phlag = np.zeros((len(encomb),len(freq)))
    tlag = np.zeros((len(encomb),len(freq)))
    psd_ci = np.zeros((len(encomb),len(freq)))
    psd_ref = np.zeros((len(encomb),len(freq)))

    for i in range(len(encomb)): # note that encomb values j,k = 0 correspond to the disk, values j, k > 0 correspond to
      #  an index (j,k) -1 in the 'energy' array of monochromatic energies for PL fluxes
        j = encomb[i,0]
        k = encomb[i,1]
        ci_irf_supp = []
        for nsupp, irfs in enumerate(ci_irfs_list):
            [ci_mean,ci_outer,ci_irf] = irfs
            print("For IRF suppression radius ",nsupp)
            print("CI mean, ref mean, CI outer, ref outer : ", ci_mean[j], ci_mean[k], ci_outer[j], ci_outer[k])
            ci_irf_supp.append([ci_irf[j,:],ci_irf[k,:]])
            
        wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd = calc_cross_psd_supp(freq,dt,ci_irf_supp,var_suppress_params,\
                    rad,irf_nbins,deltau_scale,i_rsigmax,lfreq,q,rms)

        phlag[i,:] = np.angle(wt_cross_spec)
        tlag[i,:] = np.divide(phlag[i,:],2*np.pi*freq)
        psd_ci[i,:] = wt_pow_spec_ci/np.square(ci_mean0[j])
        psd_ref[i,:] = wt_pow_spec_ref/np.square(ci_mean0[k])

        print("Calculated for energies ",encomb[i,:])
        
    print("#######################################")
    
    return freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf,\
        ci_irfs_list


def calc_cross_psd_supp(freq,dt,ci_irf_supp,var_suppress_params,rad,irf_nbins,deltau_scale,i_rsigmax,f_pk,q,rms):
    """ Produces cross/power-spectra. Takes as input a pair of linear-binned IRFs, the radial binning
    information and the Lorentzian signal parameters per radial bin, then produces the combined weighted
    power spectra and cross-spectra by summing for each radial bin the calculated power spectra and
    cross-spectra, weighted by the driving signal PSD for that bin. """
    
    i_supp = 0
    r_supp = var_suppress_params[i_supp][0]
    [ci_irf,ref_irf] = ci_irf_supp[i_supp]
    nirf = len(ref_irf)
    wt_cross_spec = np.zeros(nirf//2)
    wt_pow_spec_ref = np.zeros(nirf//2)
    wt_pow_spec_ci = np.zeros(nirf//2)
    mod_sig_psd = np.zeros(nirf//2)
    ref_irf_dum = np.multiply(ref_irf,dt)
    ci_irf_dum = np.multiply(ci_irf,dt)
    for i in range(i_rsigmax,-1,-1):
        if (i == i_rsigmax or deltau_scale[i+1] > 0.):
            if (i < i_rsigmax):
                if rad[i] < r_supp:
                    i_supp = i_supp+1
                    [ci_irf,ref_irf] = ci_irf_supp[i_supp]
                    ref_irf_dum = np.multiply(ref_irf,dt)
                    ci_irf_dum = np.multiply(ci_irf,dt)
                    r_supp = var_suppress_params[i_supp][0]
                irfbin_start = int(np.sum(irf_nbins[i+1:i_rsigmax+1]) - irf_nbins[i+1])
                irfbin_stop = int((irfbin_start + irf_nbins[i+1])-1)
                ref_irf_dum[:irfbin_stop+1] = 0.
                ci_irf_dum[:irfbin_stop+1] = 0.
            if (rms[i] > 0.):
                lor_psd = lorentz_q(freq, f_pk[i], q[i], rms[i])
                ref_fft = spfft.fft(ref_irf_dum)[1:nirf//2+1]
                ci_fft = spfft.fft(ci_irf_dum)[1:nirf//2+1]
                cross_spec = np.multiply(np.conj(ci_fft),ref_fft)
                wt_cross_spec = wt_cross_spec + np.multiply(lor_psd[0,:],cross_spec)
                pow_spec_ref = np.square(np.abs(ref_fft))
                pow_spec_ci = np.square(np.abs(ci_fft))
                wt_pow_spec_ref = wt_pow_spec_ref + np.multiply(lor_psd[0,:],pow_spec_ref)
                wt_pow_spec_ci = wt_pow_spec_ci + np.multiply(lor_psd[0,:],pow_spec_ci)
                mod_sig_psd = mod_sig_psd + lor_psd
    return wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd



######################################################################################################################
############### The functions below are used to generate lags and power spectra using numerical ######################
############### Monte Carlo simulation of light curves ###############################################################
######################################################################################################################


def run_simulation(rad,iseed_start,nsims,nbins,expquery,multquery,dt,deltau_scale,nirf,segnbin,bfac,lfreq,q,rms,gamma_par,\
                   e_seed,energy,encomb,i_rsigmax,irf_nbins,irf_binedgefrac,ldisk_disp,lseed_disp,lheat,rcor,\
                   disktocor_frac,cortodisk_frac,thermal_frac):
    # Runs simulation of propagating fluctuations following the approach of Arevalo & Uttley (2006), assuming no 
    # viscous diffusion effects.  The function runs nsims simulations and makes binned frequency-dependent spectral-
    # timing products as output.  Besides the simulation information, the function takes the disk, seed and heating
    # luminosity impulse responses as input, and convolves them with the simulated mass accretion rate time-series to 
    # produce light curves of these quantities and gamma, which are used to determine the flux light curves.
    
    # First make linear-rebinned IRFs for the disk, seed and coronal heating luminosities.    
    disk_disp_irf, disk_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,ldisk_disp,deltau_scale,nirf)
    seed_disp_irf, seed_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,lseed_disp,deltau_scale,nirf)
    heat_irf, heat_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,lheat,deltau_scale,nirf)
    
    # imin and imax set the limits of the simulated light curves to be used for analysis, to avoid the region of 
    # overlap in the convolution.
    imin = 2*int(np.sum(irf_nbins[:i_rsigmax+1]))
    imax = int(nbins - imin)
    # Set up the spectral timing quantities to be summed over simulations
    sum_pow_spec_ref = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_ci = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_refsq = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_cisq = np.zeros((len(encomb),segnbin//2))
    sum_cross_spec = np.zeros((len(encomb),segnbin//2),dtype=complex)
    nsegtot = 0 # Total number of light curve segments used to make the spectral-timing products
    print("Minimum and maximum bins: ",imin, imax)
    for i in range(nsims):
        iseed = iseed_start + i # Use a different but reproducible seed for each simulated light curve
        
        # We simulate the mdot, disk, seed, heating and photon index light curves using sim_lcs
        lc_disk_disp, lc_disk_rev, lc_seed_disp, lc_seed_rev, lc_heat,\
        lc_gamma, mdot = sim_lcs(rad,iseed,dt,nbins,expquery,multquery,lfreq,q,rms,gamma_par,i_rsigmax,irf_nbins,\
            deltau_scale,rcor,disktocor_frac,cortodisk_frac,thermal_frac,disk_disp_irf,disk_outer,seed_disp_irf,\
                                 seed_outer,heat_irf,heat_outer)
        
        print("Simulated ",i+1," light curves")
        print("rms values, disk disp: ",np.std(lc_disk_disp)," seed disp: ",np.std(lc_seed_disp)," heat: ",np.std(lc_heat))
        print("disk due to heating: ",np.std(lc_disk_rev)," seed due to heating: ",np.std(lc_seed_rev))
        print("Total mean luminosity (absorption corrected):",\
            np.mean(lc_disk_disp+lc_disk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lc_seed_disp+lc_heat+lc_seed_rev/thermal_frac)))

        N_0 = (1-np.sum(cortodisk_frac))*\
            np.multiply((lc_seed_disp[imin:imax] + lc_seed_rev[imin:imax]),lc_gamma[imin:imax]-1)/e_seed**2
        for ii in range(len(encomb)): # note that encomb values j,k = 0 correspond to the disk, values j, k > 0 correspond to
      #  an index (j,k) -1 in the 'energy' array of monochromatic energies for PL fluxes
            j = encomb[ii,0]
            k = encomb[ii,1]
            if (j == 0):  # define channel of interest lcs
                lc_ci = lc_disk_disp[imin:imax] + lc_disk_rev[imin:imax]
            else:
                lc_ci = np.multiply(N_0,np.power(energy[j-1]/e_seed,-1.*lc_gamma[imin:imax]))

            if (k == 0):  # define channel of interest lcs
                lc_ref = lc_disk_disp[imin:imax] + lc_disk_rev[imin:imax]
            else:
                lc_ref = np.multiply(N_0,np.power(energy[k-1]/e_seed,-1.*lc_gamma[imin:imax]))

            # Now measure the cross-spectrum and psds for the two bands
            pow_spec_ref, pow_spec_ci, pow_spec_refsq, pow_spec_cisq, cross_spec, nseg = \
                                                                make_cs_psd(segnbin,dt,lc_ref,lc_ci)
                                                             
            sum_pow_spec_ref[ii,:] = sum_pow_spec_ref[ii,:] + pow_spec_ref/nsims # Record running averages
            sum_pow_spec_refsq[ii,:] = sum_pow_spec_refsq[ii,:] + pow_spec_refsq/nsims
            sum_pow_spec_ci[ii,:] = sum_pow_spec_ci[ii,:] + pow_spec_ci/nsims # Channel-of-interest power spectrum
            sum_pow_spec_cisq[ii,:] = sum_pow_spec_cisq[ii,:] + pow_spec_cisq/nsims        
            sum_cross_spec[ii,:] = sum_cross_spec[ii,:] + cross_spec/nsims # And the cross-spectrum
            if (ii == 0):
                nsegtot = nsegtot + nseg # Keep track of total number of segments (just do this once per lc simulation)

    for i in range(len(encomb)):
        # Now bin up the average PSD and cross-spectrum in frequency (from the latter get the lags)
        bin_freq, bin_psd_ref, bin_psd_referr, bin_psd_ci, bin_psd_cierr, bin_phlag, bin_tlag = \
                bin_spectime_psd(nsegtot,segnbin,bfac,dt,sum_pow_spec_ref[i,:],sum_pow_spec_ci[i,:],\
                                 sum_pow_spec_refsq[i,:],sum_pow_spec_cisq[i,:],sum_cross_spec[i,:])
        if (i == 0): # Set up arrays to contain each of the spectral-timing products for all energies
            bin_psd_ciarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_cierrarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_refarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_referrarr = np.zeros((len(encomb),len(bin_freq)))            
            bin_phlagarr = np.zeros((len(encomb),len(bin_freq)))
            bin_tlagarr = np.zeros((len(encomb),len(bin_freq)))
        bin_psd_ciarr[i,:] = bin_psd_ci
        bin_psd_cierrarr[i,:] = bin_psd_cierr
        bin_psd_refarr[i,:] = bin_psd_ref
        bin_psd_referrarr[i,:] = bin_psd_referr
        bin_phlagarr[i,:] = bin_phlag
        bin_tlagarr[i,:] = bin_tlag            
            
    return bin_freq, bin_psd_refarr, bin_psd_referrarr, bin_psd_ciarr, bin_psd_cierrarr, bin_phlagarr, bin_tlagarr
    

def run_simulation_lcout(rad,iseed_start,nsims,nbins,expquery,multquery,dt,deltau_scale,nirf,segnbin,bfac,lfreq,q,rms,gamma_par,\
                   e_seed,energy,encomb,i_rsigmax,irf_nbins,irf_binedgefrac,ldisk_disp,lseed_disp,lheat,rcor,\
                   disktocor_frac,cortodisk_frac,thermal_frac):
    # Runs simulation of propagating fluctuations following the approach of Arevalo & Uttley (2006), assuming no
    # viscous diffusion effects.  The function runs nsims simulations and makes binned frequency-dependent spectral-
    # timing products as output.  Besides the simulation information, the function takes the disk, seed and heating
    # luminosity impulse responses as input, and convolves them with the simulated mass accretion rate time-series to
    # produce light curves of these quantities and gamma, which are used to determine the flux light curves.
    
    lc_out = np.zeros((len(energy)+1,nbins))
    
    # First make linear-rebinned IRFs for the disk, seed and coronal heating luminosities.
    disk_disp_irf, disk_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,ldisk_disp,deltau_scale,nirf)
    seed_disp_irf, seed_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,lseed_disp,deltau_scale,nirf)
    heat_irf, heat_outer = linear_rebin_irf(dt,i_rsigmax,irf_nbins,irf_binedgefrac,lheat,deltau_scale,nirf)
    
    # imin and imax set the limits of the simulated light curves to be used for analysis, to avoid the region of
    # overlap in the convolution.
    imin = 2*int(np.sum(irf_nbins[:i_rsigmax+1]))
    imax = int(nbins - imin)
    # Set up the spectral timing quantities to be summed over simulations
    sum_pow_spec_ref = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_ci = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_refsq = np.zeros((len(encomb),segnbin//2))
    sum_pow_spec_cisq = np.zeros((len(encomb),segnbin//2))
    sum_cross_spec = np.zeros((len(encomb),segnbin//2),dtype=complex)
    nsegtot = 0 # Total number of light curve segments used to make the spectral-timing products
    print("Minimum and maximum bins: ",imin, imax)
    for i in range(nsims):
        iseed = iseed_start + i # Use a different but reproducible seed for each simulated light curve
        
        # We simulate the mdot, disk, seed, heating and photon index light curves using sim_lcs
        lc_disk_disp, lc_disk_rev, lc_seed_disp, lc_seed_rev, lc_heat,\
        lc_gamma, mdot = sim_lcs(rad,iseed,dt,nbins,expquery,multquery,lfreq,q,rms,gamma_par,i_rsigmax,irf_nbins,\
            deltau_scale,rcor,disktocor_frac,cortodisk_frac,thermal_frac,disk_disp_irf,disk_outer,seed_disp_irf,\
                                 seed_outer,heat_irf,heat_outer)
        
        print("Simulated ",i+1," light curves")
        print("rms values, disk disp: ",np.std(lc_disk_disp)," seed disp: ",np.std(lc_seed_disp)," heat: ",np.std(lc_heat))
        print("disk due to heating: ",np.std(lc_disk_rev)," seed due to heating: ",np.std(lc_seed_rev))
        print("Total mean luminosity (absorption corrected):",\
            np.mean(lc_disk_disp+lc_disk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lc_seed_disp+lc_heat+lc_seed_rev/thermal_frac)))

        N_0 = (1-np.sum(cortodisk_frac))*\
            np.multiply((lc_seed_disp[imin:imax] + lc_seed_rev[imin:imax]),lc_gamma[imin:imax]-1)/e_seed**2
        for ii in range(len(encomb)): # note that encomb values j,k = 0 correspond to the disk, values j, k > 0 correspond to
      #  an index (j,k) -1 in the 'energy' array of monochromatic energies for PL fluxes
            j = encomb[ii,0]
            k = encomb[ii,1]
            if (j == 0):  # define channel of interest lcs
                lc_ci = lc_disk_disp[imin:imax] + lc_disk_rev[imin:imax]
            else:
                lc_ci = np.multiply(N_0,np.power(energy[j-1]/e_seed,-1.*lc_gamma[imin:imax]))

            if (k == 0):  # define channel of interest lcs
                lc_ref = lc_disk_disp[imin:imax] + lc_disk_rev[imin:imax]
            else:
                lc_ref = np.multiply(N_0,np.power(energy[k-1]/e_seed,-1.*lc_gamma[imin:imax]))
            if i == 0:
                lc_out[j,:len(lc_ci)] = lc_ci
                lc_out[k,:len(lc_ref)] = lc_ref

            # Now measure the cross-spectrum and psds for the two bands
            pow_spec_ref, pow_spec_ci, pow_spec_refsq, pow_spec_cisq, cross_spec, nseg = \
                                                                make_cs_psd(segnbin,dt,lc_ref,lc_ci)
                                                             
            sum_pow_spec_ref[ii,:] = sum_pow_spec_ref[ii,:] + pow_spec_ref/nsims # Record running averages
            sum_pow_spec_refsq[ii,:] = sum_pow_spec_refsq[ii,:] + pow_spec_refsq/nsims
            sum_pow_spec_ci[ii,:] = sum_pow_spec_ci[ii,:] + pow_spec_ci/nsims # Channel-of-interest power spectrum
            sum_pow_spec_cisq[ii,:] = sum_pow_spec_cisq[ii,:] + pow_spec_cisq/nsims
            sum_cross_spec[ii,:] = sum_cross_spec[ii,:] + cross_spec/nsims # And the cross-spectrum
            if (ii == 0):
                nsegtot = nsegtot + nseg # Keep track of total number of segments (just do this once per lc simulation)

    for i in range(len(encomb)):
        # Now bin up the average PSD and cross-spectrum in frequency (from the latter get the lags)
        bin_freq, bin_psd_ref, bin_psd_referr, bin_psd_ci, bin_psd_cierr, bin_phlag, bin_tlag = \
                bin_spectime_psd(nsegtot,segnbin,bfac,dt,sum_pow_spec_ref[i,:],sum_pow_spec_ci[i,:],\
                                 sum_pow_spec_refsq[i,:],sum_pow_spec_cisq[i,:],sum_cross_spec[i,:])
        if (i == 0): # Set up arrays to contain each of the spectral-timing products for all energies
            bin_psd_ciarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_cierrarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_refarr = np.zeros((len(encomb),len(bin_freq)))
            bin_psd_referrarr = np.zeros((len(encomb),len(bin_freq)))
            bin_phlagarr = np.zeros((len(encomb),len(bin_freq)))
            bin_tlagarr = np.zeros((len(encomb),len(bin_freq)))
        bin_psd_ciarr[i,:] = bin_psd_ci
        bin_psd_cierrarr[i,:] = bin_psd_cierr
        bin_psd_refarr[i,:] = bin_psd_ref
        bin_psd_referrarr[i,:] = bin_psd_referr
        bin_phlagarr[i,:] = bin_phlag
        bin_tlagarr[i,:] = bin_tlag
            
    return bin_freq, bin_psd_refarr, bin_psd_referrarr, bin_psd_ciarr, bin_psd_cierrarr, bin_phlagarr, bin_tlagarr, lc_out


def sim_lcs(rad,iseed,dt,nbins,expquery,multquery,lfreq,q,rms,gamma_par,i_rsigmax,irf_nbins,deltau_scale,rcor,\
            disktocor_frac,cortodisk_frac,thermal_frac,disk_disp_irf,disk_outer,seed_disp_irf,seed_outer,
            heat_irf,heat_outer):

    np.random.seed(iseed)
    mdot = np.zeros(nbins)+1.0
    deltam = np.zeros(nbins)
    lc_disk_disp = np.zeros(nbins)+disk_outer
    lc_seed_disp = np.zeros(nbins)+seed_outer
    lc_heat = np.zeros(nbins)+heat_outer
    for i in range(len(deltau_scale)-1,-1,-1):
        if (i == i_rsigmax or (i < i_rsigmax and deltau_scale[i+1] > 0)):
            print("Simulating radius ",i,", r = ",rad[i],", deltau = ",deltau_scale[i])
            if (rms[i] > 0.):
                lor_params = [lfreq[i],q[i],rms[i]]
                if (expquery == 'y'):
                    exp_params = np.copy(lor_params)
                    rmssq = lor_params[2]**2.
                    linrmssq = np.log(rmssq+1.0)
                    exp_params[2] = exp_params[2]*np.sqrt(linrmssq/rmssq)
                    deltam = tksim(dt, nbins, lorentz_q, exp_params)
                    deltam = np.exp(deltam)
                else:
                    lor_params = [lfreq[i],q[i],rms[i]]
                    deltam = 1. + tksim(dt, nbins, lorentz_q, lor_params)
                deltam = deltam/np.mean(deltam)
                print("Fractional rms of simulated light curve = ", np.std(deltam))
            else:
                np.ndarray.fill(deltam, 1.0)
            print("Mdot lengths",len(mdot),len(deltam))
            if (multquery == 'y'):
                mdot = np.multiply(mdot,deltam)
            else:
                mdot = mdot + (deltam-1.0)
        
            if (deltau_scale[i] > 0):
                irfbin_start = int(np.sum(irf_nbins[i:i_rsigmax+1]) - irf_nbins[i])
                irfbin_stop = int((irfbin_start + irf_nbins[i])-1)
            else:
                irfbin_start = irfbin_stop+1
                irfbin_stop = irfbin_start
            nirfbins = (irfbin_stop - irfbin_start)+1
            conv_disk = np.zeros(nbins)
            conv_seed = np.zeros(nbins)
            conv_heat = np.zeros(nbins)
            conv_disk[:nirfbins] = disk_disp_irf[irfbin_start:irfbin_stop+1]*dt
            conv_seed[:nirfbins] = seed_disp_irf[irfbin_start:irfbin_stop+1]*dt
            conv_heat[:nirfbins] = heat_irf[irfbin_start:irfbin_stop+1]*dt
            lc_disk_disp = lc_disk_disp + np.real(spfft.ifft(np.multiply(spfft.fft(mdot),spfft.fft(conv_disk))))
            lc_seed_disp = lc_seed_disp + np.real(spfft.ifft(np.multiply(spfft.fft(mdot),spfft.fft(conv_seed))))
            lc_heat = lc_heat + np.real(spfft.ifft(np.multiply(spfft.fft(mdot),spfft.fft(conv_heat))))
            mdot = np.roll(mdot,int(nirfbins))

    f_rev = np.sum(cortodisk_frac) # The fraction of coronal luminosity that is reprocessed by the disk
    f_return = np.sum(cortodisk_frac*disktocor_frac) # the fraction of coronal luminosity that is reprocessed by the disk
    # and returned to the corona
    lc_seed_rev = thermal_frac*f_return*(lc_seed_disp+lc_heat)/(1-f_return)
    lc_disk_rev = thermal_frac*(f_rev-f_return)*(lc_seed_disp+lc_heat)/(1-f_return)

    lc_gamma = gamma_par[0]*np.power((np.divide((lc_seed_disp + lc_seed_rev),lc_heat)),gamma_par[1])
        
    return lc_disk_disp, lc_disk_rev, lc_seed_disp, lc_seed_rev, lc_heat, lc_gamma, mdot
        

def tksim(dt, ntimes, psmodel, pspar):
# Implementation of Timmer & Koenig method to simulate a noise process for an arbitrary power-spectral shape 
# (Timmer & Koenig, 1995, A&A, 300, 707), based on an original Python implementation by Dimitrios Emmanoulopoulos
    nfreq = ntimes//2
    ft_full = np.full(ntimes,complex(0.0,0.0))
    f = np.arange(1.0, nfreq+1, 1.0)/(ntimes*dt)
    modpow = np.multiply(psmodel(f, *pspar),nfreq/dt)
    gdev1=np.random.normal(size=(1,nfreq))
    gdev2=np.random.normal(size=(1,nfreq))
    ft_re=np.multiply(np.sqrt(np.multiply(modpow,0.5)),np.reshape(gdev1,nfreq))  
    ft_im=np.multiply(np.sqrt(np.multiply(modpow,0.5)),np.reshape(gdev2,nfreq))
    ft_pos = ft_re + ft_im*1j
    ft_neg = np.conj(ft_pos) # sets negative frequencies as complex conjugate of positive (for real-valued LC)
    # From documentation for numpy.fft: a[0] should contain the zero frequency term,
    # a[1:n//2] should contain the positive-frequency terms,
    # a[n//2 + 1:] should contain the negative-frequency terms, in increasing order starting from the most
    # negative frequency. For an even number of input points, A[n//2] represents the sum of the values at the
    # positive and negative Nyquist frequencies, as the two are aliased together.
    ft_full[1:nfreq+1] = ft_full[1:nfreq+1]+ft_pos[0,:]
    ft_full[nfreq+1:] = ft_full[nfreq+1:]+ft_neg[0,-2::-1]
    ft_full[nfreq]=complex(ft_full.real[nfreq],0.0) # For symmetry need to make Nyquist freq value real - note that
    # neglecting to do this causes a small imaginary component to appear in the inverse FFT
    ift_full=np.fft.ifft(ft_full)                            
    lc=ift_full.real
    return lc


def make_cs_psd(segnbin,dt,lc_ref,lc_ci):
    
    nseg = int(len(lc_ref)/segnbin)
    pow_spec_ref = np.zeros(segnbin//2)
    pow_spec_ci = np.zeros(segnbin//2)
    pow_spec_refsq = np.zeros(segnbin//2)
    pow_spec_cisq = np.zeros(segnbin//2)
    cross_spec = np.zeros(segnbin//2,dtype=complex)
    A_rms = 2*dt/segnbin
    mean_power = np.zeros(nseg)
    mean_ref = np.zeros(nseg)
    for i in range(nseg):        
        jmin = i*segnbin
        jmax = (i*segnbin)+segnbin
        refsqmean = np.square(np.mean(lc_ref[jmin:jmax]))
        cisqmean = np.square(np.mean(lc_ci[jmin:jmax]))
        ref_fft = spfft.fft(lc_ref[jmin:jmax])[1:(segnbin//2)+1]
        ci_fft = spfft.fft(lc_ci[jmin:jmax])[1:(segnbin//2)+1]
        cross_spec = cross_spec + A_rms*np.multiply(np.conj(ci_fft),ref_fft)/np.sqrt(refsqmean*cisqmean)
        pow_spec_ref = pow_spec_ref + A_rms*np.square(np.abs(ref_fft))/refsqmean
        pow_spec_ci = pow_spec_ci + A_rms*np.square(np.abs(ci_fft))/cisqmean
        pow_spec_refsq = pow_spec_refsq + np.square(A_rms*np.square(np.abs(ref_fft))/refsqmean)
        pow_spec_cisq = pow_spec_cisq + np.square(A_rms*np.square(np.abs(ci_fft))/cisqmean)
    pow_spec_ref = pow_spec_ref/nseg
    pow_spec_ci = pow_spec_ci/nseg
    pow_spec_refsq = pow_spec_refsq/nseg
    pow_spec_cisq = pow_spec_cisq/nseg    
    cross_spec = cross_spec/nseg
    return pow_spec_ref, pow_spec_ci, pow_spec_refsq, pow_spec_cisq, cross_spec, nseg
        

def bin_spectime_psd(nsegtot,segnbin,bfac,dt,sum_psd_ref,sum_psd_ci,sum_psd_refsq,sum_psd_cisq,sum_cross_spec):
    
    minfreq = 1./(dt*segnbin) # the minimum Fourier frequency
    freq = np.linspace(minfreq, len(sum_psd_ref)*minfreq, len(sum_psd_ref))

    fstart = freq[0]
    istart = 0
    j = -1
    nfbins = np.zeros(1)
    # Now assign individual Fourier frequencies to bins
    for i in range (0,len(freq)): # Now bin in frequency        
        if (freq[i] > fstart*bfac):
            j = j + 1
            if (j == 0):
                nfbins[0] = i - istart
            else:
                nfbins = np.append(nfbins,(i-istart))
            istart = i
            fstart = freq[i]

    # Now bin up
    bin_freq = np.zeros(j+1)
    bin_psd_ref = np.zeros(j+1)
    bin_psd_ci = np.zeros(j+1)
    bin_psd_referr = np.zeros(j+1)
    bin_psd_cierr = np.zeros(j+1)
    bin_cross_spec = np.zeros(j+1,dtype=complex)
    for i in range(j+1):
        imax = int(np.sum(nfbins[:i+1]))
        imin = int(imax - nfbins[i])
        bin_freq[i] = np.mean(freq[imin:imax])
        bin_psd_ref[i] = np.mean(sum_psd_ref[imin:imax])
        bin_psd_ci[i] = np.mean(sum_psd_ci[imin:imax])
        bin_psd_referr[i] = np.mean(sum_psd_refsq[imin:imax])-np.square(bin_psd_ref[i])
        bin_psd_cierr[i] = np.mean(sum_psd_cisq[imin:imax])-np.square(bin_psd_ci[i])
        bin_cross_spec[i] = np.mean(sum_cross_spec[imin:imax])
#        print("now: ",i,imin,imax,freq[imin],freq[imax-1],bin_freq[i],bin_cross_spec[i])

    ntotal = nsegtot*nfbins
    bin_psd_referr = np.sqrt(bin_psd_referr/(ntotal-1))
    bin_psd_cierr = np.sqrt(bin_psd_cierr/(ntotal-1))
    bin_phlag = np.angle(bin_cross_spec)
    bin_tlag = np.divide(bin_phlag,2*np.pi*bin_freq)
    
    return bin_freq, bin_psd_ref, bin_psd_referr, bin_psd_ci, bin_psd_cierr, bin_phlag, bin_tlag


######################################################################################################################
################################# Miscellaneous plotting and scripting functions #####################################
######################################################################################################################


def run_stprod_mono(rin,rout,nrad,rcor,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac):


    rsigmax = (freq_rsigmax*t_scale*disk_tau_par[0]*10**(-1.*disk_tau_par[1]))**(-1/(1.5+disk_tau_par[1]))
    print("Minimum frequency =",freq_rsigmax,"Hz and maximum signal radius = ",rsigmax,"Rg")
    rout = max([2*rsigmax,rout])
        
    #First set up radial grid
    rad_edge1 = np.logspace(np.log10(rin), np.log10(rsigmax), nrad+1)  # set up radial bin edges
    rad_edge2 = np.logspace(np.log10(rsigmax), np.log10(rout), nrad+1)
    rad_edge = np.append(rad_edge1,rad_edge2[1:])
    nrad = nrad*2
    
    #First set up radial grid
    rad = np.sqrt(rad_edge[1:]*rad_edge[:-1])  # Define radial bin centres as geometric mean of bin edges
    rad_area = np.pi*(np.square(rad_edge[1:])-np.square(rad_edge[:-1]))
    # Now reset coronal outer radius (i.e. disk inner radius) to be at the nearest radial bin edges
    rcor, i_rcor = find_nearest(rad_edge,rcor)
    disp_frac, seed_frac_flow, heat_frac_flow = \
        calc_dispfrac(rad,rad_area,rin,rcor,seedff_norm,seedff_ind,heatff_norm,heatff_ind)
    # Normalised viscous dissipation and parameters for seed and heating fractions within coronal flow.
    # Reset maximum signal radius to nearest radial bin value
    rsigmax, i_rsigmax = find_nearest(rad,rsigmax)
    print("Coronal radius reset to nearest radial bin edge: ",rcor)
    print("Maximum signal radius reset to: ",rsigmax)
    # Calculate illumination of corona by disk and vice-versa
    omega_cor, disktocor_frac, cortodisk_frac = calc_illumination_fracs(rad,rad_area,cor_geometry,geopar)
    print(np.sum(disp_frac[i_rcor:]*disktocor_frac[i_rcor:])/np.sum(disp_frac[i_rcor:])," of the disk flux is intercepted by the corona")
    print(np.sum(cortodisk_frac[i_rcor:])," of the coronal flux is intercepted by the disk")
    # Now calculate radial dependence of timing parameters and calculate propagation delays
    tau, lfreq, q, rms = calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, disk_tau_par, cor_tau_par,\
                                        lor_model, lor_par)
    deltau = calc_propagation_params(rad,rad_edge,rcor,disk_prop_par,cor_prop_par)
    # Calculate the radially-dependent direct disk, seed photon and coronal heating luminosities
    ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev =\
        calc_radial_time_response(rad,rcor,disp_frac,disktocor_frac,cortodisk_frac,thermal_frac,seed_frac_flow,heat_frac_flow)
    print("Seed energy :",e_seed)
    print("Dissipation disk luminosity: ",np.sum(ldisk_disp))
    print("Dissipation seed luminosity: ",np.sum(lseed_disp))
    print("Coronal heating luminosity: ",np.sum(lheat))
    print("Disk luminosity due to heating by corona: ",np.sum(ldisk_rev))
    print("Seed luminosity due to heating by corona: ",np.sum(lseed_rev))
    print("Total luminosity: ",np.sum(ldisk_disp+ldisk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lseed_disp+lheat+lseed_rev/thermal_frac)))
    
    # Calculate IRFs for the mono-energetic case
    gamma_mean, gamma_irf, flux_irf_mono, disk_irf_mono, seed_irf_mono =calc_irfs_mono(gamma_par,e_seed,ens_mono,\
                                ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev)
    print("Mean gamma is:", gamma_mean)
    freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf,\
        ci_irf, ci_mean, ci_outer =\
    calculate_stprod_mono(nirf_mult,ens_mono,encomb,flux_irf_mono,disk_irf_mono,gamma_irf,
                              deltau,min_deltau_frac,i_rsigmax,lfreq,q,rms,t_scale)
    
    return gamma_mean, freq, phlag, tlag, psd_ci, psd_ref, tau, i_rcor, i_rsigmax, deltau, ldisk_disp,\
    lseed_disp, lheat, ldisk_rev, lseed_rev, flux_irf_mono, mod_sig_psd


def run_stprod_disken(rin,rout,nrad,rcor,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac,band_mult):


    rsigmax = (freq_rsigmax*t_scale*disk_tau_par[0]*10**(-1.*disk_tau_par[1]))**(-1/(1.5+disk_tau_par[1]))
    print("Minimum frequency =",freq_rsigmax,"Hz and maximum signal radius = ",rsigmax,"Rg")
    rout = max([2*rsigmax,rout])
        
    #First set up radial grid
    rad_edge1 = np.logspace(np.log10(rin), np.log10(rsigmax), nrad+1)  # set up radial bin edges
    rad_edge2 = np.logspace(np.log10(rsigmax), np.log10(rout), nrad+1)
    rad_edge = np.append(rad_edge1,rad_edge2[1:])
    nrad = nrad*2
    
    #First set up radial grid
    rad = np.sqrt(rad_edge[1:]*rad_edge[:-1])  # Define radial bin centres as geometric mean of bin edges
    rad_area = np.pi*(np.square(rad_edge[1:])-np.square(rad_edge[:-1]))
    # Now reset coronal outer radius (i.e. disk inner radius) to be at the nearest radial bin edges
    rcor, i_rcor = find_nearest(rad_edge,rcor)
    disp_frac, seed_frac_flow, heat_frac_flow = \
        calc_dispfrac(rad,rad_area,rin,rcor,seedff_norm,seedff_ind,heatff_norm,heatff_ind)
    # Normalised viscous dissipation and parameters for seed and heating fractions within coronal flow.
    # Reset maximum signal radius to nearest radial bin value
    rsigmax, i_rsigmax = find_nearest(rad,rsigmax)
    print("Coronal radius reset to nearest radial bin edge: ",rcor)
    print("Maximum signal radius reset to: ",rsigmax)
    # Calculate illumination of corona by disk and vice-versa
    omega_cor, disktocor_frac, cortodisk_frac = calc_illumination_fracs(rad,rad_area,cor_geometry,geopar)
    print(np.sum(disp_frac[i_rcor:]*disktocor_frac[i_rcor:])/np.sum(disp_frac[i_rcor:])," of the disk flux is intercepted by the corona")
    print(np.sum(cortodisk_frac[i_rcor:])," of the coronal flux is intercepted by the disk")
    # Now calculate radial dependence of timing parameters and calculate propagation delays
    tau, lfreq, q, rms = calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, disk_tau_par, cor_tau_par,\
                                        lor_model, lor_par)
    deltau = calc_propagation_params(rad,rad_edge,rcor,disk_prop_par,cor_prop_par)
    # Calculate the radially-dependent direct disk, seed photon and coronal heating luminosities
    ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev =\
        calc_radial_time_response(rad,rcor,disp_frac,disktocor_frac,cortodisk_frac,thermal_frac,seed_frac_flow,heat_frac_flow)
    print("Seed energy :",e_seed)
    print("Dissipation disk luminosity: ",np.sum(ldisk_disp))
    print("Dissipation seed luminosity: ",np.sum(lseed_disp))
    print("Coronal heating luminosity: ",np.sum(lheat))
    print("Disk luminosity due to heating by corona: ",np.sum(ldisk_rev))
    print("Seed luminosity due to heating by corona: ",np.sum(lseed_rev))
    print("Total luminosity: ",np.sum(ldisk_disp+ldisk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lseed_disp+lheat+lseed_rev/thermal_frac)))

    ens_mono_pl = ens_mono[ens_mono > 0]
    ens_mono_disk = np.abs(ens_mono[ens_mono < 0])
    # Calculate disk luminosity and PL flux IRFs
    gamma_mean, gamma_irf, flux_irf_mono_pl, disklum_irf_mono, seed_irf_mono =calc_irfs_mono(gamma_par,e_seed,ens_mono_pl,\
                                ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev)
    print("Mean gamma is:", gamma_mean)

    flux_irf_mono = np.zeros((len(ens_mono),len(ldisk_disp)))
    flux_disk_const = np.zeros((len(ens_mono_disk),len(ldisk_disp)))
    flux_irf_mono[len(ens_mono_disk):,:] = flux_irf_mono_pl
    for i, disk_en in enumerate(ens_mono_disk):
        disk_eband = [disk_en,disk_en*band_mult[i]]
        band_frac, kT_rad, disk_irf_mono, lseed_band, ldiskband_disp, ldiskband_rev = calc_disk_band(disp_frac,cortodisk_frac,
                                    disktocor_frac,lseed_disp,lheat,thermal_frac,rcor,rad,rad_area,disk_eband,kTvar=True)
        band_frac2, kT_rad2, disk_irf_mono2, lseed_band2, ldiskband_disp2, ldiskband_rev2 = calc_disk_band(disp_frac,cortodisk_frac,
                                    disktocor_frac,lseed_disp,lheat,thermal_frac,rcor,rad,rad_area,disk_eband,kTvar=False)
        flux_irf_mono[i,:] = disk_irf_mono
        flux_disk_const[i,:] = disk_irf_mono2
    freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf,\
        ci_irf, ci_mean, ci_outer =\
    calculate_stprod_mono(nirf_mult,ens_mono,encomb,flux_irf_mono,disklum_irf_mono,gamma_irf,
                              deltau,min_deltau_frac,i_rsigmax,lfreq,q,rms,t_scale)
    for i, i_ens in enumerate(encomb):
        if (ens_mono[i_ens[0]-1] < 0):
            psd_ci[i,:] = psd_ci[i,:]/(np.sum(flux_disk_const[i_ens[0]-1]))**2
        if (ens_mono[i_ens[1]-1] < 0):
            psd_ref[i,:] = psd_ref[i,:]/(np.sum(flux_disk_const[i_ens[1]-1]))**2
            
    
    return gamma_mean, freq, phlag, tlag, psd_ci, psd_ref, tau, i_rcor, i_rsigmax, deltau, ldisk_disp,\
    lseed_disp, lheat, ldisk_rev, lseed_rev, flux_irf_mono, mod_sig_psd


def plot_lags_psds(lags_list, psds_list, vlines_list, axis_names, freqlim, tlaglim, tlaglim_in, psdlim,
    inset_yticks, leg_title=None, figfile=None):
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8),sharex=True,gridspec_kw={'height_ratios':[5,3]})
    fig.subplots_adjust(hspace=0)
    ax3 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax1, [0.65,0.62,0.32,0.37])
    ax3.set_axes_locator(ip)
    for x in lags_list:
        freq = x[0]
        tlag = x[1]
        colour_val = x[2]
        ls_val = x[3]
        label_val = x[4]
        ax1.plot(freq,tlag,color=colour_val,linewidth=3,linestyle=ls_val,label=label_val)
        ax3.plot(freq,tlag,color=colour_val,linewidth=3,linestyle=ls_val)
    for x in psds_list:
        freq = x[0]
        psd = x[1]
        colour_val = x[2]
        ls_val = x[3]
        label_val = x[4]
        ax2.plot(freq,psd,color=colour_val,linewidth=3,linestyle=ls_val,label=label_val)
    if (vlines_list != None):
        for line in vlines_list:
            for ax in (ax1, ax2, ax3):
                ax.axvline(line[0],color=line[1],alpha=0.5,linewidth=3,linestyle=line[2])
    ax3.set_xlim(freqlim)
    ax3.set_ylim(tlaglim_in)
    ax3.set_xscale('log')
    ax3.set_yscale('linear')
    ax3.axhline(color='grey',linestyle='dotted')
    ax3.set_xticks([0.01,0.1,1,10])
    ax3.set_yticks(inset_yticks)
    ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    ax2.legend(fontsize=12,title=leg_title,title_fontsize=12)
    ax2.set_xlabel(axis_names[0],fontsize=14)
    ax1.set_ylabel(axis_names[1],fontsize=14)
    ax2.set_ylabel(axis_names[2],fontsize=14)
    ax3.set_xlabel(axis_names[0],fontsize=10)
    ax3.set_ylabel(axis_names[1],fontsize=10)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_xlim(freqlim)
    ax1.set_ylim(tlaglim)
    ax2.set_ylim(psdlim)
    for ax in (ax1, ax2, ax3):
        ax.tick_params(axis='x',labelsize=12, which='both', direction='in', top=True)
        ax.tick_params(axis='y',labelsize=12, which='both', direction='in', right=True)
        ax.tick_params(which='major', length=10)
        ax.tick_params(which='minor', length=6)
    if (figfile != None):
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()

    return


def geom_observables_calc(rin,rout,nrad,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,geom_list,
                    thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac):

    gamma_vals = np.zeros(len(geom_list))
    f_s_rev = np.zeros(len(geom_list))
    f_d_rev = np.zeros(len(geom_list))
    obsmax_disklag = np.zeros(len(geom_list))
    obsmin_disklag = np.zeros(len(geom_list))
    obsmax_pllag = np.zeros(len(geom_list))
    obsfreq_cross = np.zeros(len(geom_list))
    psd_ratio_pllo = np.zeros(len(geom_list))
    psd_ratio_plhi = np.zeros(len(geom_list))
    rcor = np.zeros(len(geom_list))
    hcor = np.zeros(len(geom_list))
    cone_angle = np.zeros(len(geom_list))
    freqmin_disklag = np.zeros(len(geom_list))
    cor_freq = np.zeros(len(geom_list))
    nyq_freq = np.zeros(len(geom_list))
    spectime_out = np.zeros((len(geom_list),500,8))
    
    j = 0
    for i, x in enumerate(geom_list):
        cor_geometry = x[0]
        rcor[i] = x[1]
        hcor[i] = x[2]
        cone_angle[i] = x[3]
        
        if (hcor[i] > 0.):
            rtop = rcor[i] + hcor[i]*np.tan(np.pi*(cone_angle[i]/180))
            geopar = [rcor[i],hcor[i],rtop,100,400]
        else:
            geopar = [rcor[i]]
        j = j + 1
        print("*******************************************************************")
        print("Calculating ",j," out of ",len(geom_list)," for geometry ",cor_geometry," with geometrical parameters: ",geopar)
        print("*******************************************************************")
        gamma_mean, freq, phlag, tlag, psd_ci, psd_ref, tau, i_rcor, i_rsigmax, deltau, ldisk_disp, lseed_disp,\
            lheat, ldisk_rev, lseed_rev, flux_irf_mono, mod_sig_psd = \
        run_stprod_mono(rin,rout,nrad,rcor[i],freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac)
        
        f_s_rev[i] = np.sum(lseed_rev)/np.sum(lseed_rev+lseed_disp)
        f_d_rev[i] = np.sum(ldisk_rev)/np.sum(ldisk_rev+ldisk_disp)
        gamma_vals[i] = gamma_mean
        obsmax_disklag[i] = np.amax(tlag[0,:])
        obsmax_pllag[i] = np.amax(tlag[1,:])
        obsmin_disklag[i] = np.amin(tlag[0,:])
        freqrcor, i_freqrcor = find_nearest(freq,1./(t_scale*tau[i_rcor]))
        imin = np.argmin(tlag[0,:])
        if ((np.amax(tlag[0,:imin]) > 0) & (np.amin(tlag[0,:imin]) < 0)):
            interplag = spinterp.interp1d(tlag[0,:imin],freq[:imin])
            obsfreq_cross[i] = interplag(0)
            print("Interpolating crossing frequency:",obsfreq_cross[i]," up to max soft lag frequency ",freq[imin],", coronal freq:",freqrcor)
        else:
            obsfreq_cross[i] = np.nan
            print("No crossing frequency interpolation possible:",obsfreq_cross[i]," up to max soft lag frequency ",freq[imin],", coronal freq:",freqrcor)
        psd_ratio_pllo[i] = psd_ref[1,i_freqrcor]/psd_ci[1,i_freqrcor]
        psd_ratio_plhi[i] = psd_ref[2,i_freqrcor]/psd_ci[2,i_freqrcor]
        nyq_freq[i] = freq[-1]
        cor_freq[i] = freqrcor
        freqmin_disklag[i] = freq[imin]
        print("Parameters: ",gamma_vals[i],obsmax_pllag[i],obsmin_disklag[i],obsfreq_cross[i])
        indices = np.unique(np.logspace(np.log2(1),np.log2(len(freq)),500,base=2,dtype=int))-1
        spectime_out[i,:len(indices),0] = freq[indices]
        spectime_out[i,:len(indices),1] = tlag[0,indices]
        spectime_out[i,:len(indices),2] = tlag[1,indices]
        spectime_out[i,:len(indices),3] = tlag[2,indices]
        spectime_out[i,:len(indices),4] = psd_ref[0,indices]
        spectime_out[i,:len(indices),5] = psd_ci[0,indices]
        spectime_out[i,:len(indices),6] = psd_ci[1,indices]
        spectime_out[i,:len(indices),7] = psd_ci[2,indices]
        print(rcor[i],hcor[i],cone_angle[i],f_s_rev[i],f_d_rev[i],gamma_vals[i],obsmax_disklag[i],obsmax_pllag[i],\
            obsmin_disklag[i],obsfreq_cross[i],psd_ratio_pllo[i],psd_ratio_plhi[i],freqmin_disklag[i],\
            cor_freq[i],nyq_freq[i])
        
    geom_stparams = np.column_stack((rcor,hcor,cone_angle,f_s_rev,f_d_rev,gamma_vals,obsmax_disklag,obsmax_pllag,\
        obsmin_disklag,obsfreq_cross,psd_ratio_pllo,psd_ratio_plhi,freqmin_disklag,cor_freq,nyq_freq))
        
    return geom_stparams, spectime_out
    
    
def geom_observables_calc_test(rin,rout,nrad,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,geom_list,
                    thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac):

    gamma_vals = np.zeros(len(geom_list))
    f_s_rev = np.zeros(len(geom_list))
    f_d_rev = np.zeros(len(geom_list))
    obsmax_disklag = np.zeros(len(geom_list))
    obsmin_disklag = np.zeros(len(geom_list))
    obsmax_pllag = np.zeros(len(geom_list))
    obsfreq_cross = np.zeros(len(geom_list))
    psd_ratio_pllo = np.zeros(len(geom_list))
    psd_ratio_plhi = np.zeros(len(geom_list))
    rcor = np.zeros(len(geom_list))
    hcor = np.zeros(len(geom_list))
    cone_angle = np.zeros(len(geom_list))
    
    j = 0
    for i, x in enumerate(geom_list):
        cor_geometry = x[0]
        rcor[i] = x[1]
        hcor[i] = x[2]
        cone_angle[i] = x[3]
        
        if (hcor[i] > 0.):
            rtop = rcor[i] + hcor[i]*np.tan(np.pi*(cone_angle[i]/180))
            geopar = [rcor[i],hcor[i],rtop,100,400]
        else:
            geopar = [rcor[i]]
        j = j + 1
        print("*******************************************************************")
        print("Calculating ",j," out of ",len(geom_list)," for geometry ",cor_geometry," with geometrical parameters: ",geopar)
        print("*******************************************************************")
        gamma_mean, freq, phlag, tlag, psd_ci, psd_ref, tau, i_rcor, i_rsigmax, deltau, ldisk_disp, lseed_disp,\
            lheat, ldisk_rev, lseed_rev, flux_irf_mono, mod_sig_psd = \
        run_stprod_mono(rin,rout,nrad,rcor[i],freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac)
        
        f_s_rev[i] = np.sum(lseed_rev)/np.sum(lseed_rev+lseed_disp)
        f_d_rev[i] = np.sum(ldisk_rev)/np.sum(ldisk_rev+ldisk_disp)
        gamma_vals[i] = gamma_mean
        obsmax_disklag[i] = np.amax(tlag[0,:])
        obsmax_pllag[i] = np.amax(tlag[1,:])
        obsmin_disklag[i] = np.amin(tlag[0,:])
        freqrcor, i_freqrcor = find_nearest(freq,1./(t_scale*tau[i_rcor]))
        if ((np.amax(tlag[0,:i_freqrcor]) > 0) & (np.amin(tlag[0,:i_freqrcor]) < 0)):
            interplag = spinterp.interp1d(tlag[0,:i_freqrcor],freq[:i_freqrcor])
            obsfreq_cross[i] = interplag(0)
            print("Interpolating:",obsfreq_cross[i])
        else:
            obsfreq_cross[i] = np.nan
            print("No interpolation possible:",obsfreq_cross[i])
        psd_ratio_pllo[i] = psd_ref[1,i_freqrcor]/psd_ci[1,i_freqrcor]
        psd_ratio_plhi[i] = psd_ref[2,i_freqrcor]/psd_ci[2,i_freqrcor]
        print("Parameters: ",gamma_vals[i],obsmax_pllag[i],obsmin_disklag[i],obsfreq_cross[i])
        if (hcor[i] == 30):
            testlag, testfreq = tlag[0,:i_freqrcor], freq[:i_freqrcor]
#        else:
#            testlag, testfreq = None, None
        
    geom_stparams = np.column_stack((rcor,hcor,cone_angle,f_s_rev,f_d_rev,gamma_vals,obsmax_disklag,obsmax_pllag,\
        obsmin_disklag,obsfreq_cross,psd_ratio_pllo,psd_ratio_plhi))
        
    return geom_stparams, testlag, testfreq
    
def plot_stparams_grid(input_lists,xcol,ycol,xlim_vals,ylim_vals,xlabel_val,ylabel_val,legend_on,legend_loc,fig_title,
                       figfile):
    plt.figure(figsize=(6,5))

    if (xlim_vals[1] == 'log'):
        xval_range = np.diff(np.log10(np.array(xlim_vals[0])))
    else:
        xval_range = np.diff(np.array(xlim_vals[0]))
    if (ylim_vals[1] == 'log'):
        yval_range = np.diff(np.log10(np.array(ylim_vals[0])))
    else:
        yval_range = np.diff(np.array(ylim_vals[0]))
        

    for inlist in input_lists:
        geom_pars = inlist[0]
        txt_list = inlist[1]
        colour_val = inlist[2]
        hcor_list = inlist[4]
        angle_list = inlist[5]
        if txt_list != None:
            txt_vals = txt_list[0]
            txt_off_scale = txt_list[1]
            txt_off = txt_list[2]
            txt_hjust = txt_list[3]
            txt_vjust = txt_list[4]
        # If spherical then plot single line
        if hcor_list == None:
            label_val = inlist[3]
            plt.plot(geom_pars[:,xcol],geom_pars[:,ycol],linewidth=3,color=colour_val,label=label_val)
            if txt_list != None:
                for i, radius in enumerate(txt_vals):
                    selected = geom_pars[geom_pars[:,0] == radius]
                    if (xlim_vals[1] == 'log'):
                        xtxt = selected[0,xcol]*10**(txt_off[i,0]*xval_range)
                    else:
                        xtxt = selected[0,xcol]+txt_off[i,0]*xval_range
                    if (ylim_vals[1] == 'log'):
                        ytxt = selected[0,ycol]*10**(txt_off[i,1]*yval_range)
                    else:
                        ytxt = selected[0,ycol]+txt_off[i,1]*yval_range

                    plt.text(xtxt,ytxt,str(int(round(radius)))+" $R_{g}$",color=colour_val,ha=txt_hjust[i],
                        va=txt_vjust[i])
        else:
            n = len(geom_pars)//2
            # First the constant hcor lines
            for i, hcor in enumerate(hcor_list):
                geom_pars2 = geom_pars[:n,:]
                selected = geom_pars2[geom_pars2[:,1] == hcor]
                if i == 0:
                    label_val = inlist[3]
                else:
                    label_val = None
                plt.plot(selected[:,xcol],selected[:,ycol],linewidth=3,color=colour_val,label=label_val)

                if txt_list != None:
                    if txt_vals[i] != 'off':
                        if txt_vals[i] == 'lower':
                            ival = 0
                            if (xlim_vals[1] == 'log'):
                                gridvec_x = np.log10(selected[0,xcol]/selected[3,xcol])/xval_range
                            else:
                                gridvec_x = (selected[0,xcol]-selected[3,xcol])/xval_range
                            if (ylim_vals[1] == 'log'):
                                gridvec_y = np.log10(selected[0,ycol]/selected[3,ycol])/yval_range
                            else:
                                gridvec_y = (selected[0,ycol]-selected[3,ycol])/yval_range
                        else:
                            ival = -1
                            if (xlim_vals[1] == 'log'):
                                gridvec_x = np.log10(selected[-1,xcol]/selected[-4,xcol])/xval_range
                            else:
                                gridvec_x = (selected[-1,xcol]-selected[-4,xcol])/xval_range
                            if (ylim_vals[1] == 'log'):
                                gridvec_y = np.log10(selected[-1,ycol]/selected[-4,ycol])/yval_range
                            else:
                                gridvec_y = (selected[-1,ycol]-selected[-4,ycol])/yval_range
                        
                        gridvec_xn = gridvec_x/np.sqrt(gridvec_x**2+gridvec_y**2)
                        gridvec_yn = gridvec_y/np.sqrt(gridvec_x**2+gridvec_y**2)

                        if (xlim_vals[1] == 'log'):
                            xtxt = selected[ival,xcol]*10**((txt_off_scale[i]*gridvec_xn+txt_off[i,0])*xval_range)
                        else:
                            xtxt = selected[ival,xcol]+(txt_off_scale[i]*gridvec_xn+txt_off[i,0])*xval_range
                        if (ylim_vals[1] == 'log'):
                            ytxt = selected[ival,ycol]*10**((txt_off_scale[i]*gridvec_yn+txt_off[i,1])*yval_range)
                        else:
                            ytxt = selected[ival,ycol]+(txt_off_scale[i]*gridvec_yn+txt_off[i,1])*yval_range
                        plt.text(xtxt,ytxt,str(round(hcor))+" $R_{g}$",color=colour_val,ha=txt_hjust[i],va=txt_vjust[i])
                            
            # Now the constant cone angle lines
            for i, angle in enumerate(angle_list):
                geom_pars2 = geom_pars[n:,:]
                selected = geom_pars2[geom_pars2[:,2] == angle]
                plt.plot(selected[:,xcol],selected[:,ycol],linewidth=3,color=colour_val,label=None)
                if txt_list != None:
                    j = i + len(hcor_list)
                    if txt_vals[j] != 'off':
                        if txt_vals[j] == 'lower':
                            ival = 0
                            if (xlim_vals[1] == 'log'):
                                gridvec_x = np.log10(selected[0,xcol]/selected[3,xcol])/xval_range
                            else:
                                gridvec_x = (selected[0,xcol]-selected[3,xcol])/xval_range
                            if (ylim_vals[1] == 'log'):
                                gridvec_y = np.log10(selected[0,ycol]/selected[3,ycol])/yval_range
                            else:
                                gridvec_y = (selected[0,ycol]-selected[3,ycol])/yval_range
                        else:
                            ival = -1
                            if (xlim_vals[1] == 'log'):
                                gridvec_x = np.log10(selected[-1,xcol]/selected[-4,xcol])/xval_range
                            else:
                                gridvec_x = (selected[-1,xcol]-selected[-4,xcol])/xval_range
                            if (ylim_vals[1] == 'log'):
                                gridvec_y = np.log10(selected[-1,ycol]/selected[-4,ycol])/yval_range
                            else:
                                gridvec_y = (selected[-1,ycol]-selected[-4,ycol])/yval_range
                                
                        gridvec_xn = gridvec_x/np.sqrt(gridvec_x**2+gridvec_y**2)
                        gridvec_yn = gridvec_y/np.sqrt(gridvec_x**2+gridvec_y**2)
                    
                        if (xlim_vals[1] == 'log'):
                            xtxt = selected[ival,xcol]*10**((txt_off_scale[j]*gridvec_xn+txt_off[j,0])*xval_range)
                        else:
                            xtxt = selected[ival,xcol]+(txt_off_scale[j]*gridvec_xn+txt_off[j,0])*xval_range
                        if (ylim_vals[1] == 'log'):
                            ytxt = selected[ival,ycol]*10**((txt_off_scale[j]*gridvec_yn+txt_off[j,1])*yval_range)
                        else:
                            ytxt = selected[ival,ycol]+(txt_off_scale[j]*gridvec_yn+txt_off[j,1])*yval_range
                        plt.text(xtxt,ytxt,str(round(angle))+"$^{\circ}$",color=colour_val,ha=txt_hjust[j],va=txt_vjust[j])
 
                
    if legend_on == True:
        plt.legend(fontsize=10,loc=legend_loc,handlelength=1,handletextpad=0.5)
    plt.title(fig_title, fontsize=14)
    plt.xlabel(xlabel_val,fontsize=12)
    plt.ylabel(ylabel_val,fontsize=12)
    plt.xlim(xlim_vals[0])
    plt.ylim(ylim_vals[0])
    plt.xscale(xlim_vals[1])
    plt.yscale(ylim_vals[1])
    if figfile != None:
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()
    return


def plot_stparams_grid_old(input_lists,xcol,ycol,xlim_vals,ylim_vals,xlabel_val,ylabel_val,legend_on,legend_loc,fig_title,
                       figfile):
    plt.figure(figsize=(6,5))
    for inlist in input_lists:
        geom_pars = inlist[0]
        txt_list = inlist[1]
        colour_val = inlist[2]
        hcor_list = inlist[4]
        angle_list = inlist[5]
        if txt_list != None:
            txt_vals = txt_list[0]
            txt_off = txt_list[1]
            txt_hjust = txt_list[2]
            txt_vjust = txt_list[3]
        # If spherical then plot single line
        if hcor_list == None:
            label_val = inlist[3]
            plt.plot(geom_pars[:,xcol],geom_pars[:,ycol],linewidth=3,color=colour_val,label=label_val)
            if txt_list != None:
                for i, radius in enumerate(txt_vals):
                    selected = geom_pars[geom_pars[:,0] == radius]
                    plt.text(selected[0,xcol]+txt_off[i,0],selected[0,ycol]+txt_off[i,1],
                             str(int(round(radius)))+" $R_{g}$",color=colour_val,ha=txt_hjust[i],va=txt_vjust[i])
        else:
            n = len(geom_pars)//2
            # First the constant hcor lines
            for i, hcor in enumerate(hcor_list):
                geom_pars2 = geom_pars[:n,:]
                selected = geom_pars2[geom_pars2[:,1] == hcor]
                if i == 0:
                    label_val = inlist[3]
                else:
                    label_val = None
                plt.plot(selected[:,xcol],selected[:,ycol],linewidth=3,color=colour_val,label=label_val)

                if txt_list != None:
                    if txt_vals[i] != 'off':
                        if txt_vals[i] == 'lower':
                            plt.text(selected[0,xcol]+txt_off[i,0],selected[0,ycol]+txt_off[i,1],
                                 str(round(hcor))+" $R_{g}$",color=colour_val,ha=txt_hjust[i],va=txt_vjust[i])
                        if txt_vals[i] == 'upper':
                            plt.text(selected[-1,xcol]+txt_off[i,0],selected[-1,ycol]+txt_off[i,1],
                                 str(round(hcor))+" $R_{g}$",color=colour_val,ha=txt_hjust[i],va=txt_vjust[i])
                            
            # Now the constant cone angle lines
            for i, angle in enumerate(angle_list):
                geom_pars2 = geom_pars[n:,:]
                selected = geom_pars2[geom_pars2[:,2] == angle]
                plt.plot(selected[:,xcol],selected[:,ycol],linewidth=3,color=colour_val,label=None)
                if txt_list != None:
                    j = i + len(hcor_list)
                    if txt_vals[j] != 'off':
                        if txt_vals[j] == 'lower':
                            plt.text(selected[0,xcol]+txt_off[j,0],selected[0,ycol]+txt_off[j,1],
                                 str(round(angle))+"$^{\circ}$",color=colour_val,ha=txt_hjust[j],va=txt_vjust[j])
                        if txt_vals[j] == 'upper':
                            plt.text(selected[-1,xcol]+txt_off[j,0],selected[-1,ycol]+txt_off[j,1],
                                 str(round(angle))+"$^{\circ}$",color=colour_val,ha=txt_hjust[j],va=txt_vjust[j])
                
    if legend_on == True:
        plt.legend(fontsize=10,loc=legend_loc,handlelength=1,handletextpad=0.5)
    plt.title(fig_title, fontsize=14)
    plt.xlabel(xlabel_val,fontsize=12)
    plt.ylabel(ylabel_val,fontsize=12)
    plt.xlim(xlim_vals[0])
    plt.ylim(ylim_vals[0])
    plt.xscale(xlim_vals[1])
    plt.yscale(ylim_vals[1])
    if figfile != None:
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()
    return


def compare_irfs(rin,rout,nrad,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,geom_list,
                    thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac,
                    figfile,leg_title):
    
    minbin = 1.0
    plt.figure()
    for i, x in enumerate(geom_list):
        cor_geometry = x[0]
        geopar = x[1]
        rcor = geopar[0]
        colour_val = x[2]
        label_val = x[3]
        sigfreq = x[4]
        rcor_col = x[5]
    
        print("*******************************************************************")
        print("Calculating for geometrical parameters: ",geopar)
        print("*******************************************************************")
        gamma_mean, rad, tau, i_rcor, i_rsigmax, deltau, ldisk_disp, lseed_disp,\
        lheat, ldisk_rev, lseed_rev, flux_irf_mono, lfreq = run_setup_mono(rin,rout,nrad,rcor,
                    freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac)
        
        disk_irf = ldisk_disp + ldisk_rev
        seed_irf = lseed_disp + lseed_rev
        heat_irf = lheat
        tau_edges, disk_irf_rb = irf_plot_rebin(disk_irf,deltau,minbin,i_rcor,i_rsigmax)
        tau_edges, seed_irf_rb = irf_plot_rebin(seed_irf,deltau,minbin,i_rcor,i_rsigmax)
        tau_edges, heat_irf_rb = irf_plot_rebin(heat_irf,deltau,minbin,i_rcor,i_rsigmax)
        
        dummy_data = (tau_edges[1:] + tau_edges[:-1])/2.
        dtau = -0.5*(tau_edges[:-1]+tau_edges[1:])
        disk_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*disk_irf_rb,
                    histtype='step',linestyle='solid',linewidth=2,color=colour_val,label=label_val)
        seed_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*seed_irf_rb,
                    histtype='step',linestyle='dotted',linewidth=2,color=colour_val)
        heat_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*heat_irf_rb,
                    histtype='step',linestyle='dashed',linewidth=2,color=colour_val)
        if (np.sum(deltau[:i_rcor]) == 0.):
            plt.scatter(np.mean(tau_edges[-2:]),disk_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='solid',color=colour_val)
            plt.scatter(np.mean(tau_edges[-2:]),seed_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='dotted',color=colour_val)
            plt.scatter(np.mean(tau_edges[-2:]),heat_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='dashed',color=colour_val)

        if sigfreq != None:
            for freq_radius in sigfreq:
                if freq_radius == 0:
                    radius = rcor
                    if (1/(tau[i_rcor]*t_scale)) >= 10:
                        freqtxt = str(round(1/(tau[i_rcor]*t_scale)))+" Hz"
                    elif (1/(tau[i_rcor]*t_scale)) < 1:
                        freqtxt = str(round(1/(tau[i_rcor]*t_scale),2))+" Hz"
                    else:
                        freqtxt = str(round(1/(tau[i_rcor]*t_scale),1))+" Hz"
                else:
                    radius = (freq_radius*t_scale*disk_tau_par[0]*10**(-1.*disk_tau_par[1]))**(-1/(1.5+disk_tau_par[1]))
                    freqtxt = str(freq_radius)+" Hz"
                rsig, i_rsig = find_nearest(rad,radius)
                tausig = np.sum(deltau[:i_rsig])
                plt.axvline(-1.*tausig,color='lightgray',linestyle='dashed',linewidth=2)
                radtxt = str(round(radius))+" $R_{g}$"
                plt.text(-1.*tausig,0.43,radtxt,fontsize=12,horizontalalignment='center')
                plt.text(-1.*tausig,0.34,freqtxt,fontsize=10,horizontalalignment='center')
        if (rcor_col != None):
            if rcor_col == False:
                rcor_colour = 'gray'
            else:
                rcor_colour = colour_val
            tausig = np.sum(deltau[:i_rcor])
            plt.axvline(-1.*tausig,color=rcor_colour,linestyle='dashdot',alpha=0.5,linewidth=2)
            radtxt = str(round(rcor))+" $R_{g}$"
            if (1/(tau[i_rcor]*t_scale)) >= 10:
                freqtxt = str(round(1/(tau[i_rcor]*t_scale)))+" Hz"
            elif (1/(tau[i_rcor]*t_scale)) < 1:
                freqtxt = str(round(1/(tau[i_rcor]*t_scale),2))+" Hz"
            else:
                freqtxt = str(round(1/(tau[i_rcor]*t_scale),1))+" Hz"
            plt.text(-1.*tausig,0.43,radtxt,color=rcor_colour,fontsize=12,horizontalalignment='center')
            plt.text(-1.*tausig,0.34,freqtxt,color=rcor_colour,fontsize=10,horizontalalignment='center')
            
        
        
    plt.yscale('log')
    plt.xscale('symlog', linthresh=100., linscale=0.5)
    plt.xlim(-1*np.sum(deltau[i_rcor:i_rsigmax+1]),20.)
    plt.ylim(2e-3,0.3)
    plt.xticks([-1e5,-1e4,-1e3,-1e2],fontsize=11)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='x',labelsize=12, which='both', direction='in')
    plt.tick_params(axis='y',labelsize=12, which='both', direction='in', right=True)
    plt.xlabel('Delay ($R_g/c$)', fontsize=13)
    plt.ylabel(r'Response $\times c/R_{g} \times$ Delay', fontsize=13)
    plt.legend(fontsize=10,loc='upper left',title=leg_title)
    if figfile != None:
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()
    return
    
    
def irf_plot_rebin(irf,deltau,minbin,i_rcor,i_rsigmax):
    if (np.sum(deltau[:i_rcor]) == 0.):
        tau_edges = np.zeros((i_rsigmax-i_rcor)+3)
        irf_rb = np.zeros(len(tau_edges)-1)
        tau_edges[-1:-3:-1] = [minbin,0.]
        tau_edges[-3::-1] = -1.*np.cumsum(deltau[i_rcor:i_rsigmax+1])
        irf_rb[-1] = np.sum(irf[0:i_rcor])/minbin
        irf_rb[-2::-1] = irf[i_rcor:i_rsigmax+1]/np.abs(np.diff(tau_edges[-2::-1]))
    else:
        tau_edges = np.zeros(i_rsigmax+2)
        irf_rb = np.zeros(len(tau_edges)-1)
        tau_edges[-2::-1] = -1.*np.cumsum(deltau[:i_rsigmax+1])
        irf_rb = irf[i_rsigmax::-1]/np.diff(tau_edges)
    return tau_edges, irf_rb


def run_setup_mono(rin,rout,nrad,rcor,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac):
    
    rsigmax = (freq_rsigmax*t_scale*disk_tau_par[0]*10**(-1.*disk_tau_par[1]))**(-1/(1.5+disk_tau_par[1]))
    print("Minimum frequency =",freq_rsigmax,"Hz and maximum signal radius = ",rsigmax,"Rg")
    rout = max([2*rsigmax,rout])
        
    #First set up radial grid
    rad_edge1 = np.logspace(np.log10(rin), np.log10(rsigmax), nrad+1)  # set up radial bin edges
    rad_edge2 = np.logspace(np.log10(rsigmax), np.log10(rout), nrad+1)
    rad_edge = np.append(rad_edge1,rad_edge2[1:])
    nrad = nrad*2
    
    rad = np.sqrt(rad_edge[1:]*rad_edge[:-1])  # Define radial bin centres as geometric mean of bin edges
    rad_area = np.pi*(np.square(rad_edge[1:])-np.square(rad_edge[:-1]))
    # Now reset coronal outer radius (i.e. disk inner radius) to be at the nearest radial bin edges
    rcor, i_rcor = find_nearest(rad_edge,rcor)
    disp_frac, seed_frac_flow, heat_frac_flow = \
        calc_dispfrac(rad,rad_area,rin,rcor,seedff_norm,seedff_ind,heatff_norm,heatff_ind)
    # Normalised viscous dissipation and parameters for seed and heating fractions within coronal flow.
    # Reset maximum signal radius to nearest radial bin value
    rsigmax, i_rsigmax = find_nearest(rad,rsigmax)
    print("Coronal radius reset to nearest radial bin edge: ",rcor)
    print("Maximum signal radius reset to: ",rsigmax)
    # Calculate illumination of corona by disk and vice-versa
    omega_cor, disktocor_frac, cortodisk_frac = calc_illumination_fracs(rad,rad_area,cor_geometry,geopar)
    print(np.sum(disp_frac[i_rcor:]*disktocor_frac[i_rcor:])/np.sum(disp_frac[i_rcor:])," of the disk flux is intercepted by the corona")
    print(np.sum(cortodisk_frac[i_rcor:])," of the coronal flux is intercepted by the disk")
    # Now calculate radial dependence of timing parameters and calculate propagation delays
    tau, lfreq, q, rms = calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, disk_tau_par, cor_tau_par,\
                                        lor_model, lor_par)
    deltau = calc_propagation_params(rad,rad_edge,rcor,disk_prop_par,cor_prop_par)
    # Calculate the radially-dependent direct disk, seed photon and coronal heating luminosities
    ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev =\
        calc_radial_time_response(rad,rcor,disp_frac,disktocor_frac,cortodisk_frac,thermal_frac,seed_frac_flow,heat_frac_flow)
    print("Seed energy :",e_seed)
    print("Dissipation disk luminosity: ",np.sum(ldisk_disp))
    print("Dissipation seed luminosity: ",np.sum(lseed_disp))
    print("Coronal heating luminosity: ",np.sum(lheat))
    print("Disk luminosity due to heating by corona: ",np.sum(ldisk_rev))
    print("Seed luminosity due to heating by corona: ",np.sum(lseed_rev))
    print("Total luminosity: ",np.sum(ldisk_disp+ldisk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lseed_disp+lheat+lseed_rev/thermal_frac)))

        # Calculate IRFs for the mono-energetic case
    gamma_mean, gamma_irf, flux_irf_mono, disk_irf_mono, seed_irf_mono =\
                    calc_irfs_mono(gamma_par,e_seed,ens_mono,ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev)
    
    print("Mean gamma is:", gamma_mean)
    
    return gamma_mean, rad, tau, i_rcor, i_rsigmax, deltau, ldisk_disp,\
    lseed_disp, lheat, ldisk_rev, lseed_rev, flux_irf_mono, lfreq

    

def calc_irf_centroid(irf,deltau):
    # Calculates IRF centroid lag seen by each radial bin
    irftau = np.where(deltau > 0, (np.cumsum(deltau)-deltau/2)*irf,np.cumsum(deltau)*irf)
    centroid = np.where(np.cumsum(irf) > 0, np.cumsum(irftau)/np.cumsum(irf),0.)
    return centroid


def compare_disk_seed_centroids(rin,rout,nrad,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,geom_list,
                    t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac,
                    figfile,leg_title,sigrad,ylim_vals,xlim_vals):
        
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,8),sharex=True)
    fig.subplots_adjust(hspace=0)
    for i, x in enumerate(geom_list):
        cor_geometry = x[0]
        geopar = x[1]
        rcor = geopar[0]
        thermal_frac = x[2]
        colour_val = x[3]
        label_val = x[4]
        lsty_val = x[5]
        
        print("*******************************************************************")
        print("Calculating for geometrical parameters: ",geopar)
        print("*******************************************************************")
        gamma_mean, rad, tau, i_rcor, i_rsigmax, deltau, ldisk_disp,\
            lseed_disp, lheat, ldisk_rev, lseed_rev, flux_irf, lfreq = \
        run_setup_mono(rin,rout,nrad,rcor,freq_rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par,gamma_par,e_seed,ens_mono,encomb,nirf_mult,min_deltau_frac)
                        
        deltau_scale = deltau*t_scale
        disk_irf = ldisk_disp + ldisk_rev
        seed_irf = lseed_disp + lseed_rev
        disk_cent = calc_irf_centroid(disk_irf,deltau_scale)
        seed_cent = calc_irf_centroid(seed_irf,deltau_scale)
        flux_cent = calc_irf_centroid(flux_irf[1,:],deltau_scale)
        gamma_ind = gamma_par[1]
        u = gamma_mean*gamma_ind*(np.log(ens_mono[1]/e_seed)-(gamma_mean-1)**(-1))
        lfrac = np.where(np.cumsum(seed_irf) > 0., 1./(np.cumsum(seed_irf)/np.sum(seed_irf)), 0.)
        if i < len(geom_list)//2:
            ax1.plot(lfreq,lfreq*(disk_cent-flux_cent),linewidth=2,color=colour_val,
                 linestyle=lsty_val,label=label_val)
        else:
            ax2.plot(lfreq,lfreq*(disk_cent-flux_cent),linewidth=2,color=colour_val,
                 linestyle=lsty_val,label=label_val)

    for i, ax in enumerate((ax1, ax2)):
        ax.axhline(color='grey',linestyle='dotted')
        if sigrad != None:
            for radius in sigrad:
                rsig, i_rsig = find_nearest(rad,radius)
                tausig = np.sum(deltau_scale[:i_rsig])
                ax.axvline(lfreq[i_rsig],color='gray',linestyle='dashed',linewidth=2)
                if ax == ax1:
                    radtxt = str(round(radius))+" $R_{g}$"
                    ax.text(lfreq[i_rsig],ylim_vals[1]*1.05,radtxt,fontsize=12,horizontalalignment='center')
            ax.set_ylabel(r'$\nu_{\rm signal} \times \left(\langle \tau_{\rm soft}\rangle - \langle \tau_{\rm d}\rangle \right)$', fontsize=14)
            ax.tick_params(axis='x',labelsize=12, which='both', direction='in', top=True)
            ax.tick_params(axis='y',labelsize=12, which='both', direction='in', right=True)
            ax.set_ylim(*ylim_vals)
            ax.legend(fontsize=11,title=leg_title[i])
    ax2.set_xlabel(r'$\nu_{\rm fluc}$ (Hz)',fontsize=14)
    ax2.set_xscale('log')
    ax2.set_xlim(*xlim_vals)
    if figfile != None:
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()
    return



def setup_params(rin,rout,nrad,rcor,rsigmax,seedff_norm,seedff_ind,heatff_norm,heatff_ind,cor_geometry,
                    geopar,thermal_frac,t_scale,disk_tau_par,cor_tau_par,lor_model,lor_par,
                    disk_prop_par,cor_prop_par):

    print("#######################################")
    print("Setting up disk and corona parameters")
    #First set up radial grid
    rad_edge = np.logspace(np.log10(rin), np.log10(rout), nrad+1)  # set up radial bin edges
    rad = np.sqrt(rad_edge[1:]*rad_edge[:-1])  # Define radial bin centres as geometric mean of bin edges
    rad_area = np.pi*(np.square(rad_edge[1:])-np.square(rad_edge[:-1]))
    disp_frac, seed_frac_flow, heat_frac_flow = \
        calc_dispfrac(rad,rad_edge,rad_area,seedff_norm,seedff_ind,heatff_norm,heatff_ind) 
    # Normalised viscous dissipation and parameters for seed and heating fractions within coronal flow.
    # Now reset coronal outer radius (i.e. disk inner radius) to be at the nearest radial bin edges
    rcor, i_rcor = find_nearest(rad_edge,rcor)
    # Reset maximum signal radius to nearest radial bin value
    rsigmax, i_rsigmax = find_nearest(rad,rsigmax)
    print("Coronal radius reset to nearest radial bin edge: ",rcor)
    print("Maximum signal radius reset to: ",rsigmax)
    # Calculate illumination of corona by disk and vice-versa
    print("Coronal geometry = ",str(cor_geometry)," with parameters ",geopar)
    omega_cor, disktocor_frac, cortodisk_frac = calc_illumination_fracs(rad,rad_area,cor_geometry,geopar)
    print(np.sum(disp_frac[i_rcor:]*disktocor_frac[i_rcor:])/np.sum(disp_frac[i_rcor:])," of the disk flux\
          is intercepted by the corona")
    print(np.sum(cortodisk_frac[i_rcor:])," of the coronal flux is intercepted by the disk")
    # Now calculate radial dependence of timing parameters and calculate propagation delays
    tau, lfreq, q, rms = calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, disk_tau_par, cor_tau_par,\
                                        lor_model, lor_par)
    deltau = calc_propagation_params(rad,rad_edge,rcor,i_rcor,disk_prop_par,cor_prop_par)
    # Calculate the radially-dependent direct disk, seed photon and coronal heating luminosities
    ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev =\
        calc_radial_time_response(rad,i_rcor,disp_frac,disktocor_frac,
                                      cortodisk_frac,thermal_frac,seed_frac_flow,heat_frac_flow)

    print("Dissipation disk luminosity: ",np.sum(ldisk_disp))
    print("Dissipation seed luminosity: ",np.sum(lseed_disp))
    print("Coronal heating luminosity: ",np.sum(lheat))
    print("Disk luminosity due to heating by corona: ",np.sum(ldisk_rev))
    print("Seed luminosity due to heating by corona: ",np.sum(lseed_rev))
    print("Total luminosity: ",np.sum(ldisk_disp+ldisk_rev/thermal_frac+(1-np.sum(cortodisk_frac))*(lseed_disp+lheat+lseed_rev/thermal_frac)))
    print("#######################################")    

    return rad, rad_edge, rad_area, deltau, i_rcor, i_rsigmax, tau, lfreq, q, rms, ldisk_disp, lseed_disp,\
        lheat, ldisk_rev, lseed_rev, disktocor_frac, cortodisk_frac


def setup_mono_irfs(rad,gamma_par,e_seed,ens_mono,ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev):

    print("#######################################")
    print("Calculating mono-energetic IRFs")
    # Calculate IRFs for the mono-energetic case
    gamma_mean, gamma_irf, flux_irf_mono, disk_irf_mono, seed_irf_mono =\
            calc_irfs_mono(rad,gamma_par,e_seed,ens_mono,ldisk_disp,lseed_disp,lheat,ldisk_rev,lseed_rev)
    print("Mean gamma is:", gamma_mean)
    print("#######################################")
    return gamma_mean, gamma_irf, flux_irf_mono, disk_irf_mono, seed_irf_mono

def plot_irfs(input_list,figfile,leg_title,sigrad):
    
    minbin = 1.0
    plt.figure()
    for x in input_list:
        [irfs_list,rad,deltau,i_rcor,i_rsigmax,colour_val,label_val] = x
        [disk_irf,seed_irf,heat_irf] = irfs_list
        tau_edges, disk_irf_rb = irf_plot_rebin(disk_irf,deltau,minbin,i_rcor,i_rsigmax)
        tau_edges, seed_irf_rb = irf_plot_rebin(seed_irf,deltau,minbin,i_rcor,i_rsigmax)
        tau_edges, heat_irf_rb = irf_plot_rebin(heat_irf,deltau,minbin,i_rcor,i_rsigmax)
        dummy_data = (tau_edges[1:] + tau_edges[:-1])/2.
        dtau = -0.5*(tau_edges[:-1]+tau_edges[1:])
        disk_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*disk_irf_rb,
                    histtype='step',linestyle='solid',linewidth=2,color=colour_val,label=label_val)
        seed_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*seed_irf_rb,
                    histtype='step',linestyle='dotted',linewidth=2,color=colour_val)
        heat_bins, edges, patches = plt.hist(dummy_data, bins=tau_edges, weights=np.abs(dtau)*heat_irf_rb,
                    histtype='step',linestyle='dashed',linewidth=2,color=colour_val)
        if (np.sum(deltau[:i_rcor]) == 0.):
            plt.scatter(np.mean(tau_edges[-2:]),disk_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='solid',color=colour_val)
            plt.scatter(np.mean(tau_edges[-2:]),seed_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='dotted',color=colour_val)
            plt.scatter(np.mean(tau_edges[-2:]),heat_irf_rb[-1],marker='_',s=200,
                linewidth=3,linestyle='dashed',color=colour_val)
        
    plt.yscale('log')
    plt.xscale('symlog', linthresh=100., linscale=0.5)
    plt.xlim(-1.*np.sum(deltau[i_rcor:i_rsigmax+1]),20.)
    plt.ylim(2e-3,0.3)
    plt.xticks([-1e5,-1e4,-1e3,-1e2],fontsize=11)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='x',labelsize=12, which='both', direction='in')
    plt.tick_params(axis='y',labelsize=12, which='both', direction='in', right=True)
    plt.xlabel('Delay ($R_g/c$)', fontsize=13)
    plt.ylabel(r'Response $\times c/R_{g} \times$ Delay', fontsize=13)
    if sigrad != None:
        for radius in sigrad:
            rsig, i_rsig = find_nearest(rad,radius)
            tausig = np.sum(deltau[:i_rsig])
            plt.axvline(-1.*tausig,color='gray',linestyle='dashed',linewidth=2)
            radtxt = str(round(radius))+" $R_{g}$"
            plt.text(-1.*tausig,0.34,radtxt,fontsize=12,horizontalalignment='center')

    plt.legend(fontsize=10,loc='upper left',title=leg_title)
    if figfile != None:
        plt.savefig(figfile,bbox_inches='tight')
    plt.show()
    return

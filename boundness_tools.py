# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:54:51 2019

@author: elcastel
"""
import numpy as np
from scipy.spatial import cKDTree 
import scipy.interpolate as interp
import scipy.stats as stats
from scipy.stats import anderson
import pandas as pd

"""
This file contains all of the functions needed to determine whether or not two groups are gravitationally bounded. Examples of using
these functions are shown in test_all_mock_boundness.py and test_part_mock_boundness.py. to use these functions, you need 
the following information:
    - czs, ras, decs, and masses of the two groups you are testing the boundness of. One of these groups will be
    the 'target' and one will be the 'group'. 
    - To use groups that are closest to one another within the catalog, use the nearest_neighbor function. the num_neighbors argument 
    can be the maximum group size + 1, so that every galaxy has a neighbor that is outside its group. 
    - get_target is a function that will return all of the info needed to test boundness for the target group. the function takes
    df, which is the entire dataframe you're working with; tarname, which is the galaxy ID of the target galaxy (one galaxy at a time
    is the target galaxy, and every galaxy in that galaxy's group is in the target group.); neighbor_dist and neighbor_ind, which are 
    the results of the nearest_neighbor function.
"""

def get_skycoords(cz, ra, dec):
    """gets x, y, z coordinates of each galaxy in group and central cz of a gruop. 
    Inputs are arrays of members czs, ras decs
    
    Uses conversion from spherical to cartesian to get x,y, and z coordinates. Good for use in nearest neighbor search! 
    """
    # convert ra and dec to radians
    ra = ra * np.pi / 180
    dec = dec * np.pi / 180
    ra = np.array(ra)
    dec = np.array(dec)

    # define variables
    H0 = 70.
    theta = ra
    phi = np.pi/2 - dec
    r = cz/H0 #Hubble distance

    # convert from spherical coordinates to cartesian
    xvals = r *  np.sin(phi) * np.cos(theta)
    yvals = r * np.sin(phi) * np.sin(theta)
    zvals = r * np.cos(phi)
    czcen = np.mean(cz)
    return xvals, yvals, zvals, czcen
    


    
    
def nearest_neighbor(ras, decs, czs, masses, num_neighbors):
    """Takes coordinates of the galaxies you are testing and the number of 
    nearest neighbors you desire. Returns distance and indices to the num_neighbors
    nearest neighbors for each galaxy in the group. Note that if the galaxy's nearest
    neighbor is not included in the group being tested, it can't be found. Neighbors
    are all within the group."""

    #coordinates in Mpc, includes peculiar motion
    xvals, yvals, zvals, czcen = get_skycoords(czs, ras, decs)
    
    #Find 3D neighbors
    coords = np.array([xvals, yvals, zvals]).T #transpose to make ordered pair coords
    kdt = cKDTree(coords)
    nn_dists, nn_inds = kdt.query(coords, k=num_neighbors) #gives distances and indices of the kth closest neighbors 

    nn_dists = nn_dists[:,1:] #delete first column which is a self match
    nn_inds = nn_inds[:,1:]
        
    return nn_dists, nn_inds

def get_proj_distance(ras1, decs1, ras2, decs2, czs1, czs2):
    H0 = 100. #h
    ra1 = np.array(np.mean(ras1))* np.pi/180
    dec1 = np.array(np.mean(decs1))* np.pi/180
    ra2 = np.array(np.mean(ras2))* np.pi/180
    dec2 = np.array(np.mean(decs2))* np.pi/180
    cz1 = np.array(np.mean(czs1))
    cz2 = np.array(np.mean(czs2))
    
    theta = 2*np.arcsin(np.sqrt((np.sin((dec2-dec1)/2))**2 + np.cos(dec1)*np.cos(dec2)*(np.sin((ra2-ra1)/2)**2)))
    distance = 1/H0 * (cz1 + cz2) * np.sin(theta/2)
    return distance

def getmhoffset(delta1, delta2, borc1, borc2, cc):
    """Adapted from Katie's code, using eqns from "Sample Variance Considerations for Cluster Surveys," Hu & Kravtsov (2003) ApJ, 584, 702
    (astro-ph/0203169)
    delta1 is overdensity of input, delta2 is overdensity of output -- for mock, delta1 = 200. for RESOLVE, delta1=280
    borc = 1 if wrt background density, borc = 0 if wrt critical density
    cc is concentration of halo- use cc=6
    """
    if borc1 == 0:
        delta1 = delta1/0.3
    if borc2 == 0:
        delta2 = delta2/0.3
    
    xin = 1./cc
    f_1overc = (xin)**3. * (np.log(1. + (1./xin)) - (1. + xin)**(-1.))
    f1 = delta1/delta2 * f_1overc
    a1=0.5116
    a2=-0.4283
    a3=-3.13e-3
    a4=-3.52e-5
    p = a2 + a3*np.log(f1) + a4*(np.log(f1))**2.
    x_f1 = (a1*f1**(2.*p) + (3./4.)**2.)**(-1./2.) + 2.*f1
    r2overr1 = x_f1
    m1overm2 = (delta1/delta2) * (1./r2overr1)**3. * (1./cc)**3.
    
    return m1overm2

# =============================================================================
# def Rvir(logmh, h, delta_mean=337):
#     rhocrit = 2.77e11*h**2. #Msun/Mpc^3
#     omega_m = 0.3
#     rhomean = omega_m*rhocrit
#     Rvir = ((3.0*10.**(logmh)/(4.*np.pi*delta_mean*rhomean))**(1./3.))
#     return Rvir
# =============================================================================
  
def virvelocity(logmh):
  G = 4.32e-9
  M = 10**logmh
  rvir = Rvir(logmh, 1., 337.)
  cvir = 11. * (M/4e12)**(-0.13)
  Ac = np.log(1+cvir) - cvir/(1+cvir)
  Vmax2 = (G*M/rvir)*0.216*cvir/Ac
  return np.sqrt(Vmax2)

def circle(r, centerx=0, centery=0):
    theta = np.linspace(0, 2*np.pi, 1000)
    x = centerx + r * np.cos(theta)
    y = centery + r * np.sin(theta)
    return x,y
  
def angular_size(R, czcen, H0=100):
  theta = R*H0*180/czcen/np.pi
  return theta

    
def angular_separation(ra1, dec1, ra2, dec2):
    """Compute the angular separation bewteen two galaxies using the Haversine formula.
    Arguments: 
      ra1, dec1, ra2, dec2 (float): each sky coordinate in decimal degrees
    Returns: angular separation of two coordinate pairs in the sky in radians (type float)
    """
    
    phi1 = ra1 * np.pi/180.
    theta1 = np.pi/2. - (dec1 * np.pi/180.)
    phi2 = ra2 * np.pi/180.
    theta2 = np.pi/2. - (dec2 * np.pi/180.)
    yangle = 2*np.arcsin(np.sqrt(np.sin((theta2-theta1)/2.0)**2.0 + np.sin(theta2)*np.sin(theta1)*np.sin((phi2-phi1)/2.0)**2.0))
    havtheta = np.float128(yangle)
    return havtheta
  
def ellipse_coords(centerx, centery, lenx, leny):
  """
  ----
  Draw an ellipse of a specified major/minor axis size at a given plotting center.
  Arguments:
    centerx (type float): center coordinate of ellipse on horizontal axis.
    centery (typle float): center coordinate of ellipse on vertical axis.
    lenx (type float): semi-length of ellipse in x 
    leny (type float): semi-length of ellipse in y
  Returns:
    xvals (iterable): x-coordinates of ellipse
    yvals (iterable): y-coordinates of ellipse
  ----
  """
  angle = np.linspace(0,2*np.pi,3000)
  xvals = [(centerx + lenx*np.cos(theta)) for theta in angle]
  yvals = [(centery + leny*np.sin(theta)) for theta in angle]
  return xvals, yvals
  
# =============================================================================
# def Rproj(grp0, h):
#     """Adapted from Katie's IDL code, which was used to calculate grprproj for FoF groups in RESOLVE"""
#     ras = np.array(grp0.radeg)* np.pi/180 #in radians
#     decs = np.array(grp0.dedeg)* np.pi/180 #in radians
#     czs = np.array(grp0.cz)
#     grpn = np.shape(grp0)[0] 
#     H0 = 100*h #km/s/Mpc
#     grpdec = np.mean(decs)
#     grpra = np.mean(ras)
#     grpcz = np.mean(czs)
#     theta = 2*np.arcsin(np.sqrt((np.sin((decs-grpdec)/2))**2 + np.cos(decs)*np.cos(grpdec)*(np.sin((ras-grpra)/2)**2)))
#     theta = np.array(theta)
#     rproj = theta*grpcz/H0
#     sortorder = np.argsort(rproj)
#     rprojsrt = rproj[sortorder]
#     rprojval = np.arange(0 , grpn)/(grpn-1.) #array from 0 to 1, use for interpolation
#     #print(rprojsrt, rprojval, grpn)
#     f = interp.interp1d(rprojval, rprojsrt) 
#     rproj75 = f(0.75)
#     return rproj75
# =============================================================================
def Rvir(logmh, h, delta_mean=337.):
    """From equations listed in RESOLVE glossary by logMh entry"""
    rhocrit = 2.77e11*h**2. #Msun/Mpc^3
    omega_m = 0.3
    rhomean = omega_m*rhocrit
    Rvir = ((3.0*10.**(logmh)/(4.*np.pi*delta_mean*rhomean))**(1./3.))
    return Rvir

def Rproj(grp0, h, ras=0, decs=0, czs=0):
    """Adapted from Katie's IDL code, which was used to calculate grprproj for FoF groups in RESOLVE. The try/except is to allow this code to work for the mock catalog or another dataset where the variable names aren't radeg, dedeg, and czs like they are in RESOLVE and ECO. For these, use optional arguments to define variable names"""
    if np.shape(grp0)[0] > 2:
        try:
      	    ras = np.array(grp0.radeg)* np.pi/180 #in radians
            decs = np.array(grp0.dedeg)* np.pi/180 #in radians
            czs = np.array(grp0.cz)
        except AttributeError:
            ras = np.array(ras)
            decs = np.array(decs)
            czs = np.array(czs)
            #print('helllllooooo')
            #print(type(ras), type(decs), type(czs))
        grpn = np.shape(grp0)[0] 
        H0 = 100*h #km/s/Mpc
        grpdec = np.mean(decs)
        grpra = np.mean(ras)
        grpcz = np.mean(czs)
        theta = 2*np.arcsin(np.sqrt((np.sin((decs-grpdec)/2))**2 + np.cos(decs)*np.cos(grpdec)*(np.sin((ras-grpra)/2)**2)))
        theta = np.array(theta)
        rproj = theta*grpcz/H0
        sortorder = np.argsort(rproj)
        rprojsrt_y = rproj[sortorder]
        rprojval_x = np.arange(0 , grpn)/(grpn-1.) #array from 0 to 1, use for interpolation
        f = interp.interp1d(rprojval_x, rprojsrt_y) 
        rproj75 = f(0.75)
    else:
        rproj75 = 0.0
        rprojval_x = 0.0
        rprojsrt_y = 0.0
    return rproj75

def calculate_compactness(grp, h):
    if np.shape(grp)[0] > 2:
        logmh_grp337 = np.array(grp.logmh_bound)[0]
        #mh280overmh337 = getmhoffset(280., 337., 1., 1., 6. )
        #logmh_grp337 = np.log10(10.**logmh_grp280/mh280overmh337)
        Rproj_grp = Rproj(grp, h)
        Rvir_grp = Rvir(logmh_grp337, h)
        compactness = Rvir_grp / Rproj_grp
    else:
        return 0.0
    return compactness
    
def calculate_gas_content(grp):
    gas = 10**np.array(grp.logmgas)
    gas = np.sum(gas)
    stars = 10**np.array(grp.logmstar)
    stars = np.sum(stars)
    gastostellar = gas/stars
    gastostellar = np.log10(gastostellar)
    return gastostellar

def color_gap(grp):
    if np.array(grp.grpn_bound)[0] > 1:
        gal_mags = np.array(grp.absrmag)
        sort = np.argsort(gal_mags)
        mags_sort = gal_mags[sort] #sorts from smallest to largest magnitude, meaning brightest to dimmest 
        brightest_mag = mags_sort[0]
        second_brightest_mag = mags_sort[1]
        brightest_gal = grp.loc[grp.absrmag == brightest_mag]
        second_brightest_gal = grp.loc[grp.absrmag == second_brightest_mag]
        brightest_ur = np.array(brightest_gal.modelu_r)[0]
        second_brightest_ur = np.array(second_brightest_gal.modelu_r)[0]
        colorgap = brightest_ur - second_brightest_ur
        return colorgap
    
def AD_test(grp):
    czs = grp.cz
    czs = np.array(czs)
    number        = 5            #min number of galaxies in a group
    n = int(np.shape(czs)[0])
    if n > number:
        #put observed redshifts in increasing order for each group
        sort = np.argsort(czs)
        czs = czs[sort]
            
        #do AD test
        andtest = anderson(czs,dist='norm')
        A2notstar = andtest[0]
        A2 = A2notstar*(1 + 0.75/n + 2.25/(n**2))
            
        a = 3.6789468
        b = 0.1749916  #from Hou et al 2009
        alpha = a * np.exp(-A2/b)  #alpha < 5% not gaussian- only 1 galaxy not gaussian
        
        return alpha
    
def urcolor(grp):
     ur = np.array(grp.modelu_r)
     avgur = np.mean(ur)
     return avgur
 
def sfr(grp):
    sfr = np.array(grp.sfr_nuv)
    sfravg = np.mean(sfr)
    return sfravg

def large_scale_environment(grp):
    try:
        dens = np.array(grp.den1mpc)
        avgdens = np.mean(dens)
        if avgdens > 0:
            return avgdens
        else:
            return 'no env data'
    except AttributeError:
        return 'no env data'
      
def get_nn_grp_info(df, nn_ind, boundorfof):
    """Set up nearest neighbor galaxy and get info on group.The 'boundorfof' is 'bound' if you want
    to test boundness to the nearest bound group and 'fof' if you want to test boundness to the nearest fof group."""
    nngal = df.iloc[nn_ind]

    if boundorfof == 'fof':
        nn_grpID = np.array(nngal.grp)
        group = df.loc[df.grp == nn_grpID]
        Ngrp = max(np.array(group.grpn))
        
        ras_grp = group.radeg
        decs_grp = group.dedeg
        czs_grp = group.cz
    
        masses_grp = 10**group.logmstar + 10**group.logmgas + (10**group.logmh)/Ngrp
        
        ras_grp = np.array(ras_grp)
        decs_grp = np.array(decs_grp)
        czs_grp = np.array(czs_grp)
        masses_grp = np.array(masses_grp)
    elif boundorfof == 'bound':
        nn_grpID = np.array(nngal.grp_bound)
        group = df.loc[df.grp_bound == nn_grpID]
        Ngrp = max(np.array(group.grpn_bound))
        
        ras_grp = group.radeg
        decs_grp = group.dedeg
        czs_grp = group.cz
    
        masses_grp = 10**group.logmstar + 10**group.logmgas + (10**group.logmh)/Ngrp
        
        ras_grp = np.array(ras_grp)
        decs_grp = np.array(decs_grp)
        czs_grp = np.array(czs_grp)
        masses_grp = np.array(masses_grp)       
    
    return group, nn_grpID, ras_grp, decs_grp, czs_grp, masses_grp, Ngrp

def get_target(df, tar_name, neighbor_dist, neighbor_ind):
    """Set up target group and get info on nearest neighbor of target galaxy. 
    One galaxy at a time is the targal. This function returns info on nearest neighbor
    of that galaxy, and info on the coordinates, mass, etc. of the group that the targal is in.
    The group can have one or more members."""
    df_indOrig = df.reset_index()
    tar_ind = df_indOrig.loc[df_indOrig['name']==tar_name].index[0] #use in .iloc
    

    targal = df.loc[df.name == tar_name]
    
    tar_grpID = np.array(targal.grp)

    df_grpIDs = df.set_index('grp')    
    targals = df_grpIDs.loc[tar_grpID ]
    tar_ngrp = np.shape(targals)[0]
    ra_tars = targals.radeg
    dec_tars = targals.dedeg
    cz_tars = targals.cz

    mass_tars = 10**targals.logmstar + 10**targals.logmgas + (10**targals.logmh)/tar_ngrp
    
    ra_tars = np.array(ra_tars)
    dec_tars = np.array(dec_tars)
    cz_tars = np.array(cz_tars)
    mass_tars = np.array(mass_tars)
    tar_grpID = np.array(tar_grpID)
    
    nn_dist = neighbor_dist[tar_ind]
    nn_ind = neighbor_ind[tar_ind] 

    return targal, targals, tar_grpID, tar_ngrp, ra_tars, dec_tars, cz_tars, mass_tars, nn_dist, nn_ind

def inverse_transform_sampling(data, n_bins=40, n_samples=1000):
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interp.interp1d(cum_values, bin_edges)
    r = np.random.rand(n_samples)
    return inv_cdf(r)

def test_boundness(ra_tars, dec_tars, cz_tars, mass_tars, ras_grp, decs_grp, czs_grp, masses_grp, proj_effect_corr = 'yes'):
    """This function takes the coordinates and masses of the two groups whose boundness is being tested.
    Returns 'yes' if groups are bounded and 'no' if they are not. """

    raddist = pd.read_csv('eco_3Ddistanceratios_duplicatesmarked_102220_wlogmh_grpn.csv')
    raddist = raddist.loc[raddist.duplicate_distance == 'no']

    veldist = pd.read_csv('eco_3Ddistance_vel_ratios_duplicatesmarked_102620_wlogmh_grpn.csv')
    veldist = veldist.loc[veldist.duplicate_ratio == 'no']

    rads = raddist.R3D_over_Rproj_obs
    vels = veldist.V_real_over_V_obs
    
    czcen_tar = np.mean(cz_tars)
    cz_cen = np.mean(czs_grp)
    G = 4.301e-9 #km^2/s^2 Mpc/Msun


    """Calculate escape velocity from point mass located at the center of mass at the position of the target.
    Distance is projected distance onto the sky (2D distance in x and y)"""
    distance = get_proj_distance(ra_tars, dec_tars, ras_grp, decs_grp, cz_tars, czs_grp)
    vgrpgrp_los = np.abs(cz_cen - czcen_tar)
    
    """If a ratio of R_grp-grp (3D) / R_proj (2D) is given, use to correct for projection effects but multiplying distance*ratio"""
    if proj_effect_corr == 'yes':
        rgrpgrp_dist = rads * distance # rads = 3D/2D distance distribution, distance = 2D distance
        vesc_dist = np.sqrt(2*G*np.sum(masses_grp)/rgrpgrp_dist)
        vgrpgrp_dist = vels * vgrpgrp_los # vels = 3D/los velocity distribution, vgrpgrp_los = los rel. vel. 
        vesc_grp = np.mean(vesc_dist)
        vgrpgrp = np.mean(vgrpgrp_dist)
        nsamples = 10000
        vgrpgrp_sample = inverse_transform_sampling(vgrpgrp_dist, n_samples=nsamples, n_bins=6000) # vpec
        vesc_sample = inverse_transform_sampling(vesc_dist, n_samples=nsamples) # vesc
        # get probability numerically
        numboundprob =  np.sum(np.float128(np.greater_equal(vesc_sample, vgrpgrp_sample)))/np.float128(nsamples)
    else: 
        vesc_grp = np.sqrt(2*G*np.sum(masses_grp)/distance)

        """Find the peculiar velocity of the target. Account for unknonwn direction by reducing the magnitude by sqrt(3) """
        vgrpgrp = vgrpgrp_los * np.sqrt(3)


    """Test if target halo is bounded to halo(s)"""
    bounded = 0
    not_bounded = 0
    ratio = vgrpgrp / vesc_grp
    
    # if vesc_grp > vgrpgrp:
    if numboundprob > 0.9:
        bounded += 1
        bound = 'yes'
    else:
        bound = 'no'
        not_bounded += 1

    return bound, numboundprob, vgrpgrp

"""def test_boundness_3d(ra_tars, dec_tars, distobs_tar, czobs_tars, czreal_tars, vptan_tars, mass_tars, ras_grp, decs_grp, distobs_grp, czsobs_grp, czsreal_grp, vptans_grp, masses_grp):
        print(vesc_grp, np.mean(vesc_dist))
        print(vgrpgrp, np.mean(vgrpgrp_dist))
    #This function takes the coordinates and masses of the two groups whose boundness is being tested.
    #Returns 'yes' if groups are bounded and 'no' if they are not. 

    czcenobs_tar = np.mean(czobs_tars)
    cz_cenobs = np.mean(czsobs_grp)
    czcenreal_tar = np.mean(czreal_tars)
    cz_cenreal = np.mean(czsreal_grp)

    rdistcen_tar = np.mean(distobs_tar)
    rdist_cen = np.mean(distobs_grp)
    
    vptancen_tar = np.mean(vptan_tars)
    vptan_cen = np.mean(vptans_grp)

    G = 4.301e-9 #km^2/s^2 Mpc/Msun

    #Get coordinates of target
#    cz_tars = np.array([cz_tars]).T
#    ra_tars = np.array([ra_tars]).T
#    dec_tars = np.array([dec_tars]).T

    #Calculate escape velocity from point mass located at the center of mass at the position of the target.
    #Distance is projected distance onto the sky (2D distance in x and y)
    distanceonsky = get_proj_distance(ra_tars, dec_tars, ras_grp, decs_grp, czreal_tars, czsreal_grp)
    distancelos = np.abs(rdistcen_tar - rdist_cen)
    distance = np.sqrt(distanceonsky**2 + distancelos**2)

    vesc_grp = np.sqrt(2*G*np.sum(masses_grp)/distance)

    #print('Escape velocity from group: ' + str(vesc_grp) + ' km/s')

    #Find the peculiar velocity of the target. 
    vzgrpgrp = np.abs(czcenobs_tar - czcenreal_tar) - (cz_cenobs - cz_cenreal)
    vptangrpgrp = np.abs(vptancen_tar - vptan_cen)
    vgrpgrp = np.sqrt(vzgrpgrp**2 + vptangrpgrp**2)

    #Test if target group is bounded to neighbor group
    bounded = 0
    not_bounded = 0
    ratio = vgrpgrp / vesc_grp
    if vesc_grp > vgrpgrp:

        bounded += 1
        bound = 'yes'

    else:
        bound = 'no'
        not_bounded += 1

    return bound, ratio, vgrpgrp"""

def test_boundness_3d(ra_tars, dec_tars, distobs_tar, czreal_tars, mass_tars, vx_tars, vy_tars, vz_tars, ras_grp, decs_grp, distobs_grp, czsreal_grp, masses_grp, vx_grp, vy_grp, vz_grp):
    """This function takes the coordinates and masses of the two groups whose boundness is being tested.
    Returns 'yes' if groups are bounded and 'no' if they are not. """

    czcenreal_tar = np.mean(czreal_tars)
    cz_cenreal = np.mean(czsreal_grp)

    rdistcen_tar = np.mean(distobs_tar)
    rdist_cen = np.mean(distobs_grp)

    grpvx = np.mean(vx_tars)
    grpvy = np.mean(vy_tars)
    grpvz = np.mean(vz_tars)
    
    nnvx = np.mean(vx_grp)
    nnvy = np.mean(vy_grp)
    nnvz = np.mean(vz_grp)

    G = 4.301e-9 #km^2/s^2 Mpc/Msun

    """Get coordinates of target"""
#    cz_tars = np.array([cz_tars]).T
#    ra_tars = np.array([ra_tars]).T
#    dec_tars = np.array([dec_tars]).T

    """Calculate escape velocity from point mass located at the center of mass at the position of the target.
    Distance is projected distance onto the sky (2D distance in x and y)"""
    distanceonsky = get_proj_distance(ra_tars, dec_tars, ras_grp, decs_grp, czreal_tars, czsreal_grp)
    distancelos = np.abs(rdistcen_tar - rdist_cen)
    distance = np.sqrt(distanceonsky**2 + distancelos**2)

    vesc_grp = np.sqrt(2*G*np.sum(masses_grp)/distance)

    #print('Escape velocity from group: ' + str(vesc_grp) + ' km/s')

    """Find the peculiar velocity of the target. """
    vgrpgrp = np.sqrt((grpvx - nnvx)**2. + (grpvy - nnvy)**2 + (grpvz - nnvz)**2)

    """Test if target group is bounded to neighbor group"""
    bounded = 0
    not_bounded = 0
    ratio = vgrpgrp / vesc_grp
    if vesc_grp > vgrpgrp:

        bounded += 1
        bound = 'yes'

    else:
        bound = 'no'
        not_bounded += 1

    return bound, ratio, vgrpgrp

def crossing_time(grp):
    """
    grp is pandas df of group or bound multi-group system you want crossing time of
    <r> is the average projected distance of group members from
    the group centre of mass.
    <|v|> is the average speed of group members relative to the
    group centre of mass
    t_c = <r>/<|v|>
    """
    rprojs = []
    relvs = []
    try:
        grpra = np.mean(np.array(grp.radeg))
        grpdec = np.mean(np.array(grp.dedeg))
        grpcz = np.mean(np.array(grp.cz))
    except AttributeError:
        grpra = np.mean(np.array(grp.ra))
        grpdec = np.mean(np.array(grp.dec))
        grpcz = np.mean(np.array(grp.cz))
    for i in range(len(grp)):
        thisgal = grp.iloc[i]
        try:
            rproj = get_proj_distance(thisgal.radeg, thisgal.dedeg, grpra, grpdec, thisgal.cz, grpcz)
        except AttributeError:
            rproj = get_proj_distance(thisgal.ra, thisgal.dec, grpra, grpdec, thisgal.cz, grpcz)
        rprojs.append(rproj*1e3) #in kpc
        speed = np.array(grp.cz.iloc[i])
        relv = np.abs(speed - grpcz)
        if relv == 0:
            relv += 0.1
        relvs.append(relv)

    if np.mean(relvs) > 0:
        t_c = np.mean(rprojs)/np.mean(relvs)
        if t_c != None:
            return t_c # in Gyr, multiply by 1e9 to get in years 
        elif t_c == None:
            print('uh oh, none')

def color_gap(grp):
    if np.array(grp.boundN)[0] > 1:
        """ gal_urcolor = np.array(grp.modelu_r)
        colordisp = np.std(gal_urcolor)
        return colordisp"""
        print(grp)
        gal_mags = np.array(grp.absrmag)
        sort = np.argsort(gal_mags)
        mags_sort = gal_mags[sort] #sorts from smallest to largest magnitude, meaning brightest to dimmest 
        brightest_mag = mags_sort[0]
        second_brightest_mag = mags_sort[1]
        try:
            brightest_gal = grp.loc[grp.absrmag == brightest_mag]
            second_brightest_gal = grp.loc[grp.absrmag == second_brightest_mag]
        except AttributeError:
            brightest_gal = grp.loc[grp.M_r == brightest_mag]
            second_brightest_gal = grp.loc[grp.M_r == second_brightest_mag]
        try:
            brightest_ur = np.array(brightest_gal.modelu_r)[0]
            second_brightest_ur = np.array(second_brightest_gal.modelu_r)[0]
        except AttributeError:
            brightest_ur = np.array(brightest_gal.u_r)[0]
            second_brightest_ur = np.array(second_brightest_gal.u_r)[0]
        colorgap = brightest_ur - second_brightest_ur
        return colorgap
    else:
        return 0

def AD_test(grp):
    from scipy.stats import anderson
    czs = grp.cz
    czs = np.array(czs)
    number        = 5            #min number of galaxies in a group
    n = int(np.shape(czs)[0])
    if n > number:
        #put observed redshifts in increasing order for each group
        sort = np.argsort(czs)
        czs = czs[sort]

        #do AD test
        andtest = anderson(czs,dist='norm')
        A2notstar = andtest[0]
        A2 = A2notstar*(1 + 0.75/n + 2.25/(n**2))

        a = 3.6789468
        b = 0.1749916  #from Hou et al 2009
        alpha = a * np.exp(-A2/b)  #alpha < 5% not gaussian- only 1 galaxy not gaussian

        return alpha
    else:
        return 0.
                                                                       
def calculate_gas_content(grp):
    try:
        gas = 10**np.array(grp.logmgas)
    except AttributeError:
        gas = np.array(grp.mhi)
    gas = np.sum(gas)
    stars = 10**np.array(grp.logmstar)
    stars = np.sum(stars)
    stars = np.log10(stars)
    gas = np.log10(gas)
    return gas, stars

from sklearn.cluster import KMeans
def scrambled(orig):
    dest = orig[:]
    np.random.shuffle(dest)
    return dest

def DS_test(grp):
    czs = grp.cz
    try:
        ras = grp.radeg
        decs = grp.dedeg
    except AttributeError:
        ras = grp.ra
        decs = grp.dec
    H0 = 70 #h
    ras = np.array(ras)
    decs = np.array(decs)
    czs = np.array(czs)
    del_sim = np.array([[0.]*100]*len(czs))
    delta = [0.] * len(czs)
    #DELSIM= [0.] * len(czs)
    #p_valueuntitled2 = [0.] * len(czs)
    n = int(np.shape(czs)[0])
    #print(n)
    size=11
    if n >= size:
        #loop through each galaxy in group
        for j in range(0,len(ras)):
            angles = np.sqrt(((ras - ras[j])*np.cos((np.pi/180.)*decs[j]))**2 + (decs - decs[j])**2)
            dists = (np.pi/180.)* angles * czs[j] / H0
            close = np.argsort(dists)[1:size] #indices of closest galaxies
            sig = np.sqrt(sum((czs - np.mean(czs))**2)/(n - 1))
            czs_loc = czs[close]    #czs of closest gals
            sig_loc = np.sqrt(sum((czs_loc - np.mean(czs_loc))**2)/(n - 1))
            delta[j] = (size/sig**2) * ((np.mean(czs_loc) - np.mean(czs))**2 + (sig_loc - sig)**2)
            #simulated values
            for k in range(0,100):
                cz_scram = np.array(scrambled(list(czs)))
                angles  = np.sqrt(((ras - ras[j])*np.cos((np.pi/180.)*decs[j]))**2 + (decs - decs[j])**2)
                dists = (np.pi/180.)* angles * cz_scram[j] / H0
                close = np.argsort(dists)[1:size] #indices of closest galaxies
                sig = np.sqrt(sum((cz_scram - np.mean(cz_scram))**2)/(n - 1))
                czs_loc = cz_scram[close]    #czs of closest gals
                sig_loc = np.sqrt(sum((czs_loc - np.mean(czs_loc))**2)/(n - 1))
                del_sim[j][k] = (size/sig**2) * ((np.mean(czs_loc) - np.mean(czs))**2 + (sig_loc - sig)**2)
        DELTA = sum(np.sqrt(delta[0:n]))
        DELSIM = sum(np.sqrt(del_sim[0:n]))
        p_value = len(np.where(DELSIM > DELTA)[0])/100.
        return p_value

def r337_overlap(grp1, grp2, h):

    offset = np.log10(h)
    try:
        Rvir1 = np.array(grp1.grpR337)[0]
        Rvir2 = np.array(grp2.grpR337)[0]
    except IndexError:
        try:
            Rvir1 = np.array(grp1.grpR337)
            Rvir2 = np.array(grp2.grpR337)
        except AttributeError:
            try:
                logmh3371 = np.array(grp1.logmh)[0]
                logmh3372 = np.array(grp2.logmh)[0]
            except IndexError:
                logmh3371 = np.array(grp1.logmh)
                logmh3372 = np.array(grp2.logmh)
            Rvir1 = bd.Rvir(logmh3371, h)
            Rvir2 = bd.Rvir(logmh3372, h)
    sumrvir = Rvir1 + Rvir2

    try:
        Rgrpgrp2d = get_proj_distance(grp1.radeg, grp1.dedeg, grp2.radeg, grp2.dedeg, grp1.cz, grp2.cz)
    except AttributeError:
        Rgrpgrp2d = get_proj_distance(grp1.radeg, grp1.dedeg, grp2.radeg, grp2.dedeg, grp1.cz_obs, grp2.cz_obs)

    raddist = pd.read_csv('eco_3Ddistanceratios_duplicatesmarked_102220_wlogmh_grpn.csv')
    raddist = raddist.loc[raddist.duplicate_distance == 'no']

    rads = raddist.R3D_over_Rproj_obs

    d = Rgrpgrp2d * rads
    meddist = np.median(d)

    if meddist > sumrvir:
        overlap = 0.
    else:
        overlap = 1

    return overlap

def r337overlap_inbdsystem(bdgrp):
    fofids = np.unique(np.array(bdgrp.grp))
    overlap = 0.
    for i in range(len(fofids)):
        fof1 = bdgrp.loc[bdgrp.grp == fofids[i]]
        for j in range(len(fofids)):
            if i != j:
                fof2 = bdgrp.loc[bdgrp.grp == fofids[j]]
                overlap = r337_overlap(fof1, fof2, 0.7)
                if overlap == 1:
                    return overlap
    return overlap



# *- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:55:32 2019

@author: elcastel
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 19:57:38 2019

@author: elcastel
"""

import numpy as np
import pandas as pd
import boundness_tools as bd

def make_grp_cat(gal_cat, twoDorthreeD, absrmag, mag_floor, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, czreal=None,  vx=None, vy=None, vz=None, rdist=None):
    galsabovefloor = gal_cat.loc[gal_cat[absrmag] < mag_floor]
    grpids = np.unique(np.array(galsabovefloor[grp]))

    if twoDorthreeD == '3d':
        grps = pd.DataFrame(columns = ['grp', 'grpn', 'radeg','dedeg', 'cz_obs', 'mass', 'logmh', 'nn_inds', 'cz_real', 'vx', 'vy', 'vz', 'rdist', 'bound', 'grpn_bound', 'grp_bound', 'logmh_bound', 'Rvir', 'Rvir_bound'])
    if twoDorthreeD == '2d':
        grps = pd.DataFrame(columns = ['grp', 'grpn', 'radeg','dedeg', 'cz_obs', 'mass', 'logmh', 'nn_inds', 'bound', 'grpn_bound', 'grp_bound', 'logmh_bound', 'Rvir', 'Rvir_bound'])

    grps.grp = grpids
    
    for i in range(len(grpids)):
        thisgrp = gal_cat.loc[(gal_cat[grp] == grpids[i]) & (gal_cat[absrmag] < mag_floor)]
        grpra = np.mean(thisgrp[radeg])
        grpdec = np.mean(thisgrp[dedeg])
        grpczobs = np.mean(thisgrp[cz])
        thisgrpn = np.shape(thisgrp)[0]
        grpmass = calc_grp_mass(thisgrp, gasmasslogged, logmstar, gasmass, logmh, grpn)
        grplogmh = np.log10(sum(10**np.unique(np.array(thisgrp[logmh]))))
        grpRvir = bd.Rvir(grplogmh, 0.7)

        grps.loc[grps.grp == grpids[i], 'radeg'] = grpra
        grps.loc[grps.grp == grpids[i], 'dedeg'] = grpdec
        grps.loc[grps.grp == grpids[i], 'cz_obs'] = grpczobs
        grps.loc[grps.grp == grpids[i], 'mass'] = grpmass   
        grps.loc[grps.grp == grpids[i], 'logmh'] = grplogmh
        grps.loc[grps.grp == grpids[i], 'grpn'] = thisgrpn
        grps.loc[grps.grp == grpids[i], 'Rvir'] = grpRvir
        
        if twoDorthreeD == '3d':
            grpczreal = np.mean(thisgrp[czreal])
            grprdist = np.mean(thisgrp[rdist])
            grpvx = np.mean(thisgrp[vx])
            grpvy = np.mean(thisgrp[vy])
            grpvz = np.mean(thisgrp[vz])

            grps.loc[grps.grp == grpids[i], 'cz_real'] = grpczreal
            grps.loc[grps.grp == grpids[i], 'rdist'] = grprdist
            grps.loc[grps.grp == grpids[i], 'vx'] = grpvx
            grps.loc[grps.grp == grpids[i], 'vy'] = grpvy
            grps.loc[grps.grp == grpids[i], 'vz'] = grpvz

    ras = np.float128(np.array(grps.radeg))
    decs = np.float128(np.array(grps.dedeg))
    if twoDorthreeD == '3d':
        czs = np.float128(np.array(grps.cz_real))
    elif twoDorthreeD == '2d':
        czs = np.float128(np.array(grps.cz_obs))
    num_neighbors = 40
    neighbor_dist, neighbor_ind = bd.nearest_neighbor(ras, decs, czs, grpmass, num_neighbors)

    
    grpids = np.unique(np.array(grps.grp))
    for i in range(len(grpids)):
        grps.iat[i, 7] = neighbor_ind[i, :]
        
    grps.bound = 0
    grps.grpn_bound = grps.grpn
    grps.grp_bound = grps.grp
    grps.logmh_bound = grps.logmh
    grps.Rvir_bound = grps.Rvir
    
    return grps

def make_gal_cat(grp_cat, gal_cat, radeg, dedeg, cz, mag_floor, absrmag, grp, twoDorthreeD):
    bdids = np.array(grp_cat.grp_bound)
    for i in range(len(bdids)):
        bdid = bdids[i]
        bdgrp = grp_cat.loc[grp_cat.grp_bound == bdid]
        thisbound = np.array(bdgrp.bound)[0]
        thisgrpn_bound = np.array(bdgrp.grpn_bound)[0]
        thisgrp_bound = np.array(bdgrp.grp_bound)[0]
        thislogmh_bound = np.array(bdgrp.logmh_bound)[0]
        thisRvir = np.array(bdgrp.Rvir)[0]
        thisRvir_bound = np.array(bdgrp.Rvir_bound)[0]
        bdgrp_fofids = np.array(grp_cat.loc[grp_cat.grp_bound == bdid].grp)
        for j in range(len(bdgrp_fofids)):
            thisfofid = bdgrp_fofids[j]
            gal_cat.loc[(gal_cat[grp] == thisfofid) & (gal_cat[absrmag] < mag_floor), 'grp_bound'] = bdid
        bdgrp = gal_cat.loc[gal_cat.grp_bound == bdid]
        ras = bdgrp[radeg] * np.pi/180
        decs = bdgrp[dedeg] * np.pi/180
        czs = bdgrp[cz] * np.pi/180
        thisRproj_bound = bd.Rproj(bdgrp, 0.7, ras, decs, czs)
        
        gal_cat.loc[gal_cat.grp_bound == bdid, 'bound'] = thisbound
        gal_cat.loc[gal_cat.grp_bound == bdid, 'grpn_bound'] = thisgrpn_bound
        gal_cat.loc[gal_cat.grp_bound == bdid, 'logmh_bound'] = thislogmh_bound
        gal_cat.loc[gal_cat.grp_bound == bdid, 'Rvir'] = thisRvir
        gal_cat.loc[gal_cat.grp_bound == bdid, 'Rvir_bound'] = thisRvir_bound
        gal_cat.loc[gal_cat.grp_bound == bdid, 'Rproj_bound'] = thisRproj_bound
    return gal_cat

def calc_grp_mass(grp, gasmasslogged, logmstar, gasmass, logmh, grpn):
    if gasmasslogged == 'yes':
        grpmass = sum(10**grp[logmstar] + 10**grp[gasmass] + (10**grp[logmh]/np.array(grp[grpn])[0]))
    elif gasmasslogged == 'no':
        grpmass = sum(10**grp[logmstar] + grp[gasmass] + (10**grp[logmh]/np.array(grp[grpn])[0]))
    return grpmass

def get_gal_cat_colnames_twoDorthreeD(catpath, mag_floor,absrmag,radeg,dedeg,cz,logmstar,gasmass, gasmasslogged,logmh,grpn,grp,name, twoDorthreeD):
    """absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, and name are the names of the columns in the dataframe that
    correspond to each piece of information needed for boundness testing"""
    gal_cat = pd.read_csv(catpath) 
    #gal_cat = gal_cat.loc[gal_cat.absrmag < mag_floor]

    """Create new columns for properties of bound systems"""
    bound = [0]*np.shape(gal_cat)[0]
    gal_cat['bound'] = bound
    gal_cat['grpn_bound'] = gal_cat[grpn]
    gal_cat['grp_bound'] = gal_cat[grp]
    gal_cat['logmh_bound'] = gal_cat[logmh]
    gal_cat['Rvir'] = 0
    gal_cat['Rproj_bound'] = 0
    gal_cat['Rvir_bound'] = 0
    return gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD

    
def main():
    c = 3.0E5 #km/s

    """Read in FoF catalog"""
    ecopath = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/fofcats_61721/ECO_G3groupcatalog_052821.csv"
    mag_floor_eco  = -17.33 #ECO
    
    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(ecopath, mag_floor_eco,'absrmag','radeg','dedeg','cz','logmstar','logmgas','yes','logmh','grpn','grp','name', '2d')
    
    resolve  = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/fofcats_61721/RESOLVE_G3groupcatalog_052821.csv"
    mag_floor_resolve = -17.0 #RESOLVE
    
    
    # gal_cat,mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass, gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(resolve, mag_floor_resolve,'absrmag','radeg','dedeg','cz','logmstar','logmgas', 'yes','logmh','grpn','grp','name', '2d')

    mock = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/fofcats_61721/ECO_cat_2_Planck_memb_cat.csv"
    # mock = "/afs/cas.unc.edu/users/e/l/elcastel/MockCatalog/csvs/eco_m200_0.csv"
    mag_floor_mock = -17.33
    
    gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','M_group','g_ngal','groupid','g_galtype','2d')

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','M_group','g_ngal','groupid','g_galtype','3d')

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','loghalom','halo_ngal','haloid','g_galtype','3d')

    if twoDorthreeD == '3d':
        czreal, vx, vy, vz, rdist = 'cz_nodist', 'vx', 'vy', 'vz', 'r_dist'

    if twoDorthreeD == '2d':
        df = make_grp_cat(gal_cat, twoDorthreeD, absrmag, mag_floor, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp)
    if twoDorthreeD == '3d':  
        df = make_grp_cat(gal_cat, twoDorthreeD, absrmag, mag_floor, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, czreal=czreal, vx=vx, vy=vy, vz=vz, rdist=rdist)

    newID = max(np.array(df.grp)) + 1
    
    grpids = np.array(df.grp)
    
    """loop through each group"""
    for i in range(len(grpids)):
        """Get group and neighbor group"""
        thisgrp = df.loc[df.grp == grpids[i]]
    
        which_neighbor = 0
        thisgrp_nninds = thisgrp.nn_inds.iloc[0]

        neighborgrp_id = df.iloc[thisgrp_nninds[which_neighbor]].grp
        neighborgrp = df.loc[df.grp == neighborgrp_id]

        """Make sure group and neighbor group are not already in a bound multi-group system together"""
        while np.array(thisgrp.grp_bound)[0] == np.array(neighborgrp.grp_bound)[0]:
            which_neighbor += 1
            neighborgrp_id = df.iloc[thisgrp_nninds[which_neighbor]].grp
            neighborgrp = df.loc[df.grp == neighborgrp_id]

        """Test for boundness between group and neighbor group"""
        if twoDorthreeD == '2d':
            bound, prob, vgrpgrp = bd.test_boundness(thisgrp.radeg, thisgrp.dedeg, thisgrp.cz_obs, thisgrp.mass,neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.cz_obs, neighborgrp.mass)
        elif twoDorthreeD == '3d':
            bound, ratio, vgrpgrp = bd.test_boundness_3d(thisgrp.radeg, thisgrp.dedeg, thisgrp.rdist, thisgrp.cz_real, thisgrp.mass, thisgrp.vx, thisgrp.vy, thisgrp.vz, neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.rdist, neighborgrp.cz_real, neighborgrp.mass, neighborgrp.vx, neighborgrp.vy, neighborgrp.vz)

        """If the groups are bound, change 'bound' 'grpn_bound' and 'grp_bound' parameters and continue testing for boundness with increasingly distant neighbors until the original FoF group is not bound to it's kth nearest neighbor"""
        while bound == 'yes':
            df.loc[df.grp_bound == np.array(thisgrp.grp_bound)[0], 'bound'] = 1
            df.loc[df.grp_bound == np.array(thisgrp.grp_bound)[0], 'grp_bound'] = newID
    
            df.loc[df.grp_bound == np.array(neighborgrp.grp_bound)[0], 'bound'] = 1
            df.loc[df.grp_bound == np.array(neighborgrp.grp_bound)[0], 'grp_bound'] = newID
    
            thisgrp = df.loc[df.grp == np.array(thisgrp.grp)[0]]
    
            which_neighbor += 1
    
            neighborgrp_id = df.iloc[thisgrp_nninds[which_neighbor]].grp
            neighborgrp = df.loc[df.grp == neighborgrp_id]    
   
            if twoDorthreeD == '2d':
                bound, prob, vgrpgrp = bd.test_boundness(thisgrp.radeg, thisgrp.dedeg, thisgrp.cz_obs, thisgrp.mass,neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.cz_obs, neighborgrp.mass)
            elif twoDorthreeD == '3d':
                bound, ratio, vgrpgrp = bd.test_boundness_3d(thisgrp.radeg, thisgrp.dedeg, thisgrp.rdist, thisgrp.cz_real, thisgrp.mass, thisgrp.vx, thisgrp.vy, thisgrp.vz, neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.rdist, neighborgrp.cz_real, neighborgrp.mass, neighborgrp.vx, neighborgrp.vy, neighborgrp.vz)
        if i % 100 == 0:
            print(i)
        newID += 1
    bdids = np.unique(np.array(df.grp_bound))
    for i in range(len(bdids)):
        bdid = bdids[i]
        bdgrp = df.loc[df.grp_bound == bdid]
        this_grpn_bound = sum(np.array(bdgrp.grpn))
        logmhgrp = np.log10(sum(10**np.array(bdgrp.logmh)))
        Rvir_bdgrp = bd.Rvir(logmhgrp, 0.7)
        if i % 100 == 0:
            print(sum(np.array(bdgrp.grpn)))
            print(bdgrp.grpn)
            print(this_grpn_bound)
        df.loc[df.grp_bound == bdid, 'logmh_bound'] = logmhgrp
        df.loc[df.grp_bound == bdid, 'Rvir_bound'] = Rvir_bdgrp
        df.loc[df.grp_bound == bdid, 'grpn_bound'] = np.float128(this_grpn_bound)
            
    df.to_csv('grpcat_mock2d2_m337_62821_1.csv')

    gal_cat = make_gal_cat(df, gal_cat, radeg, dedeg, cz, mag_floor, absrmag, grp, twoDorthreeD)
    gal_cat.to_csv('galcat_mock2d2_m337_62821_1.csv')
   
main()

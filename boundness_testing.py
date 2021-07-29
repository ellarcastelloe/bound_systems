import numpy as np
import pandas as pd
import boundness_tools as bd

def make_grp_cat(gal_cat, twoDorthreeD, absrmag, mag_floor, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, czreal=None,  vx=None, vy=None, vz=None, rdist=None):
    if mag_floor == -17.0:
        galsabovefloor = gal_cat.loc[gal_cat.fl_insample == 1]
    else:
        galsabovefloor = gal_cat.loc[gal_cat[absrmag] < mag_floor]
    grpids = np.unique(np.array(galsabovefloor[grp]))

    if twoDorthreeD == '3d':
        grps = pd.DataFrame(columns = ['grp', 'grpn', 'radeg','dedeg', 'cz_obs', 'mass', 'logmh', 'nn_inds', 'cz_real', 'vx', 'vy', 'vz', 'rdist', 'boundFlag', 'boundN', 'boundID', 'boundLog337', 'grpR337'])
    if twoDorthreeD == '2d':
        grps = pd.DataFrame(columns = ['grp', 'grpn', 'radeg','dedeg', 'cz_obs', 'mass', 'logmh', 'nn_inds', 'boundFlag', 'boundN', 'boundID', 'boundLog337', 'grpR337'])

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
        grps.loc[grps.grp == grpids[i], 'grpR337'] = grpRvir

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

    grps.boundFlag = 0
    grps.boundN = grps.grpn
    grps.boundID = grps.grp
    grps.boundLog337 = grps.logmh

    return grps

def make_gal_cat(grp_cat, gal_cat, radeg, dedeg, cz, mag_floor, absrmag, grp, twoDorthreeD):
    bdids = np.array(grp_cat.boundID)
    for i in range(len(bdids)):
        bdid = bdids[i]
        bdgrp = grp_cat.loc[grp_cat.boundID == bdid]
        thisbound = np.array(bdgrp.boundFlag)[0]
        thisboundn = np.array(bdgrp.boundN)[0]
        thisbdid = np.array(bdgrp.boundID)[0]
        thislogmh_bound = np.array(bdgrp.boundLog337)[0]
        thisRvir = np.array(bdgrp.grpR337)[0]
        bdgrp_fofids = np.array(grp_cat.loc[grp_cat.boundID == bdid].grp)
        for j in range(len(bdgrp_fofids)):
            thisfofid = bdgrp_fofids[j]
            gal_cat.loc[(gal_cat[grp] == thisfofid) & (gal_cat[absrmag] < mag_floor), 'boundID'] = bdid
        bdgrp = gal_cat.loc[gal_cat.boundID == bdid]
        ras = bdgrp[radeg] * np.pi/180
        decs = bdgrp[dedeg] * np.pi/180
        czs = bdgrp[cz] * np.pi/180
        thisRproj_bound = bd.Rproj(bdgrp, 0.7, ras, decs, czs)
        thistc = bd.crossing_time(bdgrp)
        thiscolorgap = bd.color_gap(bdgrp)
        thislogg, thislogs = bd.calculate_gas_content(bdgrp)
        thisR337overlap = bd.r337overlap_inbdsystem(bdgrp)
        if thisboundn > 5:
            thisadalpha = bd.AD_test(bdgrp)
            gal_cat.loc[gal_cat.boundID == bdid, 'boundADalpha'] = thisadalpha
        if thisboundn > 11:
            thisdspval = bd.DS_test(bdgrp)
            gal_cat.loc[gal_cat.boundID == bdid, 'boundDSpval'] = thisdspval
        gal_cat.loc[gal_cat.boundID == bdid, 'boundLogG'] = thislogg
        gal_cat.loc[gal_cat.boundID == bdid, 'boundLogS'] = thislogs
        gal_cat.loc[gal_cat.boundID == bdid, 'boundURcolorgap'] = thiscolorgap
        gal_cat.loc[gal_cat.boundID == bdid, 'boundTCross'] = thistc
        gal_cat.loc[gal_cat.boundID == bdid, 'boundFlag'] = thisbound
        gal_cat.loc[gal_cat.boundID == bdid, 'boundN'] = thisboundn
        gal_cat.loc[gal_cat.boundID == bdid, 'boundLog337'] = thislogmh_bound
        gal_cat.loc[gal_cat.boundID == bdid, 'grpR337'] = thisRvir
        gal_cat.loc[gal_cat.boundID == bdid, 'boundRproj'] = thisRproj_bound
        gal_cat.loc[gal_cat.boundID == bdid, 'boundR337overlap'] = thisR337overlap
    return gal_cat

def calc_grp_mass(grp, gasmasslogged, logmstar, gasmass, logmh, grpn):
    if gasmasslogged == 'yes':
        grpmass = sum(10**grp[logmstar] + 10**grp[gasmass] + (10**grp[logmh]/np.array(grp[grpn])[0]))
    elif gasmasslogged == 'no':
        grpmass = sum(10**grp[logmstar] + grp[gasmass] + (10**grp[logmh]/np.array(grp[grpn])[0]))
    return grpmass
def get_gal_cat_colnames_twoDorthreeD(catpath, mag_floor,absrmag,radeg,dedeg,cz,logmstar,gasmass, gasmasslogged,logmh,grpn,grp,name,urcolor, twoDorthreeD):
    """absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, and name are the names of the columns in the dataframe that
    correspond to each piece of information needed for boundness testing"""
    if isinstance(catpath,  str) == True:
        gal_cat = pd.read_csv(catpath)
    else:
        gal_cat = catpath
    #gal_cat = gal_cat.loc[gal_cat.absrmag < mag_floor]
    """Create new columns for properties of bound systems"""
    bound = [0]*np.shape(gal_cat)[0]
    gal_cat['boundFlag'] = bound
    gal_cat['boundN'] = gal_cat[grpn]
    gal_cat['boundID'] = gal_cat[grp]
    gal_cat['boundLog337'] = gal_cat[logmh]
    gal_cat['grpR337'] = 0
    gal_cat['boundRproj'] = 0
    gal_cat['boundADalpha'] = 0.
    gal_cat['boundTCross'] = 0.
    gal_cat['boundURcolorgap'] = 0.
    gal_cat['boundDSpval'] = 0.
    gal_cat['boundLogG'] = 0.
    gal_cat['boundLogS'] = 0.
    gal_cat['boundR337overlap'] = 0.
    return gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name,urcolor, twoDorthreeD


def main():
    c = 3.0E5 #km/s

    """Read in FoF catalog"""
    ecopath = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/eco_basecat_wG3grps072321.csv"
    mag_floor_eco  = -17.33 #ECO

    gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, urcolor, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(ecopath, mag_floor_eco,'absrmag','radeg','dedeg','cz','logmstar','logmgas','yes','logmh','grpn','grp','name', 'modelu_r', '2d')

    resolve  = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/resolve_basecat_wG3grps072321.csv"
    mag_floor_resolve = -17.0 #RESOLVE


    # gal_cat,mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass, gasmasslogged, logmh, grpn, grp, name, urcolor, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(resolve, mag_floor_resolve,'absrmag','radeg','dedeg','cz','logmstar','logmgas', 'yes','logmh','grpn','grp','name', 'modelu_r', '2d')

    mock = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/fofcats_61721/ECO_cat_7_Planck_memb_cat.csv"
    # mock = "/afs/cas.unc.edu/users/e/l/elcastel/MockCatalog/csvs/eco_m200_0.csv"
    mag_floor_mock = -17.33

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','M_group','g_ngal','groupid','g_galtype','2d')

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','M_group','g_ngal','groupid','g_galtype','3d')

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(mock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi', 'no','loghalom','halo_ngal','haloid','g_galtype','3d')

    zackmockpath = "/afs/cas.unc.edu/users/e/l/elcastel/ECO/boundness/ECO_cat_5_Planck_memb_cat_mvir_withG3groups.csv"
    zackmock = pd.read_csv(zackmockpath)
    zackmock['g3grpntot_l'] = zackmock['g3grpngi_l'] + zackmock['g3grpndw_l']

    # gal_cat, mag_floor, absrmag, radeg, dedeg, cz, logmstar, gasmass,gasmasslogged, logmh, grpn, grp, name, twoDorthreeD = get_gal_cat_colnames_twoDorthreeD(zackmock, mag_floor_mock,'M_r','ra','dec','cz','logmstar','mhi','no','g3logmh_l','g3grpntot_l','g3grp_l','g3fc_l','2d')


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
        while np.array(thisgrp.boundID)[0] == np.array(neighborgrp.boundID)[0]:
            which_neighbor += 1
            neighborgrp_id = df.iloc[thisgrp_nninds[which_neighbor]].grp
            neighborgrp = df.loc[df.grp == neighborgrp_id]

        """Test for boundness between group and neighbor group"""
        if twoDorthreeD == '2d':
            bound, prob, vgrpgrp = bd.test_boundness(thisgrp.radeg, thisgrp.dedeg, thisgrp.cz_obs, thisgrp.mass,neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.cz_obs, neighborgrp.mass)
        elif twoDorthreeD == '3d':
            bound, ratio, vgrpgrp = bd.test_boundness_3d(thisgrp.radeg, thisgrp.dedeg, thisgrp.rdist, thisgrp.cz_real, thisgrp.mass, thisgrp.vx, thisgrp.vy, thisgrp.vz, neighborgrp.radeg, neighborgrp.dedeg, neighborgrp.rdist, neighborgrp.cz_real, neighborgrp.mass, neighborgrp.vx, neighborgrp.vy, neighborgrp.vz)

        """If the groups are bound, change 'boundFlag' 'boundN' and 'boundID' parameters and continue testing for boundness with increasingly distant neighbors until the original FoF group is not bound to it's kth nearest neighbor"""
        while bound == 'yes':
            df.loc[df.boundID == np.array(thisgrp.boundID)[0], 'boundFlag'] = 1
            df.loc[df.boundID == np.array(thisgrp.boundID)[0], 'boundID'] = newID

            df.loc[df.boundID == np.array(neighborgrp.boundID)[0], 'boundFlag'] = 1
            df.loc[df.boundID == np.array(neighborgrp.boundID)[0], 'boundID'] = newID

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
    bdids = np.unique(np.array(df.boundID))
    df['boundADalpha'] = 0.
    df['boundTCross'] = 0.
    df['boundURcolorgap'] = 0.
    df['boundLogS'] = 0.
    df['boundLogG'] = 0.
    df['boundDSpval'] = 0.
    df['boundR337overlap'] = 0.
    for i in range(len(bdids)):
        bdid = bdids[i]
        bdgrp = df.loc[df.boundID == bdid]
        thisboundn = sum(np.array(bdgrp.grpn))
        logmhgrp = np.log10(sum(10**np.array(bdgrp.logmh)))
        df.loc[df.boundID == bdid, 'boundLog337'] = logmhgrp
        df.loc[df.boundID == bdid, 'boundN'] = np.float128(thisboundn)
    df = df.rename(columns={"logmh": "log337"})


    dimdf = gal_cat.loc[gal_cat[absrmag] > mag_floor]

    cols = ['boundN', 'boundID', 'boundLog337','grpR337', 'boundRproj', 'modelu_r', 'den1mpc', 'boundADalpha', 'boundTCross','boundURcolorgap', 'log337', 'boundLogS', 'boundLogG', 'boundDSpval']

    dimdf['boundFlag'] = 0

    for i in range(len(cols)):
        thiscol = cols[i]
        dimdf[thiscol] = -99.
    gal_cat = make_gal_cat(df, gal_cat, radeg, dedeg, cz, mag_floor, absrmag, grp, twoDorthreeD)

    for i in range(len(bdids)):
        bdid = bdids[i]
        bdgrp = df.loc[df.boundID == bdid]
        adalpha = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundADalpha)[0]
        tcross = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundTCross)[0]
        colorgap = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundURcolorgap)[0]
        gas = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundLogG)[0]
        stars = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundLogS)[0]
        dspval = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundDSpval)[0]
        overlap = np.array(gal_cat.loc[gal_cat.boundID == bdid].boundR337overlap)[0]
        df.loc[df.boundID == bdid, 'boundADalpha'] = adalpha
        df.loc[df.boundID == bdid, 'boundTCross'] = tcross
        df.loc[df.boundID == bdid, 'boundURcolorgap'] = colorgap
        df.loc[df.boundID == bdid, 'boundLogS'] = stars
        df.loc[df.boundID == bdid, 'boundLogG'] = gas
        df.loc[df.boundID == bdid, 'boundDSpval'] = dspval
        df.loc[df.boundID == bdid, 'boundR337overlap'] = overlap
    df.to_csv('grpcat_eco_m337_72521_1.csv')

    gal_cat = pd.concat([gal_cat, dimdf])
    gal_cat.to_csv('galcat_eco_m337_72521_1.csv')

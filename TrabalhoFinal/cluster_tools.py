'''
Tools to analyze clusters using the GAIA DR3 data
Author: Mariane Dias - PhD Candidate at National Observatory
Rio de Janeiro - Brazil
'''
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from zero_point import zpt


def corr_plx(data,tb):

    """
    This function applies the correction needed in Parallax of GAIA DR3 data.
    Variables:
    - data: table of cluster data
    - tb: takes the values "gaia" or "vizier" which is the origin if the data table used.
    We use the method documented by Lindegren et al.(2021) https://doi.org/10.1051/0004-6361/202039653
    To apply this correction is necessary to install the module gaiadr3_zeropoint, which can be done with: 
    pip install gaiadr3-zeropoint

    """
    zpt.load_tables()
    if tb == "vizier":
        # in GAIA: phot_g_mean_mag --> in Vizier: Gmag
        gmag = data['Gmag'].data
        # nu_eff_used_in_astrometry --> nueff
        nueffused = data['nueff'].data
        # pseudocolour --> pscol
        psc = data['pscol'].data
        # ecl_lat --> ELAT
        ecl_lat = data['ELAT'].data
        # astrometric_params_solved --> Solved
        soltype = data['Solved'].data
    elif tb == "gaia":
        
        # in GAIA: phot_g_mean_mag --> in Vizier: Gmag
        gmag = data['phot_g_mean_mag'].data
        # nu_eff_used_in_astrometry --> nueff
        nueffused = data['nu_eff_used_in_astrometry'].data
        # pseudocolour --> pscol
        psc = data['pseudocolour'].data
        # ecl_lat --> ELAT
        ecl_lat = data['ecl_lat'].data
        # astrometric_params_solved --> Solved
        soltype = data['astrometric_params_solved'].data
    
    #This is made so the func get_zpt() don't fail (sources with 2-p solutions)
    valid = soltype>3
    zpvals = zpt.get_zpt(gmag[valid], nueffused[valid], psc[valid],
                 ecl_lat[valid], soltype[valid],_warnings=False)
    
    #Applying the correction and add a new column with the corrected values
    data = data[valid]
    corr_plx = data['Plx'] - zpvals
    data.add_column(corr_plx, name='Plx_corr')

    return data

def quality_filter(data,tb):
    
    """
    This function makes a filter in the data using the quantity RUWE (normalized error) that gives the quality of the data. 
    The data used in here will have ruwe<1.4, that will extract astrometric binaries.
        Variables:
        - data: table of cluster data
        - tb: takes the values "gaia" or "vizier" which is the origin if the data table used.   
    """
    if tb == "vizier":
        ruwe = data['RUWE']

    elif tb == "gaia":
        ruwe = data['ruwe']

    kk, = np.where(ruwe < 1.4)
    data = data[kk]
    return data

def movprop_filter(data,pmra_cluster,pmde_cluster, L, tb):
    data = data
    """
    This function takes the stars inside a box with size L centered in the known position of the cluster, which can be obtained in SIMBAD, for example.
    Variables:
        - data: table of cluster data
        - tb: takes the values "gaia" or "vizier" which is the origin if the data table used.
    - pmra_cluster = pmra of the cluster from SIMBAD (or other font)
        - pmde_cluster = pmra of the cluster from SIMBAD (or other font)   
    - L = size of the box centered in the cluster   
    """

    if tb == "vizier":
        pmra = data["pmRA"]
        pmde = data["pmDE"]

    elif tb == "gaia":
        pmra = data["pmra"]
        pmde = data["pmdec"]

    kk, = np.where(((pmra > pmra_cluster - L)&(pmra < pmra_cluster + L))&((pmde > pmde_cluster - L)&(pmde < pmde_cluster + L)))
    data = data[kk]
    return data

def movprop_stat_filter(data,pmra_cluster,pmde_cluster,sigma_pmra,sigma_pmde, window, tb):
    data = data
    """
    This function takes the stars inside a box with sizes window*sigma_pmra x window*sigma_pmde centered in the position of the cluster obtained by statistical data, which can be done with a previous analysis of te cluster.
    Variables:
        - data: table of cluster data
        - tb: takes the values "gaia" or "vizier" which is the origin if the data table used.
    - pmra_cluster = pmra of the cluster from previous analysis
    - pmde_cluster = pmra of the cluster from previous analysis   
    - sigma_pmra = standard deviation of the histogram of pmra made in previous analysis of the cluster
    - sigma_pmde = standard deviation of the histogram of pmde made in previous analysis
    - window = factor that we multiply deviations to make a box that encompass the cluster in proper motion
    """
    L_pmra = window*sigma_pmra
    L_pmde = window*sigma_pmde

    if tb == "vizier":
        pmra = data["pmRA"]
        pmde = data["pmDE"]

    elif tb == "gaia":
        pmra = data["pmra"]
        pmde = data["pmdec"]

    kk, = np.where(((pmra > pmra_cluster - L_pmra)&(pmra < pmra_cluster + L_pmra))&((pmde > pmde_cluster - L_pmde)&(pmde < pmde_cluster + L_pmde)))
    data = data[kk]
    return data

def quality2_filter(data):

    """
    This function makes a filter in the data using the quantity RUWE (normalized error) that gives the quality of the data. 
    The data used in here will have ruwe<1.4, that will extract astrometric binaries.
        Variables:
        - data: table of cluster data  
    """
   
    g = data['Gmag']
    ruwe = data['RUWE']
    phot_bp_rp = data["phot_bp_rp_excess_factor_corr"]	
    c0=0.0059898
    c1=8.817481E-12
    m=7.618399
    N=5
    sigma_c1=c0+c1*(g)**m
    C_ast1=np.sqrt(phot_bp_rp**2)

    kk, = np.where((ruwe < 1.4)&((C_ast1 - N*sigma_c1)<0.))
    data = data[kk]
    return data


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



	#Obtaining the values of the correction in Parallax
	zpvals = zpt.get_zpt(gmag, nueffused, psc, ecl_lat, soltype)

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
	The data used in here will have 1.0 <ruwe<1.4, that will extract astrometric binaries.
        Variables:
        - data: table of cluster data
        - tb: takes the values "gaia" or "vizier" which is the origin if the data table used.	
	"""
        if tb == "vizier":
               ruwe = data['RUWE'] 

        elif tb == "gaia":
               ruwe = data['ruwe']	
	
	kk, = np.where((ruwe > 1.0)&(ruwe < 1.4))
	
	return data[kk]


def movprop_filter(data,pmra_cluster,pmde_cluster, L, tb):	
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
	
	return data[kk]


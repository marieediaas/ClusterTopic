'''
	This script takes data from Vizier to apply a correction in the parallax to data from GaiaDR3
'''

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from zero_point import zpt

zpt.load_tables()
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

#Obtaining the values of the correction in Parallax
zpvals = zpt.get_zpt(gmag, nueffused, psc, ecl_lat, soltype)

#This is made so the func get_zpt() don't fail (sources with 2-p solutions)
valid = soltype>3
zpvals = zpt.get_zpt(gmag[valid], nueffused[valid], psc[valid],
                     ecl_lat[valid], soltype[valid],_warnings=False)

#Applying the correction and add a new column with the corrected values
data = data[valid]
corr_plx = data['Plx'] - zpvals
df2 = df.assign(Plx_corr = corr_plx)

#Saving the new table
#The "obj" must be changed to the name of the cluster studied
t.write('obj_gaia_plxcorr.fits')


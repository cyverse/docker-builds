#!/usr/bin/env python

import fitsio
import sys
from astropy.wcs import WCS
import csv
from astropy.io import fits
import os
import numpy as np


#filepath = '/data2/fantasyzhn/Globular_cluster/Galex_image/'
filepath = sys.argv[1]
filename = []
for file in os.listdir(filepath):
    if file.endswith(".fits"):
        filename.append(os.path.join(filepath, file))
        #filename.append(file)

keys = ['NAXIS','NAXIS1','NAXIS2','CDELT1','CDELT2','CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2']
fieldsname = ['filepath','Number of axis','Dimension of axis1','Dimension of axis2','Interval of axis1','Interval of axis2','Projection type of axis1','Projection type of axis2','Central pixel of axis1','Central pixel of axis2','Value of central pixel1','Value of central pixel2']
#writer = csv.DictWriter('metadata.csv',fieldnames=fieldsname,extrasaction='ignore',delimiter = ',')
#writer.writeheader()
writer=csv.writer(open('metadata.csv','wa'))
writer.writerow(fieldsname)
#for word in fieldsname:
#    writer.writerow([word])
for i, name in enumerate(filename):
    fits_index = fitsio.FITS(name)
    hdu = fits_index[0].read_header()
    keys_value = [name]
    for j in keys:
        keys_value.append(hdu[j])
    keys_value = np.array(keys_value)	
    writer.writerow(keys_value)

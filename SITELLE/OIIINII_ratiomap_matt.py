import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

#Assign corresponding flux maps
OIII5007 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII5007correctionbin2.fits')
Hbeta = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Hbetacorrectionbin2.fits')
Halpha = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Halphacorrectionbin2.fits')
NII6583 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/NII6584correctionbin2.fits')
#Retrieving the data from fits files 
data_Hbeta = Hbeta[0].data  
data_NII6583 = NII6583[0].data  
data_Halpha = Halpha[0].data 
data_OIII5007 = OIII5007[0].data 
#Calculating O3N2 metallicity
R = np.log10((data_OIII5007/data_Hbeta) * (data_Halpha/data_NII6583))
logOH = 8.73 - 0.32 * R

with fits.open('Header.fits') as hdu_list: #Using Header fits to keep header from original fits files
    hdu = hdu_list[0]
    hdu.data = logOH
    hdu.writeto('O3N2_Metallicity.fits', overwrite=True)
#Using err flux maps
OIII5007err = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII5007correctionbin2_err.fits')
Hbetaerr = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Hbetacorrectionbin2_err.fits')
Halphaerr = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Halphacorrectionbin2_err.fits')
NII6583err = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/NII6584correctionbin2_err.fits')
#Assigning variables
A = OIII5007[0].data 
B = Hbeta[0].data 
C = Halpha[0].data 
D = NII6583[0].data 

sa = OIII5007err[0].data 
sb = Hbetaerr[0].data 
sc = Halphaerr[0].data 
sd = NII6583err[0].data 
#Calculating uncertainities
xnom = np.log10((A/B) * (C/D))
ynom = 8.73 - 0.32 * xnom

xadd = np.log10(((A+sa)/(B+sb)) * ((C+sc)/(D+sd)))
yadd = 8.73 - 0.32 * xadd

xsub = np.log10(((A-sa)/(B-sb)) * ((C-sc)/(D-sd)))
ysub = 8.73 - 0.32 * xsub

deltaadd = np.abs(yadd-ynom)
deltasub = np.abs(ysub-ynom)

highvalues = np.zeros((250, 287))

#Keeping the largest difference
for x in range(0,250):
    for y in range (0,287):
        if deltaadd[x][y] >= deltasub[x][y]:
            highvalues[x][y] = deltaadd[x][y]
        else:
            highvalues[x][y] = deltasub[x][y]
            
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = highvalues
    hdu.writeto('O3N2Uncertainities.fits', overwrite=True)

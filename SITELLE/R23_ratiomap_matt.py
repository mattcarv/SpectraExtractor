import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

#Assigning corresponding flux maps 
OII3726 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OII3726correctionbin2.fits')
OIII4959 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII4959correctionbin2.fits')
OIII5007 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII5007correctionbin2.fits')
Hbeta = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Hbetacorrectionbin2.fits')
#Retrieving the data from all the fits files 
data_Hbeta = Hbeta[0].data    
data_OIII5007 = OIII5007[0].data  
data_OIII4959 = OIII4959[0].data  
data_OII3726 = OII3726[0].data 
#Calculating Metallicity
R = np.log10((data_OII3726 + data_OIII4959 + data_OIII5007) / data_Hbeta)
R = np.log10(R)
logOH = 9.265 - (0.33*R) - (0.202*(R**2)) - (0.207*(R**3)) - (0.333*(R**4))

with fits.open('Header.fits') as hdu_list: #keeping header from previous flux map
    hdu = hdu_list[0]
    hdu.data = logOH
    hdu.writeto('R23_Metallicity_new.fits', overwrite=True)  # to write just that HDU to a new file, or
OII3726err = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OII3726correctionbin2_err.fits')
OIII4959err = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII4959correctionbin2_err.fits')
OIII5007err = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OIII5007correctionbin2_err.fits')
Hbetaerr = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Hbetacorrectionbin2_err.fits')
#retrieving data and assigning variables
A = OII3726[0].data    
B = OIII4959[0].data    
C = OIII5007[0].data    
D = Hbeta[0].data    

sa = OII3726err[0].data    
sb = OIII4959err[0].data    
sc = OIII5007err[0].data    
sd = Hbetaerr[0].data    
#Calculating uncertainities
xnom = np.log10((A + B + C) / D)
xnom = np.log10(xnom)
ynom = 9.265 - (0.33*xnom) - (0.202*(xnom**2)) - (0.207*(xnom**3)) - (0.333*(xnom**4))

xadd = np.log10(((A+sa) + (B+sb) + (C+sc)) / (D+sd))
xadd = np.log10(xadd)
yadd = 9.265 - (0.33*xadd) - (0.202*(xadd**2)) - (0.207*(xadd**3)) - (0.333*(xadd**4))

xsub = np.log10(((A-sa) + (B-sb) + (C-sc)) / (D-sd))
xsub = np.log10(xsub)
ysub = 9.265 - (0.33*xadd) - (0.202*(xadd**2)) - (0.207*(xadd**3)) - (0.333*(xadd**4))

deltaadd = np.abs(yadd-ynom)
deltasub = np.abs(ysub-ynom)

highvalues = np.zeros((250, 287))

#finding the largest difference
for x in range(0,250):
    for y in range (0,287):
        if deltaadd[x][y] >= deltasub[x][y]:
            highvalues[x][y] = deltaadd[x][y]
        else:
            highvalues[x][y] = deltasub[x][y]
            
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = highvalues
    hdu.writeto('R23Uncertainities_new.fits', overwrite=True)

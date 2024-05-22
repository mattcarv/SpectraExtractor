import numpy as np
import math
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.utils.data import get_pkg_data_filename
#Add the fits files from LUCI here for corresponding wavelength
OII3726 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OII3726correctionbin2_cardelli.fits')
NII6583 = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/NII6584correctionbin2_cardelli.fits')
#Retrieving the data from all the fits files  
data_NII6583 = NII6583[0].data  
data_OII3726 = OII3726[0].data  
NIIOII = np.log10(data_NII6583/data_OII3726)
logOH = np.log10(1.54020 + (1.26602*NIIOII) + (0.167977*(NIIOII**2))) + 8.93

with fits.open('Header.fits') as hdu_list: #using the header from the previous fits files
    hdu = hdu_list[0]
    hdu.data = logOH
    hdu.writeto('NIIOII_Metallicity_cardelli.fits', overwrite=True)  # Write the N202 numbers to a ratio map
    
    
#%%
NII6583_uncer = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/NII6584correctionbin2_err_cardelli.fits')
OII3726_uncer = fits.open('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/OII3726correctionbin2_err_cardelli.fits')
#assigning variables to make equations simpler
A = NII6583[0].data
B = OII3726[0].data
sa = NII6583_uncer[0].data
sb = OII3726_uncer[0].data
#Calculating the uncertainity to a map for uncertainities
xnom = np.log(A/B)
ynom = np.log10(1.54020 + 1.26602*xnom + 0.16799*(xnom**2))+8.93

xadd = np.log((A+sa)/(B+sb))
yadd = np.log10(1.54020 + 1.26602*xadd + 0.16799*(xadd**2))+8.93

xsub = np.log((A-sa)/(B-sb))
ysub = np.log10(1.54020 + 1.26602*xsub + 0.16799*(xsub**2))+8.93

deltaadd = np.abs(yadd-ynom)
deltasub = np.abs(ysub-ynom)

highvalues = np.zeros((250, 287))

#found the largest difference and used that for the final map
for x in range(0,250):
    for y in range (0,287):
        if deltaadd[x][y] >= deltasub[x][y]:
            highvalues[x][y] = deltaadd[x][y]
        else:
            highvalues[x][y] = deltasub[x][y]

#creating the map for uncertainities            
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = highvalues
    hdu.writeto('N202Uncertainities_cardelli.fits', overwrite=True) 

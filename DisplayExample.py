import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
import scipy.constants as scp

# Load the FITS file
# hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/UGC 2885 H-Alpha Fits Files/Rubin_High_Binning_2_Halpha_Flux.fits')
hdul = fits.open('fc_58.268000+35.591940_dss_dss2ir.fits')
header = hdul[0].header
data = hdul[0].data

IR_100 = np.sum(data)
print(IR_100, 'MJy')

# Get the WCS object from the FITS header
wcs = WCS(hdul[0].header)

# Display the image and add a color bar
fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im = ax.imshow(data, cmap='inferno')
cbar = plt.colorbar(im)
cbar.set_label('Flux')

# Display the coordinates
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
# ax.coords.grid(False, color='white', ls='solid')
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')


# subplot_kw={'projection': wcs}
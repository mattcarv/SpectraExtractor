import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np

# Load the FITS file
# hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/UGC 2885 H-Alpha Fits Files/Rubin_High_Binning_2_Halpha_Flux.fits')
hdul = fits.open('fc_58.268000+35.591940_iras_100.fits')
# header = hdul[0].header
data = hdul[0].data

#hdul2 = fits.open('C:/Users/mathe/OneDrive/Documents/PROJECTUGC2885-2022/UGC 2885 H-Alpha Fits Files/moment0co10.fits')
#hdul2 = fits.open('C:/Users/mathe/Downloads/UGC2885FILES-part1/fc_58.268000+35.591940_2mass_k.fits')
# wcs = WCS(hdul[0].header)

max_val = np.nanmax(data)
total = np.nansum(np.clip(data, np.nanmin(data), max_val))
print(total)
print(np.nansum(data[data>0]))

# # Display the image and add a color bar
# fig, ax = plt.subplots(subplot_kw={'projection': wcs})
# im = ax.imshow(StarFR, cmap='inferno', vmin=0.1e-5)
# cbar = plt.colorbar(im)
# cbar.set_label('Star Formation Rate ($M_{\odot} yr^{-1}$)')

# # Display the coordinates
# ax.set_xlabel('Right Ascension')
# ax.set_ylabel('Declination')
# # ax.coords.grid(False, color='white', ls='solid')
# ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
# ax.coords['dec'].set_axislabel('Declination (J2000)')

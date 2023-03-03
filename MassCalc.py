import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
import scipy.constants as scp

# Load the FITS file
# hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/UGC 2885 H-Alpha Fits Files/Rubin_High_Binning_2_Halpha_Flux.fits')
hdul = fits.open('lumtest.fits')
header = hdul[0].header
data = hdul[0].data

# Get the WCS object from the FITS header
wcs = WCS(hdul[0].header)


# Display the image and add a color bar
# fig, ax = plt.subplots(subplot_kw={'projection': wcs})
# im = ax.imshow(data, cmap='inferno', vmin=0)
# cbar = plt.colorbar(im)
# cbar.set_label('CO Line Luminosity (K Km/s $pc^2$)')

# # Display the coordinates
# ax.set_xlabel('Right Ascension')
# ax.set_ylabel('Declination')
# # ax.coords.grid(False, color='white', ls='solid')
# ax.coords['ra'].set_axislabel('RA (J2000)')
# ax.coords['dec'].set_axislabel('Dec (J2000)')

plt.clf()

# def LineLuminosity (data):
#     k = scp.Boltzmann
#     beam = 1.588e3
#     c = scp.c
#     nu = 115.2712018000
#     D_l = 84.105
#     z = 0.01935
    
#     luminosity = 23.5 * beam * data * (D_l**2) * ((1 + z)**(-3))
    

#     return luminosity



# data2 = LineLuminosity(data)

# hdu = fits.PrimaryHDU(data2, header)
# hdu.writeto('lumtest.fits', overwrite=True)

M = 4.3 * data

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im = ax.imshow(M, cmap='jet', vmin=0)
cbar = plt.colorbar(im)
cbar.set_label('Molecular Gas Mass ($M_{\odot}$)')

# Display the coordinates
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
# ax.coords.grid(False, color='white', ls='solid')
ax.coords['ra'].set_axislabel('RA (J2000)')
ax.coords['dec'].set_axislabel('Dec (J2000)')
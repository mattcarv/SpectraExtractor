import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
import scipy.constants as scp

# Load the FITS file
# hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/UGC 2885 H-Alpha Fits Files/Rubin_High_Binning_2_Halpha_Flux.fits')
hdul = fits.open('moment0final.fits')
header = hdul[0].header
data = hdul[0].data

# Get the WCS object from the FITS header
wcs = WCS(hdul[0].header)
print(wcs)

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

# plt.clf()

def LineLuminosity (data):
    
    '''
    This function will calculate the line luminosity based on the line density
    flux (or moment-0 map).
    
    beam: Solid angle of the observation in arcsec squared
    D_l: Luminosity distance of the source in Mpc
    z: redshift of the source
    data: density flux map in K Km/s
        
    '''
    
    beam = 1.588e3
    c = scp.c
    nu = 115.2712018000
    D_l = 84.105
    z = 0.01935
    
    luminosity = 23.5 * beam * data * (D_l**2) * ((1 + z)**(-3))
    

    return luminosity

# the line luminosity will be in units of: K Km/s pc**2

lum = LineLuminosity(data)

hdu = fits.PrimaryHDU(lum, header)
hdu.writeto('lumtest.fits', overwrite=True)

hdul = fits.open('lumtest.fits')
header = hdul[0].header
data = hdul[0].data

wcs = WCS(hdul[0].header)

M = 4.3 * data

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im = ax.imshow(M, cmap='jet', vmin=0)
cbar = plt.colorbar(im)
cbar.set_label('Molecular Gas Mass ($M_{\odot}$)')

# Display the coordinates
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
# ax.coords.grid(False, color='white', ls='solid')
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

print('Max of the molecular gas mass: ', np.nanmax(M), 'solar masses')
print('Sum of the molecular gas mass: ', np.nansum(M), 'solar masses')
print('Mean of the molecular gas mass: ', np.nanmean(M), 'solar masses')
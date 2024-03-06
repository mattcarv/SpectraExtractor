import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
import numpy as np


hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/montage/UGC2885FILES-part1/projectedwise1.fits')
header = hdul[0].header
data = hdul[0].data
wcs = WCS(hdul[0].header)

hdul2 = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/montage/UGC2885FILES-part1/projectedwise4.fits')
data2 = hdul2[0].data
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10,8)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs}, 
                               figsize=(10, 8))

im1 = plt.imshow(np.log10(data), origin='lower', cmap='Greys')
im2 = plt.contour(np.log10(data2), origin='lower', cmap='PuRd', alpha=0.9, vmin=2.3218)
plt.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
plt.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
         ha='center', va='bottom', fontsize=16)
plt.grid(color='black', lw=0.5, alpha=0.5)
plt.text(24, 58, 'Flux (log DN)')

cax1 = fig.add_axes([0.165, 0.95, 0.77, 0.04])
cbar1 = fig.colorbar(im1, cax=cax1, orientation='horizontal', spacing='proportional')
cbar1.ax.tick_params(direction='out', labeltop=True, labelbottom=False, top=True, bottom=False)

cax2 = fig.add_axes([0.96, 0.1, 0.03, 0.8])
cbar2 = fig.colorbar(im2, cax=cax2, pad=0.03)
cbar2.set_label('Flux (log DN)')


ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.subplots_adjust(left=0.15,
                    bottom=0.1,
                    right=0.95,
                    top=0.9,
                    wspace=0.6,
                    hspace=0)
plt.show()


#%%

# Load FITS data
hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/CO-files-20221207T192945Z-001/moment0final.fits')
header = hdul[0].header
data = hdul[0].data
wcs = WCS(header)

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10, 8)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs}, figsize=(10, 8))

# Plot the moment map
im1 = plt.imshow(data, origin='lower', cmap='winter')

circ_radius = 13.3*u.arcsec
conv_circ = circ_radius.to(u.deg)
circ_pixel = (conv_circ/pixel_scale).value

circle = plt.Circle((58, 8), circ_pixel, color='red', fill=False, 
                    lw=1, ls='--')
ax.add_patch(circle)
plt.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
plt.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black', ha='center', va='bottom', fontsize=16)
plt.grid(color='black', lw=0.5, alpha=0.5)

cbar1 = fig.colorbar(im1, spacing='proportional')
cbar1.set_label('Integrated Intensity (K Km s$^{-1}$)')

ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.show()

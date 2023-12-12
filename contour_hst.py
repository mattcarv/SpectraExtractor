import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from PIL import Image
from scipy import ndimage
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,12)

hdul = fits.open('/home/mdocarm/Downloads/iramco10reproject.fits')
data_iram = hdul[0].data
wcs = WCS(hdul[0].header)

hdul2 = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')
data2 = hdul2[0].data

hdul3 = fits.open('/home/mdocarm/Downloads/hideidrereproject.fits')
data3 = hdul3[0].data
# rotated_hst = ndimage.rotate(data2, 44.97, reshape=False)
# rotated_radio = ndimage.rotate(data_iram, 44.97, reshape=False)

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im_png = ax.imshow(data2, origin='lower', cmap='BuPu', vmin=0.0001, vmax=0.05)
#im_fits = ax.contour(data_iram, origin='lower', cmap='winter')
im_hi = ax.contour(data3, cmap='Greens', vmin=0, vmax=500)

ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

cbar = plt.colorbar(im_hi)
cbar.set_label('Integrated Flux ($Jy\; beam^{-1}\; m \; s^{-1}$)')

plt.show()

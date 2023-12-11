from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename


hdu1 = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')[0]
hdu2 = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/U2885-HI/projected_deidre.fits')[0]

from astropy.wcs import WCS
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [18, 8]

ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(hdu1.data, origin='lower', vmax=0.12, vmin=0.001)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')


ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(hdu2.data, origin='lower', vmax=0.12, vmin=0.001)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax1.set_title('f475')
ax2.set_title('f814')
#%%

from reproject import reproject_interp
array, footprint = reproject_interp(hdu2, hdu1.header)

#%%

from astropy.wcs import WCS
import matplotlib.pyplot as plt

ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(array, origin='lower', vmax=0.12, vmin=0.001)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Reprojected MSX band E image')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu1.header))
ax2.imshow(footprint, origin='lower')
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.coords['dec'].set_axislabel_position('r')
ax2.coords['dec'].set_ticklabel_position('r')
ax2.set_title('MSX band E image footprint')

fits.writeto('/home/mdocarm/Downloads/hideidrereproject.fits', array, hdu1.header, overwrite=True)
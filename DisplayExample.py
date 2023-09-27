import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u

hdul2 = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/CO-files-20221207T192945Z-001/moment2final.fits')
data2 = hdul2[0].data
wcs = WCS(hdul2[0].header)
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10,8)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im1 = plt.imshow(data2, origin='lower', cmap='coolwarm')


ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
          ha='center', va='bottom', fontsize=16)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
cbar = plt.colorbar(im1)
cbar.set_label('Line Width $km \; s^{-1}$')
plt.show()

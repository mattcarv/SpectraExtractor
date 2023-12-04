import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from scipy import ndimage

hdul2 = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/montage/UGC2885FILES-part1/f606_rotated.fits')
data2 = hdul2[0].data
wcs = WCS(hdul2[0].header)
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,12)

# data2 = ndimage.rotate(data2, 44.97, reshape=False)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im1 = plt.imshow(data2, origin='lower', cmap='Greys_r', vmin=0.0001, vmax=0.125)

ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
          ha='center', va='bottom', fontsize=16)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

# hdu = fits.PrimaryHDU(data2, header=hdul2[0].header)
# hdu.writeto('f606_rotated.fits')

# cbar = plt.colorbar(im1)
# cbar.set_label('Flux')
# plt.grid(visible=True, axis='both')
plt.show()

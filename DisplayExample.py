import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from PIL import Image
from scipy import ndimage

hdul2 = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')
image = Image.open('HST_rgb.png')
# data2 = hdul2[0].data
wcs = WCS(hdul2[0].header)
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,12)

data2 = ndimage.rotate(image, -44.97, reshape=False)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im1 = plt.imshow(data2, origin='lower', cmap='Greys_r', vmin=0.0001, vmax=0.125)

plt.plot([250, 250 + arcmin_pixel], [320, 320], color='white', lw=2)
plt.text(250 + arcmin_pixel / 2, 290, '1 arcmin', color='white',
          ha='center', va='bottom', fontsize=16)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
# ax.invert_xaxis()
ax.invert_yaxis()
ax.set_facecolor('k')

# Pixel coordinates for the bottom right of the image
pixel_x = 7984 - 500  # Adjust the value as needed
pixel_y = 8695 - 500  # Adjust the value as needed

# Add legend with arrows for North and East at the bottom right
arrow_length = 1000
arrow_head_width = 30

# Arrow for North
ax.annotate("N", xy=(pixel_x, pixel_y), xytext=(pixel_x, pixel_y),
            arrowprops=dict(arrowstyle="->", lw=2, color='white'),
            ha='center', va='center', color='white', fontsize=14)

# Arrow for East
ax.annotate("E", xy=(pixel_x, pixel_y), xytext=(pixel_x + arrow_length, pixel_y),
            arrowprops=dict(arrowstyle="->", lw=2, color='white'),
            ha='center', va='center', color='white', fontsize=14)

# hdu = fits.PrimaryHDU(data2, header=hdul2[0].header)
# hdu.writeto('f606_rotated.fits')

# cbar = plt.colorbar(im1)
# cbar.set_label('Flux')
# plt.grid(visible=True, axis='both')

plt.show()

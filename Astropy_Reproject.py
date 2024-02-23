from astropy.io import fits

# Here I reproject all HST files to keep their coordinates consistent.

# First I reproject all HST files to the f475 coordinates (for some reason they
# had different number for axis pixels)

#Then, I reproject all 3 to the sky coordinates I use with all the other files 
# (hdu2)
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

from reproject import reproject_interp
array, footprint = reproject_interp(hdu2, hdu1.header)

#%% Showing the reprojection

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

#%% Creating the RGB composite

from astropy.visualization import make_lupton_rgb

# Read in the three images downloaded from here:
g_data = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')[0].data
r = fits.open('/home/mdocarm/Downloads/f606_reprojected_final.fits')[0].data
i = fits.open('/home/mdocarm/Downloads/f814_reprojected_final.fits')[0].data

rgb_default = make_lupton_rgb(i, r, g_data, Q=5, stretch=0.1,
                              filename="HST_rgb.png")
plt.imshow(rgb_default, origin='lower')

#%% Displaying the image.

from astropy import units as u
from PIL import Image
from scipy import ndimage
from astropy.wcs.utils import proj_plane_pixel_scales

hdul2 = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')
image = Image.open('/home/mdocarm/Downloads/HST_rgb.png')
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

# Coordinates for the arrow
pixel_x = 7984 - 1200  
pixel_y = 8695 - 500  
arrow_length = 850
arrow_head_width = 30

# Arrow for North
ax.annotate("", xy=(pixel_x, pixel_y-arrow_length), xytext=(pixel_x, pixel_y),
            arrowprops=dict(arrowstyle="->", lw=2, color='white'),
            ha='center', va='center', color='white', fontsize=14)

# Arrow for East
ax.annotate("", xy=(pixel_x-arrow_length, pixel_y), xytext=(pixel_x, pixel_y),
            arrowprops=dict(arrowstyle="->", lw=2, color='white'),
            ha='center', va='center', color='white', fontsize=14)

plt.text(6784, 7345, 'N', color='w', ha='center', va='bottom', fontsize=14)
plt.text(5900, 8295, 'E', color='w', ha='center', va='bottom', fontsize=14)

# hdu = fits.PrimaryHDU(data2, header=hdul2[0].header)
# hdu.writeto('f606_rotated.fits')

# cbar = plt.colorbar(im1)
# cbar.set_label('Flux')
# plt.grid(visible=True, axis='both')

plt.show()
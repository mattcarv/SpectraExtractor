import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u


hdul = fits.open('projectedwise1.fits')
header = hdul[0].header
data = hdul[0].data
wcs = WCS(hdul[0].header)

hdul2 = fits.open('projectedwise4.fits')
data2 = hdul2[0].data
plt.rcParams.update({'font.size': 18})

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value


fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': wcs}, 
                               figsize=(18, 8))


im1 = ax1.imshow(data, origin='lower', cmap='Greys')
# ax1.set_title('WISE band 1 ($3.6\;microns$)')
ax1.text(6, 40, '$3.6\;\mu m$', color='black', fontsize=16)
ax1.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax1.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
         ha='center', va='bottom', fontsize=16)
ax1.grid(color='black', lw=0.5, alpha=0.5)
cax1 = fig.add_axes([0.47, 0.22, 0.02, 0.56])
cbar1 = fig.colorbar(im1, cax=cax1, ticks=[200, 450, 700, 950, 1200])
cbar1.set_label('Flux (DN)')


im2 = ax2.imshow(data2, origin='lower', vmax=212, cmap='Greys')
# ax2.set_title('WISE band 4 ($22\;microns$)')
ax2.text(6, 40, '$22\;\mu m$', color='black', fontsize=16)
ax2.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax2.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
         ha='center', va='bottom', fontsize=16)
ax2.grid(color='black', lw=0.5, alpha=0.5)
cax2 = fig.add_axes([0.96, 0.22, 0.02, 0.56])
cbar2 = fig.colorbar(im2, cax=cax2, pad=0.01, 
                     ticks=[2.1e2, 2.105e2, 2.11e2, 2.115e2, 2.12e2])
cbar2.set_label('Flux (DN)')

ax1.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax1.coords['dec'].set_axislabel('Declination (J2000)')
ax2.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax2.coords['dec'].set_axislabel('Declination (J2000)')

plt.subplots_adjust(left=0.15,
                    bottom=0.1,
                    right=0.95,
                    top=0.9,
                    wspace=0.6,
                    hspace=0)
plt.show()

import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u


hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/U2885-HI/projected_deidre.fits')
header = hdul[0].header
data = hdul[0].data/1000
wcs = WCS(hdul[0].header)

hdul2 = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/CO-files-20221207T192945Z-001/moment0finalreal.fits')
data2 = hdul2[0].data
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10,8)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value

fig, ax = plt.subplots(subplot_kw={'projection': wcs}, 
                               figsize=(10, 8))

im1 = plt.imshow(data2, origin='lower', cmap='winter')
im2 = plt.contour(data, origin='lower', cmap='PuRd', alpha=0.9, vmin=0)
plt.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
plt.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
         ha='center', va='bottom', fontsize=16)
plt.grid(color='black', lw=0.5, alpha=0.5)
plt.text(18, 80, 'Integrated Flux ($K\; km\; s^{-1}$)')

cax1 = fig.add_axes([0.23, 0.95, 0.64, 0.04])
cbar1 = fig.colorbar(im1, cax=cax1, orientation='horizontal')
cbar1.ax.tick_params(direction='out', labeltop=True, labelbottom=False, top=True, bottom=False)

cax2 = fig.add_axes([0.9, 0.1, 0.025, 0.8])
cbar2 = fig.colorbar(im2, cax=cax2, pad=0.01)
cbar2.set_label('Integrated Flux ($Jy\; beam^{-1} \; m\; s^{-1}$)')


ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.subplots_adjust(left=0.15,
                    bottom=0.1,
                    right=0.95,
                    top=0.9,
                    wspace=0.6,
                    hspace=0)
plt.show()

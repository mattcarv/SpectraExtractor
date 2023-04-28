import matplotlib.pyplot as plt
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

with fits.open('Rubin_High_Binning_2_Halpha_Flux.fits') as hdul:
    data = hdul[0].data
    wcs = WCS(hdul[0].header)


pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel

arcmin = 1 * u.arcmin.to(u.deg)

arcmin_pixel = (arcmin / pixel_scale).value

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=wcs)

ax.imshow(data, origin='lower', cmap='gray', vmax=9e-18, vmin=0)
ax.plot([10, 10 + arcmin_pixel], [10, 10], color='white', lw=2)
ax.text(10 + arcmin_pixel / 2, 12, '1 arcmin', color='white',
         ha='center', va='bottom', fontsize=12)


ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.show()

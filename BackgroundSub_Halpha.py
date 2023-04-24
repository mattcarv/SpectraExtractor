import dask.array as da
import astropy.io.fits as fits
from scipy.ndimage.filters import median_filter
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS
from astropy import units as u

# Load FITS image as Dask array
with fits.open('Rubin_High_Binning_2_Halpha_Flux.fits') as hdul:
    image = da.from_array(hdul[0].data)
    wcs = WCS(hdul[0].header)

# Smooth image using median filter
smoothed = da.from_array(median_filter(image.compute(), size=3))

# Compute median background level using Dask
median = da.from_array(np.nanpercentile(smoothed.compute(), 80))

# Subtract background from image
background_subtracted = image - median

new_array = np.array(background_subtracted)
max_val = np.nanmax(new_array)
total = np.nansum(np.clip(new_array, np.nanmin(new_array), max_val))

print(total)
# max_val = da.nansum(background_subtracted).compute()
# print(max_val)

# Save the subtracted file
header = hdul[0].header
hdu = fits.PrimaryHDU(background_subtracted.compute())
hdu.header = header

# hdu.writeto("output.fits", overwrite=True)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel

arcmin = 1 * u.arcmin.to(u.deg)

arcmin_pixel = (arcmin / pixel_scale).value

# Plot both figures
fig, ax = plt.subplots(1, 2, figsize=(15, 6))

im1 = ax[0].imshow(image, cmap='viridis', vmin=0, vmax=1e-17)
ax[0].invert_yaxis()
ax[0].set_title('Original')

ax[0].plot([10, 10 + arcmin_pixel], [10, 10], color='white', lw=2)
ax[0].text(10 + arcmin_pixel / 2, 12, '1 arcmin', color='white',
         ha='center', va='bottom', fontsize=12)

# Plot second array in second subplot
im2 = ax[1].imshow(background_subtracted, cmap='viridis', vmin=0, vmax=1e-17)
ax[1].invert_yaxis()
ax[1].set_title('Subtracted using a median filter')

ax[1].plot([10, 10 + arcmin_pixel], [10, 10], color='white', lw=2)
ax[1].text(10 + arcmin_pixel / 2, 12, '1 arcmin', color='white',
         ha='center', va='bottom', fontsize=12)

# Add colorbars
cbar1 = fig.colorbar(im1, ax=ax[0])
cbar1.set_label('Flux ($erg/s/cm^{2}/A$)')
cbar2 = fig.colorbar(im2, ax=ax[1])
cbar2.set_label('Flux ($erg/s/cm^{2}/A$)')

# Show the figure
plt.show()

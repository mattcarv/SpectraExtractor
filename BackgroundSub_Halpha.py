import dask.array as da
import astropy.io.fits as fits
from scipy.ndimage.filters import median_filter
import matplotlib.pyplot as plt
import numpy as np

# Load FITS image as Dask array
with fits.open('Rubin_High_Binning_2_Halpha_Flux.fits') as hdul:
    image = da.from_array(hdul[0].data)

# Smooth image using median filter
smoothed = da.from_array(median_filter(image.compute(), size=3))

# Compute median background level using Dask
median = da.from_array(np.nanpercentile(smoothed.compute(), 75))

# Subtract background from image
background_subtracted = image - median

new_array = np.array(background_subtracted)
max_val = np.nanmax(new_array)
total = np.nansum(np.clip(new_array, 0, max_val))
print(total)
# max_val = da.nansum(background_subtracted).compute()
# print(max_val)

# Save the subtracted file
header = hdul[0].header
hdu = fits.PrimaryHDU(background_subtracted.compute())
hdu.header = header

hdu.writeto("output.fits", overwrite=True)

# Plot both figures
fig, ax = plt.subplots(1, 2, figsize=(15, 6))

im1 = ax[0].imshow(image, cmap='viridis', vmin=0, vmax=1e-17)
ax[0].invert_yaxis()
ax[0].set_title('Original')

# Plot second array in second subplot
im2 = ax[1].imshow(background_subtracted, cmap='viridis', vmin=0, vmax=1e-17)
ax[1].invert_yaxis()
ax[1].set_title('Subtracted using a median filter')

# Add colorbars
cbar1 = fig.colorbar(im1, ax=ax[0])
cbar1.set_label('Flux ($erg/s/cm^{2}/A$)')
cbar2 = fig.colorbar(im2, ax=ax[1])
cbar2.set_label('Flux ($erg/s/cm^{2}/A$)')

# Show the figure
plt.show()
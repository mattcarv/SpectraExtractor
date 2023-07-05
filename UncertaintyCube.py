import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
# Load the line width map
width_hdu = fits.open('moment2final.fits')[0]
width_map = width_hdu.data
header = width_hdu.header
wcs = WCS(width_hdu.header)

uncertainty_integrated = 0.085 * (width_map/(np.sqrt(width_map)))


fig, ax = plt.subplots(subplot_kw={'projection': wcs})
plt.imshow(uncertainty_integrated, origin='lower', cmap='Greys')

uncertainty_integrated = uncertainty_integrated[~np.isnan(uncertainty_integrated)]
squared_values = np.square(uncertainty_integrated)
sum_of_squares = np.sum(squared_values)
result = np.sqrt(sum_of_squares)

plt.text(0.05, 0.9, f'Total: {result:.5f} $K \; Km/s$', transform=plt.gca().transAxes)
plt.colorbar(label='Uncertainty (K km/s)')
plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')
plt.title('Uncertainty in Integrated Intensity')
plt.show()
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np

# Define constants
beam = 558
D_l = 84.11
z = 0.01935
conversion_factor = 4.3 * 23.5 * beam * (D_l**2) * ((1 + z)**(-3))

# Load the FITS file
hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/CO-files-20221207T192945Z-001/CO-files/moment0final.fits')
data = hdul[0].data

# Calculate line luminosity
luminosity = conversion_factor * data

# Calculate the total luminosity
total_luminosity = np.nansum(luminosity)

# Define coordinates for top and bottom features
top_feature_x = np.arange(33, 42)
top_feature_y = np.arange(32, 41)
bot_feature_x = np.arange(23, 32)
bot_feature_y = np.arange(24, 33)

# Create WCS
wcs = WCS(hdul[0].header)

# Create a figure
fig, ax = plt.subplots(subplot_kw={'projection': wcs})

# Plot the line luminosity
im = ax.imshow(luminosity, cmap='gray', vmin=0)
plt.scatter(top_feature_x, top_feature_y, c='red', marker='x')
plt.scatter(bot_feature_x, bot_feature_y, c='blue', marker='x')
cbar = plt.colorbar(im)
cbar.set_label('$H_{2}$ Mass$\;$($M_{\odot}$)')

# Display the coordinates
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

# Clear the figure
plt.clf()

# Load the uncertainty map
hdul = fits.open('/home/mdocarm/Downloads/PROJECTUGC2885-2022/CO-files-20221207T192945Z-001/unc_map.fits')
data = hdul[0].data

# Calculate uncertainty
uncertainty = conversion_factor * data

# Extract top and bottom regions from the line luminosity and uncertainty
top = np.nan_to_num(luminosity[top_feature_x, top_feature_y])
bot = np.nan_to_num(luminosity[bot_feature_y, bot_feature_x])
top_unc = np.nan_to_num(uncertainty[top_feature_x, top_feature_y])
bot_unc = np.nan_to_num(uncertainty[bot_feature_y, bot_feature_x])

# Calculate mass and mass uncertainty
mass = (top + bot) / (2 * 4450827.06)
mass_unc = (top_unc + bot_unc) / (2 * 4450827.06)

# Create distances array
distances = np.arange(len(mass)) * 2.96

# Create a scatter plot with error bars
plt.scatter(distances, mass, marker='X', label='Data')
plt.errorbar(distances, mass, yerr=mass_unc, fmt='none', color='green', 
             alpha=0.5, capsize=3, label='Error')

# Plot the mass data
plt.plot(distances, mass, color='green', linestyle='-', alpha=0.5)
plt.ylabel('$\Sigma_{mol}$ ($M_{\odot}\; pc^{-2}$)')
plt.xlabel('Distance from centre ($kpc$)')

plt.show()

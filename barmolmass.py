import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np

hdul = fits.open('lumtest.fits')
header = hdul[0].header
data = hdul[0].data

wcs = WCS(hdul[0].header)

#print('Integrated Luminosity: ', np.nansum(data), 'K Km/s $pc^2$')

M = 4.3 * data


top_feature_x = np.arange(33, 42)
top_feature_y = np.arange(32, 41)

bot_feature_x = np.arange(23, 32)
bot_feature_y = np.arange(24, 33)

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im = ax.imshow(M, cmap='gray', vmin=0)
plt.scatter(top_feature_x, top_feature_y, c='red', marker='x')
plt.scatter(bot_feature_x, bot_feature_y, c='blue', marker='x')
cbar = plt.colorbar(im)
cbar.set_label('$H_{2}$ Mass$\;$($M_{\odot}$)')

# Display the coordinates
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
# ax.coords.grid(False, color='white', ls='solid')
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
plt.clf()


top = np.nan_to_num(M[top_feature_x, top_feature_y])
bot = np.nan_to_num(M[bot_feature_y, bot_feature_x])

mass = (top+bot)/(2*4450827.06)
distances = np.arange(len(mass)) * 2.96

plt.scatter(distances, mass, marker='X')

plt.errorbar(distances, mass, yerr=0.4*mass, fmt='none', color='green', 
             alpha=0.5, capsize=3, label='Error')

plt.plot(distances, mass, color='green', linestyle='-', alpha=0.5)
plt.ylabel('$\Sigma_{mol}$ ($M_{\odot}\; pc^{-2}$)')
plt.xlabel('Distance from centre ($kpc$)')
plt.show()
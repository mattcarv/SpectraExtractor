import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

hdul =  fits.open('C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_NII6548_Flux.fits')
data = hdul[0].data
wcs = WCS(hdul[0].header)

fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im1 = plt.imshow(data, origin='lower', cmap='viridis', vmax=(.001*data.max()))

#plt.plot([250, 250 + arcmin_pixel], [320, 320], color='white', lw=2)
#plt.text(250 + arcmin_pixel / 2, 290, '1 arcmin', color='white',
#          ha='center', va='bottom', fontsize=16)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.show()

#%%------------------------------------------------------------------------------

# Plot a histogram of the pixel values
plt.hist(data.flatten(), bins=100, log=True)
plt.title('Pixel Value Histogram')
plt.xlabel('Pixel Value')
plt.ylabel('Frequency (log scale)')
plt.show()
#%%

masked_data = data<0.0001e-13

# Plot the image after masking outliers
fig, ax = plt.subplots(subplot_kw={'projection': wcs})
im2 = plt.imshow(masked_data, origin='lower', cmap='viridis')

# Set axis labels
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

# Show the image after masking outliers
plt.show()

#%%
# Define a central ellipse to compute the flux.

#%% BPT Boundaries
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

# Define the function
# (0.61 / (x - 0.47)) + 1.19
def BPT(x):
    return (0.61 / (x - 0.05)) + 1.3

# Generate x values
x_values = np.linspace(-1.6, 0, 100)

# Calculate corresponding y values
y_values = BPT(x_values)

# Plot the function
plt.plot(x_values, y_values, c='r')
plt.text(-0.75, -0.5, 'STAR-FORMING', fontsize=15)
plt.text(0.2, 1, 'AGN', fontsize=15)
plt.xlim(-1.5, 0.5)
plt.ylim(-1.2, 1.5)
plt.xlabel(r'log ($NII\lambda 6584/ H\alpha$)')
plt.ylabel(r'log ($OIII\lambda 5007/ H\beta$)')
plt.show()

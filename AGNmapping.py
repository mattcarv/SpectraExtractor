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
#%%
# BPT Boundaries

a_low = 0.034
b_low = 1.447
c_low = -0.986

a_high = 0.885
b_high = -0.792
c_high = -6.712

def bptCurveN_O(x, a, b, c):
    return a + b*x + c*(x**2)

x_values_low = np.linspace(-0.45, 0.29, 400)
x_values_high = np.linspace(-0.45, -0.12, 400)


y_low = bptCurveN_O(x_values_low, a_low, b_low, c_low)
y_high = bptCurveN_O(x_values_high, a_high, b_high, c_high)

# Plot both curves
plt.plot(x_values_low, y_low, label='Lower Limit')
plt.plot(x_values_high, y_high, label='Higher Limit')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of Lower and Higher Limits')
plt.legend()
plt.grid(True)
plt.show()

#%%

# Define the function
def calculate_y(x):
    return (0.61 / (x - 0.47)) + 1.19

# Generate x values
x_values = np.linspace(-0.45, 0.29, 100)

# Calculate corresponding y values
y_values = calculate_y(x_values)

# Plot the function
plt.plot(x_values, y_values)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()

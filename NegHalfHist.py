import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

hdul = fits.open('maskedcube.fits')
header = hdul[0].header
data = hdul[0].data
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10,8)

# Assuming you have a data cube called 'data' with shape (channels, height, width)
channel_index = 45  # Index of the desired channel

# Select the desired channel from the data cube
channel_data = data[channel_index]

# Flatten the channel data into a 1D array
pixel_values = channel_data.flatten()

# Remove non-finite values (NaN, Inf)
pixel_values = pixel_values[np.isfinite(pixel_values)]

# Plot the histogram
hist, bins, _ = plt.hist(pixel_values, bins=100, density=True, alpha=0.5)

# Fit a Gaussian distribution to the histogram
params = stats.norm.fit(pixel_values)

# Generate x-values for the fitted curve
x = np.linspace(np.min(pixel_values), np.max(pixel_values), 1000)

# Evaluate the fitted Gaussian distribution
pdf = stats.norm.pdf(x, *params)

mean = params[0]
std = params[1]
print(std)

# Plot the histogram and fitted curve
plt.plot(x, pdf, 'r-', label='Gaussian Fit; $\sigma_{CO}$=0.0012 K')
plt.xlabel('Peak Intensity (K)')
plt.ylabel('Number of Pixels')
plt.axvline(x=mean - std, color='g', linestyle='--', label='Mean - $\sigma$')
plt.axvline(x=mean + std, color='g', linestyle='--', label='Mean + $\sigma$')
plt.legend()
plt.title(f'Channel {channel_index}')
plt.show()

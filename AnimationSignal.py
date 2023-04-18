import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from astropy.io import fits
from astropy.wcs import WCS

# Read in the data cube
hdu = fits.open('maskedcube.fits')[0]


hdu2 = fits.open('moment0co10.fits')
header = hdu2[0].header
wcs = WCS(hdu2[0].header)


data = hdu.data[27:54]
# Define the animation function
def animate(i):
    if np.all(data[i] == 0):
        return

    # Clear the previous plot
    plt.clf()
    vel = np.linspace(5333.3411648850, 6008.3411648850, 27)

    # Plot the image for the i-th channel
    im = plt.imshow(data[i], origin='lower')
    
    cbar = fig.colorbar(im)
    cbar.set_label('Brightness Temperature (K)')
    plt.clim(0, 0.05)

# Set up the figure and animation
fig, ax = plt.subplots(subplot_kw={'projection': wcs})
animation = FuncAnimation(fig, animate, frames=data.shape[0], interval=250)

# Save the animation as a GIF file
animation.save('animation.gif', writer='imagemagick')

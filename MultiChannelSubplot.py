import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Read in the data cube
hdu = fits.open('maskedcube.fits')[0]

data = hdu.data[26:54]

# Define the figure and subplot layout
rows = 9  # Number of rows in the grid
cols = 3  # Number of columns in the grid

# Create the figure and subplots
fig, axes = plt.subplots(rows, cols, figsize=(11, 22))

# Flatten the axes array for easy indexing
axes = axes.flatten()

# Loop through each frame and plot it in the corresponding subplot
for i, ax in enumerate(axes):
    if np.all(data[i] == 0):
        continue

    ax.imshow(data[i], origin='lower', vmin=0, vmax=0.04)
    ax.set_title(f'Channel {i+27}: {(i*25)+5333.3} $Km/s$')
    ax.set_xticks([])
    ax.set_yticks([])

# Hide empty subplots if there are fewer frames than subplots
for j in range(len(data), rows * cols):
    axes[j].axis('off')

# Adjust the spacing between subplots
fig.tight_layout()

# Add a colorbar to the last subplot
cax = fig.add_axes([1.02, 0.1, 0.02, 0.85])
cbar = fig.colorbar(axes[-1].images[0], cax=cax)
cbar.set_label('Brightness Temperature (K)')


# Show the plot
plt.show()

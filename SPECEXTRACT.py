import numpy as np
import astropy.io.fits as fits
from numpy import savetxt

def extract_spectra(data_cube, x, y):
    
    nz, ny, nx = data_cube.shape # Extract the cube shape  
    spectrum = data_cube[:, y, x]# Define the spectrum as the z-axis
    return spectrum # Return the brightness for each channel, aka spectrum

# Load the data
data_cube = fits.getdata('ugc2885_co21_v3.fits')

# num_spectra = #


# Positions in which we are extracting the specific spectrum from
positions = [(i, j) for i in range(127) for j in range(127)]
# positions = [(30, 39)]

# Loop for each of the positions
for i, pos in enumerate(positions):
    x, y = pos
    spectrum = extract_spectra(data_cube, x, y)

    savetxt(f'spectrum{x}_{y}.txt', spectrum, delimiter=',')
    
    
#------------------------------------------------------------------------------    
# from numpy import load
# data = load('spectrum_0.npy')
# print(data)
# -----------------------------------------------------------------------



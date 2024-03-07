import warnings
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from spectral_cube import SpectralCube
plt.rcParams.update({'font.size': 18})
# RUNTIME: since it is a fairly simple code (and not a big file - meaning not 
# many of channels or pixels) it runs in approx. 5 minutes on my computer.


#%% In this block, we display the original cube and derive its noise level
#   from the standard deviation in the spectrum
warnings.filterwarnings('ignore')

filename = fits.open('ugc2885_co10_v4.fits')

cube_original = SpectralCube.read(filename, format='fits', use_dask=True)

num_slices = cube_original.shape[0]

modified_slices = []

# Correct for the beam efficiency of IRAM
for i in range(num_slices):
    slice_unmasked = cube_original.unmasked_data[i, :, :] * 1.205
    modified_slices.append(slice_unmasked)

modified_cube = SpectralCube(data=np.array(modified_slices), wcs=cube_original.wcs, 
                             mask=None, meta=None)

#modified_cube.write('modifiedcube.fits', overwrite=True)

filename = fits.open('modifiedcube.fits')

cube = SpectralCube.read(filename, format='fits', use_dask=True)

# Here I get rid of the very outstanding channels that have high noise levels
cube = cube.spectral_slab(cube.spectral_axis[13], cube.spectral_axis[100])

cube = cube.with_spectral_unit(u.km / u.s)

header = filename[0].header
wcs = WCS(filename[0].header)
peak_intensity = cube.max(axis=0)

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value
circ_radius = 13.3*u.arcsec
conv_circ = circ_radius.to(u.deg)
circ_pixel = (conv_circ/pixel_scale).value
circle = plt.Circle((62, 64), circ_pixel, color='red', fill=False, 
                    lw=1, ls='--')


plt.figure(figsize=(10,8))
ax = plt.subplot(projection=peak_intensity.wcs)
ax.add_patch(circle)
ax.grid(True, color='k', lw=0.5, alpha=0.3)
im = ax.imshow(peak_intensity.value, origin='lower', cmap='viridis', vmax=0.12)
ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
          ha='center', va='bottom', fontsize=16)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
cbar = plt.colorbar(im)
cbar.set_label('Peak (K)')
plt.show()

mad_std_spectrum = cube.mad_std(axis=(1, 2))

# Mean standard deviation throughout the whole cube:


plt.figure(figsize=(10,8))
plt.plot(mad_std_spectrum.spectral_axis.value, mad_std_spectrum.value, drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r' Noise standard deviation $\sigma$ (K)')

# # Best to extend the range to 0.
plt.ylim([0.007, 0.015])
plt.axhline(0.0096, linestyle='--', color='k', linewidth=3, 
            label='Average level of $\sigma$')
plt.legend(frameon=True)
plt.show()

# Clip values above 2-sigma
cube_sclip = cube.sigma_clip_spectrally(2) 

mad_std_spectrum_sclip = cube_sclip.mad_std(axis=(1, 2))

plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, drawstyle='steps-mid')
plt.xlabel('Velocity (km/s)')
plt.ylabel(r' Noise standard deviation $\sigma$ (K)')
print(np.mean(mad_std_spectrum_sclip))
# Best to extend the range to 0.
plt.ylim([0.007, 0.010])

plt.axhline(0.0081, linestyle='--', color='k', linewidth=3, label='Average noise level')
plt.legend(frameon=True)
plt.show()
#%% Here we display the noise map - notice the edge effects from the observation

mad_std_map_sclip = cube_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension

# mad_std_map_sclip.quicklook()

plt.figure(figsize=(10,8))
ax = plt.subplot(projection=mad_std_map_sclip.wcs)
ax.grid(True, color='k', lw=0.5, alpha=0.3)
im = ax.imshow(mad_std_map_sclip.value, origin='lower', cmap='viridis', vmax=0.12)
ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
          ha='center', va='bottom', fontsize=16)
circ_radius = 13.3*u.arcsec
conv_circ = circ_radius.to(u.deg)
circ_pixel = (conv_circ/pixel_scale).value
circle = plt.Circle((62, 64), circ_pixel, color='red', fill=False, 
                    lw=1, ls='--')
ax.add_patch(circle)
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
cbar = plt.colorbar(im)
cbar.set_label('Peak (K)')
plt.show()
#%% Here we define the masks
low_snr_mask = (cube > 2 * mad_std_map_sclip).include()
high_snr_mask = (cube > 5 * mad_std_map_sclip).include()

# I decided to use 2 and 5 for the thresholds because more strict values would
# clean out the whole cube, not leaving any signal and lower values would include
# too much noise (from edge effects).

from dask_image import ndmeasure

# Find connected structures
structure = np.ones((3, 3, 3), dtype=bool)

low_snr_mask_labels, num_labels = ndmeasure.label(low_snr_mask,
                                                  structure=structure)

# Ask dask to execute the operation
num_labels = num_labels.compute()

print(f"Initial number of regions found: {num_labels}")

# From the labels, count the number of pixels within each label.

# Count how many >6 sigma pixels (high_snr_mask) are within regions defined in low_snr_mask_labels
num_pixels_in_high_snr_mask = ndmeasure.sum_labels(high_snr_mask,
                                                    label_image=low_snr_mask_labels,
                                                    index=range(1, num_labels + 1)) # +1 offset for mask labels

# Count how many >3 sigma pixels (low_snr_mask) are within regions defined in low_snr_mask_labels.
num_pixels_in_low_snr_mask = ndmeasure.sum_labels(low_snr_mask,
                                                  label_image=low_snr_mask_labels,
                                                  index=range(1, num_labels + 1)) # +1 offset for mask labels


# To preserve the low_snr_mask, we will create a new signal mask where we will remove 
# regions that do not pass the criteria.
signal_mask = low_snr_mask

low_min_pixels = 20
high_min_pixels = 3
# Above is another variable for this code, the number of pixel selected for each
# defined region. Higher resolution allows for more signal pixels - in our case
# 20 and 3 did the job for the low and high masks.

# From here on, it is pretty much the same procedure as for the signal-masking notebook


for num, (high_pix_num, low_pix_num) in enumerate(zip(num_pixels_in_high_snr_mask, num_pixels_in_low_snr_mask)):
    if high_pix_num >= high_min_pixels and low_pix_num >= low_min_pixels:
        # This region passes the criteria. Keep it in the mask.
        continue

    # Remove regions that do not pass the criteria.
    # NOTE: enumerate will start with 0, but the mask labels start at 1
    # We apply a +1 offset to `num` to account for this.
    signal_mask[low_snr_mask_labels == num + 1] = False
    
signal_mask_labels, num_labels = ndmeasure.label(signal_mask,
                                                  structure=structure)

num_labels = num_labels.compute()

print(f"Final number of regions found: {num_labels}")

from dask_image import ndmorph
from dask import array as da

structure = np.ones((3, 3), dtype=bool)

structure_spec = np.zeros((3, 3), dtype=bool)
structure_spec[1, 1] = True

# Add 1 spectral element on each side of the spatial structure.
# np.dstack stacks the arrays along a new 3rd dimension:
structure = np.dstack([structure_spec, structure, structure_spec])

# Convert to a dask array
structure = da.from_array(structure)

signal_mask = ndmorph.binary_dilation(signal_mask, structure=structure, iterations=1)


signal_mask = signal_mask.compute()

masked_cube = cube.with_mask(signal_mask)

masked_cube.write('maskedcube.fits', overwrite=True)

# I chose to have the actual file so I could use different codes to analyze it
# instead of having to run this one every single time.

# From here on, we get the first analysis of the masked cube, with the moment 
# maps. From the moment 0 we will get the luminosity and therefore, its associated 
# mass.

peak_intensity_sigmask = masked_cube.max(axis=0)

plt.figure(figsize=(10,8))

ax = plt.subplot(projection=peak_intensity_sigmask.wcs)
ax.grid(True)
im = ax.imshow(peak_intensity_sigmask.value, origin='lower', cmap='viridis')
ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
          ha='center', va='bottom', fontsize=16)

im = ax.imshow(peak_intensity_sigmask.value, origin='lower', cmap='viridis')
cbar = plt.colorbar(im)
cbar.set_label('Peak (K)')
ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.show()

# This is an example noise-reduced spectrum of a pixel. In orange you have the
# signal above the noise.
centre_spectrum = cube[:, 33, 33]
centre_spectrum_sigmask = masked_cube[:, 33, 33]

plt.plot(centre_spectrum.spectral_axis.value,
          centre_spectrum.filled_data[:].value,
          drawstyle='steps-mid', label='Original')
plt.plot(centre_spectrum_sigmask.spectral_axis.value,
          centre_spectrum_sigmask.filled_data[:].value, drawstyle='steps-mid',
          linewidth=3, label='Signal', color='orange')

plt.legend(frameon=True)

plt.xlabel("Velocity (km/s)")
plt.ylabel('Brightness Temp. (K)')
plt.show()

# %%MOMENT MAPS -----------------------------------------------------------------
masked_moment0 = masked_cube.moment0()
masked_moment1 = masked_cube.moment1()
masked_moment2 = masked_cube.moment2()
masked_linewidth = masked_cube.linewidth_sigma()

# Getting the files for each of the moment maps. Mainly to understand some radial
# properties of the velocity and dispersion.


#masked_moment0.write('moment0final.fits', overwrite=True)
#masked_moment1.write('moment1final.fits', overwrite=True)
#masked_linewidth.write('moment2final.fits', overwrite=True)

hdul = fits.open('moment0final.fits')
header = hdul[0].header
data = hdul[0].data
wcs = WCS(hdul[0].header)

# Here I calculate the luminosity based on the moment-0 map. I use the equation from
# Solomon, 1997
def LineLuminosity (data):
    
    '''
    This function will calculate the line luminosity based on the line density
    flux (or moment-0 map).
    
    beam: Solid angle of the observation in arcsec squared
    D_l: Luminosity distance of the source in Mpc
    z: redshift of the source
    data: density flux map in K Km/s
        
    '''
    
    beam = 558
    D_l = 84.34
    z = 0.01935
    
    luminosity = 23.5 * beam * (D_l**2)* data * ((1 + z)**(-3))
    
    return luminosity
    
lum = LineLuminosity(data)    

print('Integrated Luminosity: ', np.nansum(lum), '$K Km/s $pc^2$')

# This is the alpha C0 conversion from Bolatto, 2013. I used the Galactic factor.
M = 4.3 * lum

print('Sum of the molecular gas mass: ', np.nansum(M), '$M_{\odot}]$')

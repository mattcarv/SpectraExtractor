import os
from pathlib import Path
from time import time
import warnings


import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('seaborn-colorblind')

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from spectral_cube import SpectralCube

#%%
warnings.filterwarnings('ignore')
filename = fits.open('ugc2885_co10_v4.fits')

cube = SpectralCube.read(filename, format='fits', use_dask=True)

cube = cube.spectral_slab(cube.spectral_axis[13], cube.spectral_axis[100])

cube = cube.with_spectral_unit(u.km / u.s)

header = filename[0].header
wcs = WCS(filename[0].header)
peak_intensity = cube.max(axis=0)
plt.rcParams.update({'font.size': 18})
# peak_intensity.quicklook()
# plt.imshow(peak_intensity[0,:,:])

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
arcmin = 1 * u.arcmin.to(u.deg)
arcmin_pixel = (arcmin / pixel_scale).value


# plt.figure(figsize=(10,8))
# ax = plt.subplot(projection=peak_intensity.wcs)
# im = ax.imshow(peak_intensity.value, origin='lower', cmap='viridis', vmax=0.25)
# ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#           ha='center', va='bottom', fontsize=16)
# ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
# ax.coords['dec'].set_axislabel('Declination (J2000)')
# cbar = plt.colorbar(im)
# cbar.set_label('Peak (K)')

mad_std_spectrum = cube.mad_std(axis=(1, 2))

# mean standard deviation throughout the whole cube
# print(np.mean(mad_std_spectrum))

# plt.figure(figsize=(10,8))
# plt.plot(mad_std_spectrum.spectral_axis.value, mad_std_spectrum.value, drawstyle='steps-mid')
# plt.xlabel('Velocity (km/s)')
# plt.ylabel(r' Noise standard deviation $\sigma$ (K)')

# # # Best to extend the range to 0.
# plt.ylim([0.0065, 0.012])
# plt.axhline(0.00799, linestyle='--', color='k', linewidth=3, 
#             label='Average level of $\sigma$')
# plt.legend(frameon=True)

cube_sclip = cube.sigma_clip_spectrally(3) # Clip values above 3-sigma

mad_std_spectrum_sclip = cube_sclip.mad_std(axis=(1, 2))


# plt.plot(mad_std_spectrum_sclip.spectral_axis.value, mad_std_spectrum_sclip.value, drawstyle='steps-mid')
# plt.xlabel('Velocity (km/s)')
# plt.ylabel(r' Noise standard deviation $\sigma$ (K)')

# # # Best to extend the range to 0.
# plt.ylim([0.006, 0.012])

# plt.axhline(0.0077, linestyle='--', color='k', linewidth=3, label='A priori noise expectation')
# plt.legend(frameon=True)
#%%
mad_std_map_sclip = cube_sclip.mad_std(axis=0) # Calculate sigma along the spectral dimension

# mad_std_map_sclip.quicklook()

# plt.figure(figsize=(10,8))
# ax = plt.subplot(projection=mad_std_map_sclip.wcs)
# im = ax.imshow(mad_std_map_sclip.value, origin='lower', cmap='viridis', vmax=0.12)
# ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#          ha='center', va='bottom', fontsize=16)
# ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
# ax.coords['dec'].set_axislabel('Declination (J2000)')
# cbar = plt.colorbar(im)
# cbar.set_label('Peak (K)')
#%%
low_snr_mask = (cube > 2 * mad_std_map_sclip).include()
high_snr_mask = (cube > 5 * mad_std_map_sclip).include()

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

# masked_cube.write('maskedcube.fits')

peak_intensity_sigmask = masked_cube.max(axis=0)

# plt.figure(figsize=(10,8))
# ax = plt.subplot(projection=peak_intensity_sigmask.wcs)
# im = ax.imshow(peak_intensity_sigmask.value, origin='lower', cmap='viridis')
# ax.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#          ha='center', va='bottom', fontsize=16)
# ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
# ax.coords['dec'].set_axislabel('Declination (J2000)')
# cbar = plt.colorbar(im)
# cbar.set_label('Peak (K)')

# ax = plt.subplot(projection=peak_intensity_sigmask.wcs)
# im = ax.imshow(peak_intensity_sigmask.value, origin='lower', cmap='viridis')
# cbar = plt.colorbar(im)
# cbar.set_label('Peak (K)')

# # ax.invert_yaxis()
# ax.set_ylabel('Declination')
# ax.set_xlabel('Right Ascension')

# peak_intensity_sigmask.quicklook()
# plt.clf()

# centre_spectrum = cube[:, 33, 33]
# centre_spectrum_sigmask = masked_cube[:, 33, 33]

# plt.plot(centre_spectrum.spectral_axis.value,
#           centre_spectrum.filled_data[:].value,
#           drawstyle='steps-mid', label='Original')
# plt.plot(centre_spectrum_sigmask.spectral_axis.value,
#           centre_spectrum_sigmask.filled_data[:].value, drawstyle='steps-mid',
#           linewidth=3, label='Masked', color='orange')

# plt.legend(frameon=True)

# plt.xlabel("Velocity (km/s)")
# plt.ylabel('Brightness Temp. (K)')


# MOMENT MAPS -----------------------------------------------------------------
masked_moment0 = masked_cube.moment0()
masked_moment1 = masked_cube.moment1()
masked_moment2 = masked_cube.moment2()
masked_linewidth = masked_cube.linewidth_sigma()

# masked_moment0.write('moment0final.fits')
# masked_moment1.write('moment1final.fits')
# masked_linewidth.write('moment2final.fits')
wcs = masked_moment0.wcs
# #_______________________________________________________________
# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, subplot_kw={'projection': wcs}, 
#                                     figsize=(34, 8))


# ax1 = plt.subplot(projection=masked_moment0.wcs)
#----
# im1 = ax1.imshow(masked_moment0.value, origin='lower', cmap='inferno')
# ax1.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax1.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#           ha='center', va='bottom', fontsize=16)
# cbar1 = plt.colorbar(im1)
# cbar1.set_label('Integrated Intensity (K km/s)')
# ax1.set_ylabel('Declination (J2000)')
# ax1.set_xlabel('Right Ascension (J2000)')

# #_______________________________________________________________

# ax2 = plt.subplot(projection=masked_moment1.wcs)
#----
# im2 = ax2.imshow(masked_moment1.value, origin='lower', cmap='coolwarm', vmin=5000, vmax=6200)
# ax2.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax2.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#           ha='center', va='bottom', fontsize=16)
# cbar2 = plt.colorbar(im2)
# cbar2.set_label('Centroid (km/s)')
# ax2.set_ylabel('Declination (J2000)')
# ax2.set_xlabel('Right Ascension (J2000)')

# #_______________________________________________________________


# ax3 = plt.subplot(projection=masked_linewidth.wcs)
#----
# im3 = ax3.imshow(masked_linewidth.value, origin='lower', cmap='coolwarm')
# ax3.plot([6, 6 + arcmin_pixel], [6, 6], color='black', lw=2)
# ax3.text(6 + arcmin_pixel / 2, 8, '1 arcmin', color='black',
#           ha='center', va='bottom', fontsize=16)
# cbar3 = plt.colorbar(im3)
# cbar3.set_label('Line Width (km/s)')
# ax3.set_ylabel('Declination (J2000)')
# ax3.set_xlabel('Right Ascension (J2000)')

# plt.show()
#%%
# CUBE FITTING ----------------------------------------------------------------
# y, x = 33, 33
# size = 33

# moment0_cutout = masked_moment0[y-size:y+size, x-size:x+size]
# # moment0_cutout.quicklook(use_aplpy=True)

# from astropy.modeling import models, fitting


# # Define the spatial grid for the fit centered at y, x = 32, 32
# yy, xx = moment0_cutout.spatial_coordinate_map

# # Define a single 2D Gaussian model.
# p_init_gauss2D = models.Gaussian2D(x_mean=xx[size, size], y_mean=yy[size, size],
#                                    x_stddev=8 * u.arcsec, y_stddev=8 * u.arcsec)


# # And fit with the Levenberg-Marquardt algorithm and least squares statistic.
# fit_p = fitting.LevMarLSQFitter()


# # Set NaNs to 0 for the fit.
# moment0_cutout_quant = moment0_cutout.quantity
# # moment0_cutout_quant = moment0_cutout.value
# moment0_cutout_quant[np.isnan(moment0_cutout_quant)] = 0.

# with warnings.catch_warnings():
#     # Ignore model linearity warning from the fitter
#     warnings.simplefilter('ignore')
#     p_gauss2D = fit_p(p_init_gauss2D, xx, yy, moment0_cutout_quant)
    
    
# p_gauss2D

# print(f"{p_gauss2D.x_stddev.quantity.to(u.arcsec)} by {p_gauss2D.y_stddev.quantity.to(u.arcsec)}")

# plt.figure(figsize=(18, 6))

# plt.subplot(1, 3, 1)
# plt.title("Image", fontsize=18)

# plt.imshow(moment0_cutout.value, origin='lower', cmap='inferno')
# plt.colorbar()
# plt.subplot(1, 3, 2)
# plt.title("Model", fontsize=18)

# plt.imshow(p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
# plt.colorbar()
# plt.subplot(1, 3, 3)
# plt.title("Residual", fontsize=18)

# plt.imshow(moment0_cutout.value - p_gauss2D(xx, yy).value, origin='lower', cmap='inferno', vmin=-0.7, vmax=0.7)
# plt.colorbar()


# plt.tight_layout()

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
Hbeta = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_Hbeta_Flux.fits")
Halpha = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_sincgauss_Halpha_Flux.fits")
OIII3726 = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OII3726_Flux.fits")
NII6584 = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_sincgauss_NII6583_Flux.fits")
OIII4959 = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OIII4959_Flux.fits")
OIII5007 = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OIII5007_Flux.fits")
#%%


# Step 1: Determine the observed fluxes for the wavelengths of interest
observed_OII3726 = OIII3726[0].data
observed_Halpha = Halpha[0].data
observed_Hbeta = Hbeta[0].data 
observed_NII6584 = NII6584[0].data 
observed_OIII4959 = OIII4959[0].data 
observed_OIII5007 = OIII5007[0].data
#%%
import numpy as np
from extinction import remove, fm07

def extinction(wavelength, flux):
    corrected_flux = remove(fm07(wavelength, 1.0), flux)
    return corrected_flux

extinction_corrected_OII3726 = []

for x in range(len(observed_OII3726)):
    for y in range(len(observed_OII3726[0])):
        extinction_corrected_flux = extinction(np.array([3756]), observed_OII3726[x][y])
        extinction_corrected_OII3726.append(extinction_corrected_flux)

extinction_corrected_OII3726 = np.array(extinction_corrected_OII3726).reshape(len(observed_OII3726), len(observed_OII3726[0]))

# Apply extinction correction for observed_Halpha
extinction_corrected_Halpha = []

for x in range(len(observed_Halpha)):
    for y in range(len(observed_Halpha[0])):
        extinction_corrected_flux = extinction(np.array([6562]), observed_Halpha[x][y])
        extinction_corrected_Halpha.append(extinction_corrected_flux)

# Reshape the list of arrays into a numpy array with shape (250, 287)
extinction_corrected_Halpha = np.array(extinction_corrected_Halpha).reshape(len(observed_Halpha), len(observed_Halpha[0]))

# Apply extinction correction for observed_Hbeta
extinction_corrected_Hbeta = []

for x in range(len(observed_Hbeta)):
    for y in range(len(observed_Hbeta[0])):
        extinction_corrected_flux = extinction(np.array([4861]), observed_Hbeta[x][y])
        extinction_corrected_Hbeta.append(extinction_corrected_flux)

# Reshape the list of arrays into a numpy array with shape (250, 287)
extinction_corrected_Hbeta = np.array(extinction_corrected_Hbeta).reshape(len(observed_Hbeta), len(observed_Hbeta[0]))

# Apply extinction correction for observed_NII6584
extinction_corrected_NII6584 = []

for x in range(len(observed_NII6584)):
    for y in range(len(observed_NII6584[0])):
        extinction_corrected_flux = extinction(np.array([6584]), observed_NII6584[x][y])
        extinction_corrected_NII6584.append(extinction_corrected_flux)

# Reshape the list of arrays into a numpy array with shape (250, 287)
extinction_corrected_NII6584 = np.array(extinction_corrected_NII6584).reshape(len(observed_NII6584), len(observed_NII6584[0]))

# Apply extinction correction for observed_OIII4959
extinction_corrected_OIII4959 = []

for x in range(len(observed_OIII4959)):
    for y in range(len(observed_OIII4959[0])):
        extinction_corrected_flux = extinction(np.array([4959]), observed_OIII4959[x][y])
        extinction_corrected_OIII4959.append(extinction_corrected_flux)

# Reshape the list of arrays into a numpy array with shape (250, 287)
extinction_corrected_OIII4959 = np.array(extinction_corrected_OIII4959).reshape(len(observed_OIII4959), len(observed_OIII4959[0]))

# Apply extinction correction for observed_OIII5007
extinction_corrected_OIII5007 = []

for x in range(len(observed_OIII5007)):
    for y in range(len(observed_OIII5007[0])):
        extinction_corrected_flux = extinction(np.array([5007]), observed_OIII5007[x][y])
        extinction_corrected_OIII5007.append(extinction_corrected_flux)

# Reshape the list of arrays into a numpy array with shape (250, 287)
extinction_corrected_OIII5007 = np.array(extinction_corrected_OIII5007).reshape(len(observed_OIII5007), len(observed_OIII5007[0]))

#%%

plt.imshow(observed_Halpha, vmin=1e-18, vmax=8e-16, origin='lower')
plt.imshow(extinction_corrected_Halpha, vmin=1e-18, vmax=8e-16, origin='lower')

#%%
#Creating flux maps that are corrected for extinction
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OII3726
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OII3726correctionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_Halpha
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('Halphacorrectionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_Hbeta
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('Hbetacorrectionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_NII6584
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('NII6584correctionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OIII4959
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OIII4959correctionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OIII5007
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OIII5007correctionbin2.fits', overwrite=True)  # to write just that HDU to a new file, or

#%% Same but for the uncertainty files

Hbeta_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_Hbeta_Flux_err.fits")
Halpha_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_sincgauss_Halpha_Flux_err.fits")
OIII3726_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OII3726_Flux_err.fits")
NII6584_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_sincgauss_NII6583_Flux_err.fits")
OIII4959_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OIII4959_Flux_err.fits")
OIII5007_err = fits.open("C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/Flux Maps/UGC2885_2_gaussian_OIII5007_Flux_err.fits")

observed_Hbeta_err = Hbeta_err[0].data
observed_Halpha_err = Halpha_err[0].data
observed_OII3726_err = OIII3726_err[0].data
observed_NII6584_err = NII6584_err[0].data
observed_OIII4959_err = OIII4959_err[0].data
observed_OIII5007_err = OIII5007_err[0].data

#%%

extinction_corrected_OII3726_err = []

for x in range(len(observed_OII3726_err)):
    for y in range(len(observed_OII3726_err[0])):
        extinction_corrected_flux_err = extinction(np.array([3756]), observed_OII3726_err[x][y])
        extinction_corrected_OII3726_err.append(extinction_corrected_flux_err)

extinction_corrected_OII3726_err = np.array(extinction_corrected_OII3726_err).reshape(len(observed_OII3726_err), len(observed_OII3726_err[0]))

extinction_corrected_Halpha_err = []

for x in range(len(observed_Halpha_err)):
    for y in range(len(observed_Halpha_err[0])):
        extinction_corrected_flux_err = extinction(np.array([6562]), observed_Halpha_err[x][y])
        extinction_corrected_Halpha_err.append(extinction_corrected_flux_err)

extinction_corrected_Halpha_err = np.array(extinction_corrected_Halpha_err).reshape(len(observed_Halpha_err), len(observed_Halpha_err[0]))

extinction_corrected_Hbeta_err = []

for x in range(len(observed_Hbeta_err)):
    for y in range(len(observed_Hbeta_err[0])):
        extinction_corrected_flux_err = extinction(np.array([4861]), observed_Hbeta_err[x][y])
        extinction_corrected_Hbeta_err.append(extinction_corrected_flux_err)

extinction_corrected_Hbeta_err = np.array(extinction_corrected_Hbeta_err).reshape(len(observed_Hbeta_err), len(observed_Hbeta_err[0]))

extinction_corrected_NII6584_err = []

for x in range(len(observed_NII6584_err)):
    for y in range(len(observed_NII6584_err[0])):
        extinction_corrected_flux_err = extinction(np.array([6584]), observed_NII6584_err[x][y])
        extinction_corrected_NII6584_err.append(extinction_corrected_flux_err)

extinction_corrected_NII6584_err = np.array(extinction_corrected_NII6584_err).reshape(len(observed_NII6584_err), len(observed_NII6584_err[0]))

extinction_corrected_OIII4959_err = []

for x in range(len(observed_OIII4959_err)):
    for y in range(len(observed_OIII4959_err[0])):
        extinction_corrected_flux_err = extinction(np.array([4959]), observed_OIII4959_err[x][y])
        extinction_corrected_OIII4959_err.append(extinction_corrected_flux_err)

extinction_corrected_OIII4959_err = np.array(extinction_corrected_OIII4959_err).reshape(len(observed_OIII4959_err), len(observed_OIII4959_err[0]))

extinction_corrected_OIII5007_err = []

for x in range(len(observed_OIII5007_err)):
    for y in range(len(observed_OIII5007_err[0])):
        extinction_corrected_flux_err = extinction(np.array([5007]), observed_OIII5007_err[x][y])
        extinction_corrected_OIII5007_err.append(extinction_corrected_flux_err)

extinction_corrected_OIII5007_err = np.array(extinction_corrected_OIII5007_err).reshape(len(observed_OIII5007_err), len(observed_OIII5007_err[0]))

#%%

with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OII3726_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OII3726correctionbin2_err.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_Halpha_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('Halphacorrectionbin2_err.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_Hbeta_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('Hbetacorrectionbin2_err.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_NII6584_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('NII6584correctionbin2_err.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OIII4959_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OIII4959correctionbin2_err.fits', overwrite=True)  # to write just that HDU to a new file, or
with fits.open('Header.fits') as hdu_list:
    hdu = hdu_list[0]
    hdu.data = extinction_corrected_OIII5007_err
    # or hdu.data[<some index>] = <some value> i.e. just directly modify the existing array
    hdu.writeto('OIII5007correctionbin2_err.fits', overwrite=True)

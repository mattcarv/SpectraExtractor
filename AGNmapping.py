import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

hdul =  fits.open('C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_NII6583_Flux.fits')
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


#%% Using Bavi's definition of central region
import matplotlib.patches as patches
from astropy.coordinates import SkyCoord

fig, ax = plt.subplots(subplot_kw={'projection': wcs})

ax.imshow(data, cmap='viridis', origin='lower', vmax=(.001*data.max()))

# wcs_center = SkyCoord('03h53m02.4811s', '+35d35m22.103s', frame='icrs')
# g_ell_center_small = wcs.world_to_pixel(wcs_center)
# x_center, y_center = (125, 128)
g_ell_center_small = (124, 128)
g_ell_width_small = 11.743836 * 2
g_ell_height_small = 7.3121888 * 2
angle_small = 315.71316

ellipse = patches.Ellipse(g_ell_center_small, g_ell_width_small, g_ell_height_small, angle=angle_small, fill=False, edgecolor='red', linewidth=2)
ax.add_patch(ellipse)


ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')

plt.show()

#------------------------------------------------------------------------------
x, y = np.meshgrid(np.arange(hdul[0].data.shape[1]), np.arange(hdul[0].data.shape[0]))

cos_angle = np.cos(np.radians(180. - angle_small))
sin_angle = np.sin(np.radians(180. - angle_small))

xc = x - g_ell_center_small[0]
yc = y - g_ell_center_small[1]

xct = xc * cos_angle - yc * sin_angle
yct = xc * sin_angle + yc * cos_angle

rad_cc = (xct**2 / (g_ell_width_small / 2.)**2) + (yct**2 / (g_ell_height_small / 2.)**2)

pix_ell = np.sum(rad_cc <= 1.0)
pix_ell_values = hdul[0].data[rad_cc <= 1.0]


print(f"Number of pixels inside the ellipse: {pix_ell}")
print(pix_ell_values.mean())
#%% Holwerda's definition of the central region
import matplotlib.patches as patches
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u

fig, ax = plt.subplots(subplot_kw={'projection': wcs})

ax.imshow(data, cmap='viridis', origin='lower', vmax=(.001*data.max()))

# wcs_center = SkyCoord('03h53m02.4811s', '+35d35m22.103s', frame='icrs')
# g_ell_center_small = wcs.world_to_pixel(wcs_center)
g_circ_center_small = (124, 128)

diameter_arcseconds = 6.1

pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel

diameter_degrees = diameter_arcseconds * u.arcsec.to(u.deg)

diameter_pixels = (diameter_degrees / pixel_scale).value

circle = patches.Circle((g_circ_center_small[0], g_circ_center_small[1]), radius=diameter_pixels / 2, fill=False, edgecolor='red', linewidth=2)
ax.add_patch(circle)

ax.coords['ra'].set_axislabel('Right Ascension (J2000)')
ax.coords['dec'].set_axislabel('Declination (J2000)')
plt.show()

#------------------------------------------------------------------------------
xc = x - g_circ_center_small[0]
yc = y - g_circ_center_small[1]

rad_cc_circle = (xc**2 + yc**2) <= (diameter_pixels / 2)**2

pix_circle = np.sum(rad_cc_circle)
pix_circle_values = data[rad_cc_circle]

print(f"Number of pixels inside the circle: {pix_circle}")
print(pix_circle_values.mean())
#%%

plt.hist(pix_ell_values, bins=100)
plt.xlabel('Flux')
plt.ylabel('Count')
plt.show()

#%%

plt.hist(pix_circle_values, bins=100)
plt.xlabel('Flux')
plt.ylabel('Count')
plt.show()

#%% BPT Boundaries (from Kauffmann, 2003)

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
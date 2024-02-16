import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

hdul =  fits.open('C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_SII6716_Flux.fits')
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


#%% Using Bavi's definition of central region - ELLIPSE
import matplotlib.patches as patches
from astropy.coordinates import SkyCoord

fig, ax = plt.subplots(subplot_kw={'projection': wcs})

ax.imshow(data, cmap='viridis', origin='lower', vmax=(.001*data.max()))

# wcs_center = SkyCoord('03h53m02.4811s', '+35d35m22.103s', frame='icrs')
# g_ell_center_small = wcs.world_to_pixel(wcs_center)
# x_center, y_center = (125, 128)
g_ell_center_small = (122.05, 126.35)
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
#print(pix_ell_values)

#%% Hist2D
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [15, 15]

ellipse_pixel_data = np.zeros(pix_ell, dtype=[('x', int), ('y', int), ('value', float)])
ellipse_pixel_data['x'], ellipse_pixel_data['y'] = np.where(rad_cc <= 1.0)
ellipse_pixel_data['value'] = pix_ell_values

def gaussian_2d(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2) / (2 * sigma_x**2) + (np.sin(theta)**2) / (2 * sigma_y**2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (4 * sigma_y**2)
    c = (np.sin(theta)**2) / (2 * sigma_x**2) + (np.cos(theta)**2) / (2 * sigma_y**2)
    g = offset + amplitude * np.exp(- (a * (x - xo)**2 + 2 * b * (x - xo) * (y - yo) + c * (y - yo)**2))
    return g.ravel()

initial_guess = (1, ellipse_pixel_data['x'].mean(), ellipse_pixel_data['y'].mean(),
                  8, 8, 0, 0)
popt, _ = curve_fit(gaussian_2d, (ellipse_pixel_data['x'], ellipse_pixel_data['y']),
                    ellipse_pixel_data['value'], p0=initial_guess)


x_fit = np.linspace(ellipse_pixel_data['x'].min(), ellipse_pixel_data['x'].max(), 100)
y_fit = np.linspace(ellipse_pixel_data['y'].min(), ellipse_pixel_data['y'].max(), 100)
x_fit, y_fit = np.meshgrid(x_fit, y_fit)

z_fit = gaussian_2d((x_fit, y_fit), *popt).reshape(x_fit.shape)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D surface
ax.plot_trisurf(ellipse_pixel_data['x'], ellipse_pixel_data['y'], ellipse_pixel_data['value'],
                cmap='plasma', edgecolor='k', alpha=0.5)

ax.plot_surface(x_fit, y_fit, z_fit, color='green', alpha=0.3)

# Set axis labels
ax.set_xlabel('X Pixel')
ax.set_ylabel('Y Pixel')
ax.set_zlabel('Pixel Value')

ax.view_init(elev=20, azim=30)


plt.show()

x_peak, y_peak = popt[1], popt[2]
print(f"Peak Position (x, y): ({x_peak}, {y_peak})")

#%% Holwerda's definition of the central region
import matplotlib.patches as patches
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u

fig, ax = plt.subplots(subplot_kw={'projection': wcs})

ax.imshow(data, cmap='viridis', origin='lower', vmax=(.001*data.max()))

# wcs_center = SkyCoord('03h53m02.4811s', '+35d35m22.103s', frame='icrs')
# g_ell_center_small = wcs.world_to_pixel(wcs_center)
g_circ_center_small = (122.05, 126.35)

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
circle_pixel_data = np.zeros(pix_circle, dtype=[('x', int), ('y', int), ('value', float)])
circle_pixel_data['x'], circle_pixel_data['y'] = np.where(rad_cc_circle)
circle_pixel_data['value'] = pix_circle_values

def gaussian_2d(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2) / (2 * sigma_x**2) + (np.sin(theta)**2) / (2 * sigma_y**2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (4 * sigma_y**2)
    c = (np.sin(theta)**2) / (2 * sigma_x**2) + (np.cos(theta)**2) / (2 * sigma_y**2)
    g = offset + amplitude * np.exp(- (a * (x - xo)**2 + 2 * b * (x - xo) * (y - yo) + c * (y - yo)**2))
    return g.ravel()

initial_guess = (1, circle_pixel_data['x'].mean(), circle_pixel_data['y'].mean(),
                  4, 5, 0, 0)
popt, _ = curve_fit(gaussian_2d, (circle_pixel_data['x'], circle_pixel_data['y']),
                    circle_pixel_data['value'], p0=initial_guess)


x_fit = np.linspace(circle_pixel_data['x'].min(), circle_pixel_data['x'].max(), 100)
y_fit = np.linspace(circle_pixel_data['y'].min(), circle_pixel_data['y'].max(), 100)
x_fit, y_fit = np.meshgrid(x_fit, y_fit)

z_fit = gaussian_2d((x_fit, y_fit), *popt).reshape(x_fit.shape)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D surface
ax.plot_trisurf(circle_pixel_data['x'], circle_pixel_data['y'], circle_pixel_data['value'],
                cmap='plasma', edgecolor='k', alpha=0.5)

# Plot the fitted Gaussian surface
ax.plot_surface(x_fit, y_fit, z_fit, color='green', alpha=0.3)

ax.set_xlabel('X Pixel')
ax.set_ylabel('Y Pixel')
ax.set_zlabel('Pixel Value')

ax.view_init(elev=10, azim=80)

plt.show()

x_peak, y_peak = popt[2], popt[1]
print(f"Peak Position (x, y): ({x_peak}, {y_peak})")


#%% For-loop to run through many files and compute the flux from the same region
# as defined above: BASED ON HOLWERDA'S CIRCLE
file_names = [
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/Rubin_SN2_FINAL_NOSTOCHS_forpaper_2_OIII5007_Flux.fits',
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/Rubin_SN2_FINAL_NOSTOCHS_forpaper_2_Hbeta_Flux.fits',
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_NII6583_Flux.fits',
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_Halpha_Flux.fits',
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_SII6716_Flux.fits',
    'C:/Users/mathe/Downloads/LUCI Files-20240207T161845Z-001/LUCI Files/RubinsGalaxy_withuncer_finalNOSTOCHS_forpaper_2_SII6731_Flux.fits'
]

list_mean_circ = []



def FluxCalc (data):
    g_circ_center_small = (122.05, 126.35)
    diameter_arcseconds = 6.1
    pixel_scale = proj_plane_pixel_scales(wcs)[0] * u.deg / u.pixel
    diameter_degrees = diameter_arcseconds * u.arcsec.to(u.deg)
    diameter_pixels = (diameter_degrees / pixel_scale).value
    
    x, y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    
    xc = x - g_circ_center_small[0]
    yc = y - g_circ_center_small[1]
    rad_cc_circle = (xc**2 + yc**2) <= (diameter_pixels / 2)**2
    
    # pix_circle = np.sum(rad_cc_circle)
    
    pix_circle_values = data[rad_cc_circle]
    
    return(pix_circle_values.mean())

for file_name in file_names:
    
    hdul = fits.open(file_name)
    data = hdul[0].data
    
    wcs = WCS(hdul[0].header)
    
    list_mean_circ.append(FluxCalc(data))
    
print(list_mean_circ)

logOIII_Hbeta_circ = np.log10(list_mean_circ[0]/list_mean_circ[1])
logNII_Halpha_circ = np.log10(list_mean_circ[2]/list_mean_circ[3])
logSII_Halpha_circ = np.log10((list_mean_circ[4]+list_mean_circ[5])/list_mean_circ[3])
#%%
# DO THE SAME WITH THE ELLIPSE TO COMPARE
list_mean_ell = []



def FluxCalc (data):
    g_ell_center_small = (122.05, 126.35)
    g_ell_width_small = 11.743836 * 2
    g_ell_height_small = 7.3121888 * 2
    angle_small = 315.71316
    
    x, y = np.meshgrid(np.arange(hdul[0].data.shape[1]), np.arange(hdul[0].data.shape[0]))

    cos_angle = np.cos(np.radians(180. - angle_small))
    sin_angle = np.sin(np.radians(180. - angle_small))

    xc = x - g_ell_center_small[0]
    yc = y - g_ell_center_small[1]

    xct = xc * cos_angle - yc * sin_angle
    yct = xc * sin_angle + yc * cos_angle

    rad_cc = (xct**2 / (g_ell_width_small / 2.)**2) + (yct**2 / (g_ell_height_small / 2.)**2)
    #pix_ell = np.sum(rad_cc <= 1.0)
    pix_ell_values = hdul[0].data[rad_cc <= 1.0]
    
    return(pix_ell_values.mean())

for file_name in file_names:
    
    hdul = fits.open(file_name)
    data = hdul[0].data
    
    wcs = WCS(hdul[0].header)
    
    list_mean_ell.append(FluxCalc(data))
    
print(list_mean_ell)

logOIII_Hbeta_ell = np.log10(list_mean_ell[0]/list_mean_ell[1])
logNII_Halpha_ell = np.log10(list_mean_ell[2]/list_mean_ell[3])
logSII_Halpha_ell = np.log10((list_mean_ell[4]+list_mean_ell[5])/list_mean_circ[3])

#%% BPT Boundaries (from KEWLEY, 2006) FOR OIII/HBETA VERSUS NII/HALPHA
plt.rcParams["figure.figsize"] = [10, 8]

# Define the function
# (0.61 / (x - 0.47)) + 1.19
def BPT_NII(x):
    return (0.61 / (x - 0.05)) + 1.3

def BPT_NII_Comp(x):
    return (0.61 / (x - 0.47)) + 1.19

# Generate x values
x_values = np.linspace(-1.6, 0, 100)
x_values_comp = np.linspace(-1.6, 0.25, 100)

# Calculate corresponding y values
y_values = BPT_NII(x_values)
y_values_comp = BPT_NII_Comp(x_values_comp)

# Plot the function
plt.plot(x_values, y_values, c='k', linestyle='-')
plt.plot(x_values_comp, y_values_comp, c='k', linestyle='--')
plt.scatter(logNII_Halpha_circ, logOIII_Hbeta_circ, c='r', marker='d', s=100,
            label='Circle')
plt.scatter(logNII_Halpha_ell, logOIII_Hbeta_ell, c='k', marker='d', s=100,
            label='Ellipse')
plt.text(-1.15, -0.3, 'STAR-FORMING', fontsize=15)
plt.text(0,0.8, 'AGN', fontsize=15)
plt.text(-0.2, -0.7, 'Composite', fontsize=15)
plt.xlim(-1.5, 0.5)
plt.ylim(-1.2, 1.5)
plt.xlabel(r'log ($NII\lambda 6584/ H\alpha$)')
plt.ylabel(r'log ($OIII\lambda 5007/ H\beta$)')

plt.legend(loc=2)
plt.show()

#%% (KEWLEY, 2006) FOR OIII/HBETA VERSUS SII/HALPHA

# Define the function
def BPT_SII(x):
    return (0.72 / (x - 0.32)) + 1.3

def BPT_SII_LINER(x):
    return (1.89*x) + 0.76

# Generate x values
x_values = np.linspace(-1.6, 0.15, 100)
x_values_liner = np.linspace(-0.31, 0.6, 100)
# Calculate corresponding y values
y_values = BPT_SII(x_values)
y_values_liner = BPT_SII_LINER(x_values_liner)

# Plot the function
plt.plot(x_values, y_values, c='k', linestyle='--')
plt.plot(x_values_liner, y_values_liner, c='k', linestyle='--')
plt.scatter(logSII_Halpha_circ, logOIII_Hbeta_circ, c='r', marker='d', s=100,
            label='Circle')
plt.scatter(logSII_Halpha_ell, logOIII_Hbeta_ell, c='k', marker='d', s=100,
            label='Ellipse')
plt.text(-1.15, -0.3, 'STAR-FORMING', fontsize=15)
plt.text(0.2, -0.2, 'LINER', fontsize=15)
plt.text(-0.7, 1, 'Seyfert', fontsize=15)
plt.xlim(-1.5, 0.5)
plt.ylim(-1.2, 1.5)
plt.xlabel(r'log ($SII\lambda \lambda 6716, 31/ H\alpha$)')
plt.ylabel(r'log ($OIII\lambda 5007/ H\beta$)')

plt.legend(loc=2)
plt.show()

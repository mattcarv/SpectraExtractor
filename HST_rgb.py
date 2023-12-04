import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits

# Read in the three images downloaded from here:
g = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')[0]
    
g_data = fits.open('/home/mdocarm/Downloads/f475_mosaic.fits')[0].data*2.4
r = fits.open('/home/mdocarm/Downloads/f606_reprojected_final.fits')[0].data
i = fits.open('/home/mdocarm/Downloads/f814_reprojected_final.fits')[0].data*1.1

#%%

rgb_default = make_lupton_rgb(i, r, g_data, Q=5, stretch=0.1,
                              filename="HST_rgb.png")
plt.imshow(rgb_default, origin='lower')

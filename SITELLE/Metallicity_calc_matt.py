import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
fig,ax = plt.subplots(1)
ax.set_aspect('equal')
    
def isPointInEllpise(x,y, g_ell_center, g_ell_width, g_ell_height, angle):

    g_ellipse = patches.Ellipse(g_ell_center, g_ell_width, g_ell_height, angle=angle, fill=False, edgecolor='green', linewidth=2)
    ax.add_patch(g_ellipse)

    cos_angle = np.cos(np.radians(180.-angle))
    sin_angle = np.sin(np.radians(180.-angle))

    xc = x - g_ell_center[0]
    yc = y - g_ell_center[1]

    xct = xc * cos_angle - yc * sin_angle
    yct = xc * sin_angle + yc * cos_angle 

    rad_cc = (xct**2/(g_ell_width/2.)**2) + (yct**2/(g_ell_height/2.)**2)

    return rad_cc <= 1. 
    # Set the colors. Black if outside the ellipse, green if inside
    # colors_array = np.array(['black'] * len(rad_cc))
    # colors_array[np.where()[0]] = 'green'
    
pointsInTheEllipseX_small = []
pointsInTheEllipseY_small = []

pointsInTheEllipseX_medium = []
pointsInTheEllipseY_medium = []

pointsInTheEllipseX_large = []
pointsInTheEllipseY_large = []

entire_galaxyX = []
entire_galaxyY = []

#Storing the Pixels from Each Eclipse in their respective list
for x in range(0,166):
    for y in range (0,191):
        if(isPointInEllpise(x,y, (84.255593, 84.442971), 11.743836*2, 7.3121888*2, 315.71316)):
            pointsInTheEllipseX_small.append(x)
            pointsInTheEllipseY_small.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
        elif(isPointInEllpise(x,y, (84.255593,84.442971), 36.746757*2, 18.63428*2, 315.71316)):
            pointsInTheEllipseX_medium.append(x)
            pointsInTheEllipseY_medium.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
        elif(isPointInEllpise(x,y, (90.016614,79.125106), 75.909136*2, 43.399425*2, 315.71316) and not isPointInEllpise(x,y, (48.360789,112.14037), 6.8689095*2, 7.312065*2, 0.)):
            pointsInTheEllipseX_large.append(x)
            pointsInTheEllipseY_large.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
            #add point to things that you should plot
        

#ax.scatter(pointsInTheEllipseX_small, pointsInTheEllipseY_small)
#ax.scatter(pointsInTheEllipseX_medium, pointsInTheEllipseY_medium)
#ax.scatter(pointsInTheEllipseX_large, pointsInTheEllipseY_large)
#ax.scatter(entire_galaxyX, entire_galaxyY)

#%%

#plt.show()
#Choose the ratio map you want to find the metallicity from
image_data = fits.getdata('C:/Users/mathe/Downloads/Bavi-Files/UGC-2885-LUCI-Files/R23Uncertainities_new.fits')
#Calculating metallicity of smallest ellipse

n = 0 
y = []
m = 0 
for i in range(0, len(pointsInTheEllipseX_small)):
    x = pointsInTheEllipseX_small[i]
    y = pointsInTheEllipseY_small[i]
    # There exsists some nans and inf in the flux maps that we have chosen to remove
    if np.isnan(image_data[x][y]): 
        continue
    if np.isinf(image_data[x][y]):
        continue
    else: 
        #print(image_data[x][y])
        #Adding all the values from the small ellipse
        n = n + image_data[x][y]
        m = m + 1

#Calculating metallicity
small_ellipse_metallicity = n/m 
print(small_ellipse_metallicity)
#Calculating metallicity from medium ellipse

#%%

n = 0 
y = []
m = 0 

for i in range(0, len(pointsInTheEllipseX_medium)):
    x = pointsInTheEllipseX_medium[i]
    y = pointsInTheEllipseY_medium[i]
    if np.isnan(image_data[x][y]): 
        continue
    if np.isinf(image_data[x][y]):
        continue
    if (image_data[x][y]) <= 0:
        continue
    else: 
        #print(image_data[x][y])
        n = n + image_data[x][y]
        m = m + 1

medium_ellipse_metallicity = n/m
print(medium_ellipse_metallicity)
#Calculating metalliicty from largest ellipse

#%%

n = 0 
y = []
m = 0 

for i in range(0, len(pointsInTheEllipseX_large)):
    x = pointsInTheEllipseX_large[i]
    y = pointsInTheEllipseY_large[i]
    if np.isnan(image_data[x][y]): 
        continue
    if np.isinf(image_data[x][y]):
        continue
    else: 
        #print(image_data[x][y])
        n = n + image_data[x][y]
        m = m + 1

large_ellipse_metallicity = n/m
print(large_ellipse_metallicity)

#Calculating the global metallicity 
#%%
n = 0 
y = []
m = 0 
for i in range(0, len(entire_galaxyX)):
    x = entire_galaxyX[i]
    y = entire_galaxyY[i]
    if np.isnan(image_data[x][y]): 
        continue
    if np.isinf(image_data[x][y]):
        continue
    else: 
        n = n + image_data[x][y]
        m = m + 1

entire_ellipse_metallicity = n/m
print(entire_ellipse_metallicity)

#%%


from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

hdulist = fits.open('UGC2885_bin2_2_sincgauss_Halpha_Flux.fits') 
hdu = hdulist[0] 
wcs = WCS(hdu.header) 

pointsInTheEllipseX = []
pointsInTheEllipseY = []

for x in range(0,166):
    for y in range (0,191):
        if(isPointInEllpise(x,y, (84.255593, 84.442971), 11.743836*2, 7.3121888*2, 315.71316)):
            pointsInTheEllipseX_small.append(x)
            pointsInTheEllipseY_small.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
        elif(isPointInEllpise(x,y, (84.255593,84.442971), 36.746757*2, 18.63428*2, 315.71316)):
            pointsInTheEllipseX_medium.append(x)
            pointsInTheEllipseY_medium.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
        elif(isPointInEllpise(x,y, (90.016614,79.125106), 75.909136*2, 43.399425*2, 315.71316) and not isPointInEllpise(x,y, (48.360789,112.14037), 6.8689095*2, 7.312065*2, 0.)):
            pointsInTheEllipseX_large.append(x)
            pointsInTheEllipseY_large.append(y)
            entire_galaxyX.append(x)
            entire_galaxyY.append(y)
            #add point to things that you should plot
            
            
plt.subplot(projection=wcs)
#plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
plt.imshow(hdu.data, origin='lower')
plt.grid(color='white', ls='solid')
ellipse = Ellipse((84.255593, 84.442971), 11.743836*2, 7.3121888*2, angle=315.71316, fill=False)
ellipse_medium = Ellipse((84.255593,84.442971), 36.746757*2, 18.63428*2, angle=315.71316, fill=False)
ellipse_large = Ellipse((90.016614,79.125106), 75.909136*2, 43.399425*2, angle=315.71316, fill=False)
ellipse_circle = Ellipse((48.360789,112.14037), 6.8689095*2, 7.312065*2, angle=0, fill=False)
plt.gca().add_patch(ellipse);
plt.gca().add_patch(ellipse_medium);
plt.gca().add_patch(ellipse_large);
plt.gca().add_patch(ellipse_circle);
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')

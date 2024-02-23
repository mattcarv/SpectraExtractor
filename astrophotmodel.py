import os
import astrophot as ap
import numpy as np
import torch
from astropy.io import fits
import matplotlib.pyplot as plt
from time import time

hdu = fits.open('/home/mdocarm/Downloads/HST-idne08050-idne08050-PRODUCT/HST/product/idne08050_drc.fits')
target_data = np.array(hdu[1].data, dtype=np.float64)


target = ap.image.Target_Image(
    data = target_data,
    pixelscale = 0.039, 
    zeropoint = 26.04, 
    #variance = np.ones(target_data.shape)/1e3,# set the variance for this image (in general it should be more accurate than this)
)

model3 = ap.models.AstroPhot_Model(
    name = "model with target", 
    model_type = "spline galaxy model",
    target = target,
    window = [[9, 700],[41, 580]],
)

fig6, ax6 = plt.subplots(figsize = (8,8))
ap.plots.target_image(fig6, ax6, model3.target)
ap.plots.model_window(fig6, ax6, model3)
plt.clf()

model3.initialize()

result = ap.fit.LM(model3, verbose = 1).fit()
print(result.message)

fig7, ax7 = plt.subplots(1, 2, figsize = (16,6))
ap.plots.model_image(fig7, ax7[0], model3)
ap.plots.residual_image(fig7, ax7[1], model3)
plt.show()

fig10, ax10 = plt.subplots(figsize = (8,8))
ap.plots.galaxy_light_profile(fig10, ax10, model3)
ap.plots.radial_median_profile(fig10, ax10, model3)
plt.show()

result.update_uncertainty()
for P in model3.parameters:
    print(f"parameter {P.name} is: {P.value.detach().cpu().tolist()} +- {P.uncertainty.detach().cpu().tolist()}")
    
    
fig, ax = ap.plots.covariance_matrix(result.covariance_matrix.detach().cpu().numpy(), 
                                      model3.parameters.get_vector().detach().cpu().numpy(), 
                                      model3.parameters.get_name_vector()
)
plt.show()
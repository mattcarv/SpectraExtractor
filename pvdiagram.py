import pylab as pl
import numpy as np
from astropy.visualization import quantity_support
from astropy import units as u
from astropy import wcs
import astropy.io.fits as fits

# set so that these display properly on black backgrounds
pl.rcParams['figure.facecolor']='w'

from spectral_cube import SpectralCube

from pvextractor import extract_pv_slice, Path

cube = SpectralCube.read('maskedcube.fits')
header = cube[0].header
data = cube[0].data

# pl.imshow(cube[40].value, origin='lower')



path = Path([(27,13), (53,37)], width=10*u.arcsec)

# ax = pl.subplot(111, projection=cube.wcs.celestial)
# ax.imshow(cube[40].value)
# path.show_on_axis(ax, spacing=1, color='r')

pvdiagram = extract_pv_slice(cube=cube, path=path, spacing=1)
pvdiagram

# ww = wcs.WCS(pvdiagram.header)

# ax = pl.subplot(111, projection=ww)
# im = ax.imshow(pvdiagram.data)
# cb = pl.colorbar(mappable=im)
# cb.set_label("Brightness Temperature [K]")

# ax0 = ax.coords[0]
# ax0.set_format_unit(u.arcmin)
# ax1 = ax.coords[1]
# ax1.set_format_unit(u.km/u.s)

# ax.set_ylabel("Velocity [km/s]")
# ax.set_xlabel("Offset [arcmin]")

mx = cube.max(axis=0).value

pl.figure(figsize=(12,6))
ax = pl.subplot(121, projection=cube.wcs.celestial)
ax.imshow(mx)
path.show_on_axis(ax, spacing=1, edgecolor='r', linewidth=0.75, linestyle=':')

ww = wcs.WCS(pvdiagram.header)
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")

ax = pl.subplot(122, projection=ww)
im = ax.imshow(pvdiagram.data)
ax.set_aspect(0.5)

cb = pl.colorbar(mappable=im)
cb.set_label("Brightness Temperature [K]")

ax0 = ax.coords[0]
ax0.set_format_unit(u.arcmin)
ax1 = ax.coords[1]
ax1.set_format_unit(u.km/u.s)


ax.set_ylabel("Velocity [km/s]")
ax.set_xlabel("Offset [arcmin]")



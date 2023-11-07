import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})


dat = Table.read('/home/mdocarm/Downloads/xCOLDGASS_PubCat.fits')
df = dat.to_pandas()

# SELECTION

u2885 = pd.DataFrame({'LOGMSTAR': [11.68], 'LOGSFR_BEST': [0.21], 'LOGTDEP': [11.04]})

df = df[df.LOGMH2>0]

tdep = df.LOGMH2-df.LOGSFR_BEST
df = df.assign(LOGTDEP = tdep)

df = pd.concat([df, u2885])

def curve_function(x, a, b, c):
    return a * x - b * x**2 + c * x**3

filtered_df = df[df.LOGMSTAR.apply(np.isfinite) & df.LOGSFR_BEST.apply(np.isfinite)]
filtered_df_curve = filtered_df[filtered_df.LOGMSTAR < 11.5]

x = filtered_df.LOGMSTAR
y = filtered_df.LOGSFR_BEST
params, _ = curve_fit(curve_function, filtered_df_curve.LOGMSTAR, filtered_df_curve.LOGSFR_BEST)

x_fit = np.linspace(min(x), 11.4, 100)
y_fit = curve_function(x_fit, *params)

y_dotted = curve_function(x_fit, params[0], params[1], params[2])

plt.scatter(x, y, c=filtered_df.LOGTDEP, cmap='cool', vmax=11.5)
plt.text(11.38, 0.3, 'UGC 2885', c='black', fontsize=12)
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('log $t_{dep}\; (yr)$')

print("Fitted Parameters:")
print("a:", params[0])
print("b:", params[1])
print("c:", params[2])
plt.show()
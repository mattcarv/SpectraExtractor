import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})

dat = Table.read('/home/mdocarm/Downloads/lisenfeld_ss.csv')
df = dat.to_pandas()

u2885 = pd.DataFrame({'logMstar': [11.68], 'logSFR': [0.21], 'logMmolMstar': [-0.41]})

logMmolMstar = df.logMmol-df.logMstar
df = df.assign(logMmolMstar = logMmolMstar)

df = pd.concat([df, u2885])


def curve_function(x, a, b, c):
    return a * x - b * x**2 + c * x**3

filtered_df = df[df.logMstar.apply(np.isfinite) & df.logSFR.apply(np.isfinite)]
filtered_df_curve = filtered_df[filtered_df.logMstar < 11.68]

x = filtered_df.logMstar
y = filtered_df.logSFR
params, _ = curve_fit(curve_function, filtered_df_curve.logMstar, filtered_df_curve.logSFR)

x_fit = np.linspace(min(x), 11.68, 100)
y_fit = curve_function(x_fit, *params)

y_dotted = curve_function(x_fit, params[0], params[1], params[2])

plt.scatter(df.logMstar, df.logSFR, c=df.logMmolMstar, cmap='PuBu')
plt.text(11.63, 0.25, 'UGC 2885', c='black', fontsize=12)
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('log $M_{H_2}/M_{*}$')

plt.show()
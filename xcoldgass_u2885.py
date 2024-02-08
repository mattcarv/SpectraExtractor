import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})


#dat = Table.read('/home/mdocarm/Downloads/xCOLDGASS_PubCat.fits')
dat = Table.read('C:/Users/mathe/Downloads/xcoldgass_hi')
df = dat.to_pandas()


# SELECTION

u2885 = pd.DataFrame({'LOGMSTAR': [11.68], 'LOGSFR_BEST': [0.21], 'LOGTDEP_H2': [11.04],
                      'LOGTDEP_HI': [10.6], 'Z_PP04_O3N2': [9.05], 'Z_PP04_N2': [9.21],
                      'LOGMH2': [11.27]})

df = df[df.LOGMH2>0]

tdep_h2 = df.LOGMH2-df.LOGSFR_BEST
tdep_hi = df.logMH-df.LOGSFR_BEST
df = df.assign(LOGTDEP_HI = tdep_hi)
df = df.assign(LOGTDEP_H2 = tdep_h2)

df = pd.concat([df, u2885])


def curve_function(x, a, b, c):
    return a * x - b * x**2 + c * x**3

filtered_df = df[df.LOGMSTAR.apply(np.isfinite) & df.LOGSFR_BEST.apply(np.isfinite)]
filtered_df_curve = filtered_df[filtered_df.LOGMSTAR < 11.5]

x = filtered_df_curve.LOGMSTAR
y = filtered_df_curve.LOGSFR_BEST
params, _ = curve_fit(curve_function, x, y)

x_fit = np.linspace(min(x), 11.4, 100)
y_fit = curve_function(x_fit, *params)

y_dotted = curve_function(x_fit, params[0], params[1], params[2])
#%%
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})


fig, ax = plt.subplots()
plt.scatter(x, y, c=filtered_df.LOGTDEP_H2, cmap='cool', vmax=11)
plt.text(11.38, 0.3, 'UGC 2885', c='black', 
         fontsize=12)
plt.text(11.4, 1.4, 'H$_{2}$', c='k', fontsize=20, bbox=dict(facecolor='none',
                                                             edgecolor='red'))
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar(ticks=[7.8,8.8, 9.8, 10.8])
cbar.set_label('log $t_{dep}\; (yr)$')

# print("Fitted Parameters:")
print("a:", params[0])
print("b:", params[1])
print("c:", params[2])
plt.show()


plt.scatter(x, y, c=filtered_df.LOGTDEP_HI, cmap='winter')
plt.text(11.38, 0.3, 'UGC 2885', c='black', fontsize=12)
plt.text(11.4, 1.4, 'HI', c='k', fontsize=20, bbox=dict(facecolor='none',
                                                             edgecolor='red'))
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('log $t_{dep}\; (yr)$')

# print("Fitted Parameters:")
# print("a:", params[0])
# print("b:", params[1])
# print("c:", params[2])
plt.show()
#%%
# Plotting against metallicity

plt.scatter(x, y, c=filtered_df.Z_PP04_O3N2, cmap='Reds')
plt.text(11.38, 0.3, 'UGC 2885', c='black', fontsize=12)
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('12 + $log(O/H)$')

# print("Fitted Parameters:")
# print("a:", params[0])
# print("b:", params[1])
# print("c:", params[2])
plt.clf()

plt.scatter(x, y, c=filtered_df.Z_PP04_N2, cmap='Reds')
plt.text(11.38, 0.3, 'UGC 2885', c='black', fontsize=12)
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('12 + $log(O/H)$')

# print("Fitted Parameters:")
# print("a:", params[0])
# print("b:", params[1])
# print("c:", params[2])
plt.show()

#%%
# Plotting the MGMS

def linear_function(x, a, b):
    return a * x + b

filtered_df = df[df.LOGMSTAR.apply(np.isfinite) & df.LOGSFR_BEST.apply(np.isfinite)]
filtered_df_curve = filtered_df[filtered_df.LOGMSTAR < 11.5]

x = filtered_df.LOGMSTAR
y = filtered_df.LOGMH2
params_linear = np.polyfit(filtered_df_curve.LOGMSTAR, filtered_df_curve.LOGMH2, deg=1)

x_fit = np.linspace(min(x), 11.4, 100)
y_fit = linear_function(x_fit, *params_linear)


plt.scatter(x, filtered_df.LOGMH2, c=filtered_df.LOGSFR_BEST, cmap='magma')
plt.text(11.38, 11.12, 'UGC 2885', c='black', fontsize=12)
plt.ylim(7.4, 11.4)
plt.plot(x_fit, y_fit, 'k', linewidth=2)
plt.plot(x_fit, y_dotted+0.4, 'k-.', linewidth=1)
plt.plot(x_fit, y_dotted-0.4, 'k-.', linewidth=1)
plt.ylabel('log Molecular Gas Mass ($M_{\odot}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')

cbar = plt.colorbar()
cbar.set_label('log SFR ($M_{\odot} \; yr^{-1})$')

plt.show()
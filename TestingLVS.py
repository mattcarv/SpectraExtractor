#NEDLVS SAMPLE TESTING
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from astropy.table import Table
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams.update({'font.size': 18})


dat = Table.read('/home/mdocarm/Downloads/NEDLVS_20210922.fits')
df = dat.to_pandas()

# SELECTION
df = df[df.z<0.1]
df = df[df.Mstar>1e+11]
df = df[df.SFR_W4<16]
df.Lum_Ks = df.Lum_Ks*7.37/(3.839e33)
df = df[df.Lum_Ks>3e11]
df.dropna()

df = df.sample(850, random_state=1)

plt.scatter(np.log10(df.Mstar), df.SFR_W4)
plt.plot(11.95, 1.45, 'o', c='r')
plt.text(11.95, 0.6, 'UGC 2885', c='black', fontsize=12)
plt.ylabel('SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')
plt.clf()

# LINEAR REGRESSION
from scipy import stats

res = stats.linregress(np.log10(df.Mstar), df.SFR_W4)
# print(f"R-squared: {res.rvalue**2:.6f}")
# R-squared: 0.032475

plt.scatter(np.log10(df.Mstar), df.SFR_W4)
plt.plot(np.log10(df.Mstar), res.intercept + res.slope*np.log10(df.Mstar), 
         'r', label='fitted line')
plt.plot(11.95, 1.45, 'o', c='r')
plt.text(11.95, 0.6, 'UGC 2885', c='black', fontsize=12)
plt.ylabel('SFR ($M_{\odot} \; yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')
plt.clf()

df = df[['Mstar', 'SFR_W4']]

# CLUSTERING TECHNIQUES

from sklearn.mixture import BayesianGaussianMixture

gm = BayesianGaussianMixture(n_components=2, random_state=1, 
                              covariance_type='spherical', init_params='random').fit(df)
labels = gm.predict(df)
gm.means_

from sklearn.cluster import Birch

# birch = Birch(n_clusters=2).fit(df)
# labels = birch.predict(df)

from sklearn.cluster import KMeans

# kmeans = KMeans(n_clusters=2, init='k-means++')
# kmeans.fit(df)
# labels = kmeans.labels_

plt.scatter(np.log10(df.Mstar), np.log10(df.SFR_W4/df.Mstar), c=labels)
plt.plot(11.95, -11.79, 'o', c='r')
plt.text(11.95, -11.9, 'UGC 2885', c='black', fontsize=12)
plt.ylabel('log sSFR ($yr^{-1}$)')
plt.xlabel('log Stellar Mass ($M_{\odot}$)')
plt.show()

from sklearn.metrics import silhouette_score

silhouette_score = silhouette_score(df, labels, metric='euclidean')
print(silhouette_score)
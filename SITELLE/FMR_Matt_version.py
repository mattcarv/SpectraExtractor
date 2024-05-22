import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 18})

df = pd.read_csv('C:/Users/mathe/Downloads/fmr.csv')

df.set_index(df.columns[0], inplace=True)

df.index.name = None

print(df)

#%%
import seaborn as sns

plt.figure(figsize=(24, 12))
sns.heatmap(df.iloc[7:,:], cmap='viridis_r', annot=True, fmt=".2f", cbar=False)
plt.scatter(15, 5.5, marker='*', color='lightblue', s=700)
plt.annotate('', xy=(17.3, 6.3), xytext=(15, 6.3), 
             arrowprops=dict(color='orange', arrowstyle='-|>', linewidth=5))
plt.scatter(17.5, 6.3, marker='*', color='orange', s=700)
plt.xlabel('log Stellar Mass ($M_{\odot}$)')
plt.ylabel('log SFR ($M_{\odot} yr^{-1}$)')
plt.autoscale()
plt.show()


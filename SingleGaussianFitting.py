import os
import threading
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

results = {}
sigmas = {}

def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x-mu)**2 / (2*sigma**2))

def find_peak(file_path, noise_range, plot=False):
    try:
        distribution = np.loadtxt(file_path)

        x_axis = np.linspace(4383.3411648850003, 7733.3411648850003, 136)
        x = np.arange(len(distribution))

        noise_mask = (distribution >= noise_range[0]) & (distribution <= noise_range[1])
        distribution[noise_mask] = 0

        peak = np.argmax(distribution)
        mu, sigma = peak, len(distribution) // 10
        A = np.max(distribution)
        params, _ = curve_fit(gaussian, x, distribution, p0=[mu, sigma, A])
        area = np.sum(gaussian(x, *params))

        if plot:
            plt.plot(x_axis, distribution, 'bo', label='Original Distribution')
            plt.plot(x_axis, gaussian(x, *params), 'r', label='Fitted Gaussian')
            plt.xlabel('Velocity (Km/s)')
            plt.ylabel('Flux (K)')
            plt.legend()
            plt.show()
            
        # print("mu: ", params[0])
        # print("sigma: ", params[1])
        # print("A: ", params[2], 'K')
        # print("Integrated Flux: ", area, 'K Km/s')
        results[file_path] = area
        sigmas[file_path] = params[1]
        return params[0]
    except:
        pass
    
folder_path = 'C:/Users/mathe/OneDrive/Documents/PROJECTUGC2885-2022/CO files-20221207T192945Z-001/CO files/spectra10'
files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]

valid_files = []
for file in files:
    file_path = os.path.join(folder_path, file)
    try:
        data = np.loadtxt(file_path)
        if not np.isnan(data).any():
            valid_files.append(file)
    except:

        pass
    
data = np.array(valid_files)
# print(data)

specs = []
threads = []
for d in data:
    x = threading.Thread(target=find_peak, args=(d, (-0.03, 0.01), False,))
    threads.append(x)
    
for thread in threads:
    thread.start()
    thread.join()
    

print('End processing')


# for r in results:
#     print(f"{r}: {results[r]}")
    
df = pd.DataFrame({'files': results.keys(), 'values': results.values(), 'sigmas': sigmas.values()})
df.to_csv('testfluxes.csv')
print(df)

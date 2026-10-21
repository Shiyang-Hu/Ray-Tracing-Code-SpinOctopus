import numpy as np
import matplotlib.pyplot as plt 
import os
from astropy.timeseries import LombScargle
current_dir = os.path.dirname(os.path.abspath(__file__))
times = np.loadtxt("data/time-series.dat",unpack=True)
time_min = min(times)
time_max = max(times)
count_rate_max=10
#time_max = 10000
#time_min = 0
bins =400
range = (time_min,time_max)
counts, bin_edges = np.histogram(times, bins=bins,range=range)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_widths = bin_edges[1:] - bin_edges[:-1]
frequency, power = LombScargle(bin_centers, counts/bin_widths).autopower()
fig, axes = plt.subplots(2, 1, figsize=(15, 8))
axes[0].step(bin_centers,counts/bin_widths,where="mid")
#axes[0].step(bin_centers,counts,where="mid")
#axes[0].plot(bin_centers,counts)
#axes[0].set_xlim(0,time_max)
axes[0].set_xlim(0,time_max)
axes[0].set_ylim(0,count_rate_max)
axes[0].set_xlabel("Time [M]")
axes[0].set_ylabel("Photon counting rate")


axes[1].plot(frequency, power)
axes[1].set_xlim(0,0.1)
axes[1].set_ylim(0,1)
axes[1].set_xlabel("Frequency")
axes[1].set_ylabel("Power")
plt.savefig('results/light curve.png', dpi=600, bbox_inches='tight')
plt.show()
np.savetxt("data/counts rate.dat", np.column_stack((bin_centers, counts/bin_widths)),fmt="%.8f")
np.savetxt("data/power spectrum.dat", np.column_stack((frequency, power)),fmt="%.8f")
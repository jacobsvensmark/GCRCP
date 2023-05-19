import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from PLOT_utils import centergrid

lon_arr = np.arange(0,360,5)
lat_arr = np.arange(-90,95,5)

cutoff_map = np.genfromtxt(  '../dat/jupiter_I30_E1_28_06_2021_VERTICAL_rigidity.csv',delimiter=',')
cutoff_map = np.genfromtxt(  '../dat/present_B_VERTICAL_rigidity.csv',delimiter=',')
#cutoff_map = np.genfromtxt(  '../dat/jupiter_I30_E1_28_06_2021_VERTICAL_penumbra_info.csv',delimiter=',')

fig, axes = plt.subplots(nrows=1, ncols=1)
xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)
im = axes.pcolormesh(lon_arr,lat_arr,cutoff_map,shading='auto') #add vmin=xx, vmax=yy to set range
axes.set_title('Rigidity Map')
fig.subplots_adjust(right=0.8,wspace = 0.18 )
fig.colorbar(im)

plt.show()
#plt.savefig('Map.jpg')

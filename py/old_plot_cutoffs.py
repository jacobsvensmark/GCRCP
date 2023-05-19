import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from PLOT_utils import centergrid

lon_arr = np.arange(0,360,5)
lat_arr = np.arange(-90,95,5)

def get_data(filename):
    return np.genfromtxt(filename, delimiter=',')


levels = {'TOP':[1,2,3],'MID':[1,2,3,4,5,6,7,8],'LOW':[1,2,3,4,5,6,7,8,9,10,11,12]}

cutoff_maps = {}

for level in levels:
    for az in levels[level]:
       fname = '../dat/denoised/pres_B_'+level+str(az)+'_rigidity_denoised.csv'
       cutoff_maps[level+str(az)] = get_data(fname)



fig, axes = plt.subplots(nrows=5, ncols=5)
fig.set_size_inches(14.5, 8.5)
xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)

row = 0
col = 0
for level in levels:
    for az in levels[level]:
        im = axes[col,row].pcolormesh(lon_arr,lat_arr,cutoff_maps[level+str(az)],vmin=0.0,vmax=30.0)
        axes[col,row].set_title(level+str(az))
        row = row+1
        if row==5:
            row = 0
            col = col + 1

fig.subplots_adjust(right=0.8,wspace = 0.18 )
fig.colorbar(im, ax=axes.ravel().tolist())


plt.show()

#fig2, axes2 = plt.subplots(nrows=2, ncols=3)
#fig2.set_size_inches(14.5, 8.5)
#im2 = axes2[0,0].pcolormesh(lon_arr,lat_arr,N26,vmin=0.0, vmax=20)
#im2 = axes2[0,1].pcolormesh(lon_arr,lat_arr,E26,vmin=0.0, vmax=20)
#im2 = axes2[1,0].pcolormesh(lon_arr,lat_arr,S26,vmin=0.0, vmax=20)
#im2 = axes2[1,1].pcolormesh(lon_arr,lat_arr,W26,vmin=0.0, vmax=20)
#im2 = axes2[0,2].pcolormesh(lon_arr,lat_arr,vertical,vmin=0.0, vmax=20)
#axes2[0,0].set_title("N26")
#axes2[0,1].set_title("E26")
#axes2[1,0].set_title("S26")
#axes2[1,1].set_title("W26")
#axes2[0,2].set_title("vertical")
#fig2.subplots_adjust(right=0.8,wspace = 0.18 )
#fig2.colorbar(im2, ax=axes2.ravel().tolist())
plt.show()

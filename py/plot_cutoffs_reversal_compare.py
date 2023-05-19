import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from PLOT_utils import centergrid

lon_arr = np.arange(-182.5,187.5,5)
lat_arr = np.arange(-92.5,97.5,5)
xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)

def get_data(filename):
    return np.genfromtxt(filename, delimiter=',')

time_array=("794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437")
letter_arr = ['a)','b)','c)','d)','e)','f)','g)','h)','i)']
fig, axes = plt.subplots(ncols=3,nrows=3, sharex=True, sharey=True)
fig.set_size_inches(14.5, 8.5)

cnt = 0
cutoff_map = {}
for row in range(3):
    for col in range(3):
        time = time_array[cnt]
        fname = '../dat/reversal_T'+time+'_VERTICAL_rigidity.csv'
        cutoff_map['map'+time] = np.roll(get_data(fname),36,axis=1)
        cutoff_map['map'+time] = np.c_[cutoff_map['map'+time],cutoff_map['map'+time][:,0]]
        im = axes[row,col].pcolormesh(lon_arr,lat_arr,cutoff_map['map'+time],vmin=0.0,vmax=12.0,linewidth=0,rasterized=True)
        im.set_edgecolor('face')
        axes[row,col].set_title(letter_arr[cnt]+' T={:.1f}kyr'.format(float(time)), y=1.0, pad=-14,color='white')
        axes[row,col].set_xlim([-182.5,182.5])
        axes[row,col].set_ylim([-90,90])
        axes[row,col].set_xticks(np.arange(-180,225, 90))
        axes[row,col].set_yticks(np.arange(-90,120,30))
        cnt += 1
fig.subplots_adjust(bottom=0.1,top=0.9,left=0.1,right=0.8,hspace=0.08,wspace=0.1)

cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.set_ticks(np.arange(0, 12, 1.0))
cbar.set_label("Cutoff rigidity [GV]")



lon_arr = np.arange(-182.5,187.5,5)
lat_arr = np.arange(-92.5,97.5,5)
xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)

time_array=("794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437","794.247437")
letter_arr = ['a)','b)','c)','d)','e)','f)','g)','h)','i)']
fig, axes = plt.subplots(ncols=3,nrows=3, sharex=True, sharey=True)
fig.set_size_inches(14.5, 8.5)

cnt = 0
cutoff_map = {}
for row in range(3):
    for col in range(3):
        time = time_array[cnt]
        fname = '../dat/5x5_reversal/T'+time+'_rigidity.csv'
        #cutoff_map['map'+time] = np.roll(get_data(fname),36,axis=1)
        cutoff_map['map'+time] = get_data(fname)
        cutoff_map['map'+time] = np.c_[cutoff_map['map'+time],cutoff_map['map'+time][:,0]]
        im = axes[row,col].pcolormesh(lon_arr,lat_arr,cutoff_map['map'+time],vmin=0.0,vmax=12.0,linewidth=0,rasterized=True)
        im.set_edgecolor('face')
        axes[row,col].set_title(letter_arr[cnt]+' T={:.1f}kyr'.format(float(time)), y=1.0, pad=-14,color='white')
        axes[row,col].set_xlim([-182.5,182.5])
        axes[row,col].set_ylim([-90,90])
        axes[row,col].set_xticks(np.arange(-180,225, 90))
        axes[row,col].set_yticks(np.arange(-90,120,30))
        cnt += 1
fig.subplots_adjust(bottom=0.1,top=0.9,left=0.1,right=0.8,hspace=0.08,wspace=0.1)

cb_ax = fig.add_axes([0.93, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
cbar.set_ticks(np.arange(0, 12, 1.0))
cbar.set_label("Cutoff rigidity [GV]")

plt.show()

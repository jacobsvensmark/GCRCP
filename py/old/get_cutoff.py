import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from PLOT_utils import centergrid

nlon=72
nlat=37
lat_lower = -85
lat_higher = 85
delta = 5
# 1 lvl
#directions = ["VERTICAL"]
# 2 lvl
#direction = "vertical" # 'vertical', 'N26', 'E26', 'S26' or 'W26'
# 3 lvl
directions = ["TOP1","TOP2","TOP3",
              "MID1","MID2","MID3","MID4","MID5","MID6","MID7","MID8",
              "LOW1","LOW2","LOW3","LOW4","LOW5","LOW6","LOW7","LOW8","LOW9","LOW10","LOW11","LOW12"]
lon_arr = np.arange(0,360,delta)
lat_arr = np.arange(-90,95,delta)

# For the fourth array, the numbers mean the following:
# 1: No unallowed paths found
# 2: 
for direction in directions:
    file_prefix = "pres_B_" + direction
    cutoff = np.empty([4,37,72])
    cutoff[:] = np.nan
    for lat_idx,lat in enumerate(lat_arr):
        print(lat)
        filename = "../dat/5x5grid_3lvl/"+file_prefix +'_lat_'+ str(lat) + '.txt'
        try:
            df = pd.read_table(filename,header=None, delimiter=r"\s+")
        except:
            print("Rigidity file for latitude {} not found at {}".format(lat,filename))
            cutoff[:,lat_idx,:] = np.nan
            continue
        df.columns = ['lat','lon','E','az_or_zen','zen_or_az','unknown1','unknown2','unknown3','unknown4','allowed','unknown5']
        for lon_idx,lon in enumerate(lon_arr):
            df_cut = df[df["lon"]==lon]
            a = df_cut['allowed'].values
            b = a.astype(int)
            if len(b) == 0: # Data point missing
                cutoff[0,lat_idx,lon_idx] = np.nan
                cutoff[1,lat_idx,lon_idx] = np.nan
                cutoff[2,lat_idx,lon_idx] = np.nan
                cutoff[3,lat_idx,lon_idx] = 5 
            elif np.count_nonzero(b) == 0: # Penumbra not reached at lowest energy
                cutoff[0,lat_idx,lon_idx] = np.nan
                cutoff[1,lat_idx,lon_idx] = np.nan
                cutoff[2,lat_idx,lon_idx] = df_cut['E'].values[-1]
                cutoff[3,lat_idx,lon_idx] = 1
            elif np.count_nonzero(b) == len(b): # Penumbra starts before highest energy sampeled
                cutoff[0,lat_idx,lon_idx] = np.nan
                cutoff[1,lat_idx,lon_idx] = np.nan
                cutoff[2,lat_idx,lon_idx] = df_cut['E'].values[0]
                cutoff[3,lat_idx,lon_idx] = 4
            else:
                b[b>1] = 1
                idx1 = np.nonzero(b==1)[0][0] 
                idx2 = np.nonzero(b==0)[0][-1]
                if (idx2 == len(b)-1): # Penumbra stretches into lowest energy
                   cutoff[0,lat_idx,lon_idx] = df_cut['E'].values[idx1]
                   cutoff[1,lat_idx,lon_idx] = df_cut['E'].values[idx2]
                   cutoff[2,lat_idx,lon_idx] = (cutoff[0,lat_idx,lon_idx]+cutoff[1,lat_idx,lon_idx])/2.0
                   cutoff[3,lat_idx,lon_idx] = 2
                elif (idx1-1 == idx2): # Penumbra has 0 width
                   cutoff[0,lat_idx,lon_idx] = df_cut['E'].values[idx1]
                   cutoff[1,lat_idx,lon_idx] = df_cut['E'].values[idx2]
                   cutoff[2,lat_idx,lon_idx] = (cutoff[0,lat_idx,lon_idx]+cutoff[1,lat_idx,lon_idx])/2.0
                   cutoff[3,lat_idx,lon_idx] = 3
                else:
                   cutoff[0,lat_idx,lon_idx] = df_cut['E'].values[idx1]
                   cutoff[1,lat_idx,lon_idx] = df_cut['E'].values[idx2]
                   cutoff[2,lat_idx,lon_idx] = (cutoff[0,lat_idx,lon_idx]+cutoff[1,lat_idx,lon_idx])/2.0
                   cutoff[3,lat_idx,lon_idx] = 0
    
    xlist = centergrid(lon_arr)
    ylist = centergrid(lat_arr)
    #plt.pcolormesh(lon_arr,lat_arr,cutoff[2,:,:])
    #plt.colorbar(label='Energy in GeV?? Or is it Rigidity in GV?')
    #plt.show()
    
    import csv
    
    with open("../dat/"+file_prefix+"_rigidity.csv","w+") as my_csv:
        csvWriter = csv.writer(my_csv,delimiter=',')
        csvWriter.writerows(cutoff[2,:,:])
    
    with open("../dat/"+file_prefix+"_penumbra_info.csv","w+") as my_csv:
        csvWriter = csv.writer(my_csv,delimiter=',')
        csvWriter.writerows(cutoff[3,:,:])

print("done") 

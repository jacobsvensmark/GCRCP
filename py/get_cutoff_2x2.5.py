iport numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from PLOT_utils import centergrid

magneticfield="present"

lat_array_str=(  "-89.50"  ,"-88.00", "-86.00", "-84.00", "-82.00", "-80.00", "-78.00", "-76.00",
           "-74.00"  ,"-72.00", "-70.00", "-68.00", "-66.00", "-64.00", "-62.00", "-60.00", 
           "-58.00"  ,"-56.00", "-54.00", "-52.00", "-50.00", "-48.00", "-46.00", "-44.00", 
           "-42.00"  ,"-40.00", "-38.00", "-36.00", "-34.00", "-32.00", "-30.00", "-28.00",
           "-26.00"  ,"-24.00", "-22.00", "-20.00", "-18.00", "-16.00", "-14.00", "-12.00",
           "-10.00"  ,"-8.00",  "-6.00",  "-4.00",  "-2.00",  "0.00",   "2.00",   "4.00",
           "6.00"    ,"8.00",   "10.00",  "12.00",  "14.00",  "16.00",  "18.00",  "20.00",
           "22.00"   ,"24.00",  "26.00",  "28.00",  "30.00",  "32.00",  "34.00",  "36.00",
           "38.00"   ,"40.00",  "42.00",  "44.00",  "46.00",  "48.00",  "50.00",  "52.00",
           "54.00"   ,"56.00",  "58.00",  "60.00",  "62.00",  "64.00",  "66.00",  "68.00",
           "70.00"   ,"72.00",  "74.00",  "76.00",  "78.00",  "80.00",  "82.00",  "84.00",
           "86.00"   ,"88.00",  "89.50")

lon_array_str=("-180.00","-177.50","-175.00","-172.50","-170.00","-167.50","-165.00","-162.50",
               "-160.00","-157.50","-155.00","-152.50","-150.00","-147.50","-145.00","-142.50",
               "-140.00","-137.50","-135.00","-132.50","-130.00","-127.50","-125.00","-122.50",
               "-120.00","-117.50","-115.00","-112.50","-110.00","-107.50","-105.00","-102.50",
               "-100.00","-97.50","-95.00","-92.50","-90.00","-87.50","-85.00","-82.50",
               "-80.00","-77.50","-75.00","-72.50","-70.00","-67.50","-65.00","-62.50",
               "-60.00","-57.50","-55.00","-52.50","-50.00","-47.50","-45.00","-42.50",
               "-40.00","-37.50","-35.00","-32.50","-30.00","-27.50","-25.00","-22.50",
               "-20.00","-17.50","-15.00","-12.50","-10.00","-7.50","-5.00","-2.50",
               "0.00","2.50","5.00","7.50","10.00","12.50","15.00","17.50",
               "20.00","22.50","25.00","27.50","30.00","32.50","35.00","37.50",
               "40.00","42.50","45.00","47.50","50.00","52.50","55.00","57.50",
               "60.00","62.50","65.00","67.50","70.00","72.50","75.00","77.50",
               "80.00","82.50","85.00","87.50","90.00","92.50","95.00","97.50",
               "100.00","102.50","105.00","107.50","110.00","112.50","115.00","117.50",
               "120.00","122.50","125.00","127.50","130.00","132.50","135.00","137.50",
               "140.00","142.50","145.00","147.50","150.00","152.50","155.00","157.50",
               "160.00","162.50","165.00","167.50","170.00","172.50","175.00","177.50")

cutoff = np.empty([4,nlat,nlon])
cutoff[:] = np.nan
for lat_idx,lat in enumerate(lat_array_str):
    print(lat)
    for lon_idx,lon in enumerate(lon_array_str):
        filename = "../dat/2x2.5_grid_model_output/"+magneticfield +'_VERTICAL_lat_'+ lat + '_lon_'+lon+'.txt'
        try:
            df = pd.read_table(filename,header=None, delimiter=r"\s+")
        except:
            print("Rigidity file for latitude {} not found at {}".format(lat,filename))
            cutoff[:,lat_idx,:] = np.nan
            continue
        df.columns = ['lat','lon','E','az_or_zen','zen_or_az','unknown1','unknown2','unknown3','unknown4','allowed','unknown5']
        df_cut = df[df["lon"]==float(lon)]
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

import csv

with open("../dat/"+magneticfield+"_rigidity_2x2.5.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(cutoff[2,:,:])

with open("../dat/"+magneticfield+"_penumbra_info_2x2.5.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(cutoff[3,:,:])

print("done") 

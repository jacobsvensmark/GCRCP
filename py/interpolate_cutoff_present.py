import numpy as np
import matplotlib.pyplot as plt
from PLOT_utils import centergrid
from scipy.interpolate import interp2d,RectBivariateSpline,SmoothBivariateSpline,griddata
import csv

lon_arr = np.array([ -180, -177.5, -175, -172.5, -170, -167.5, -165, -162.5, -160, -157.5,
                     -155, -152.5, -150, -147.5, -145, -142.5, -140, -137.5, -135, -132.5,
                     -130, -127.5, -125, -122.5, -120, -117.5, -115, -112.5, -110, -107.5,
                     -105, -102.5, -100, -97.5, -95, -92.5, -90, -87.5, -85, -82.5, -80,
                     -77.5, -75, -72.5, -70, -67.5, -65, -62.5, -60, -57.5, -55, -52.5, -50,
                     -47.5, -45, -42.5, -40, -37.5, -35, -32.5, -30, -27.5, -25, -22.5, -20,
                     -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15,
                     17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50,
                     52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5, 75, 77.5, 80, 82.5, 85,
                     87.5, 90, 92.5, 95, 97.5, 100, 102.5, 105, 107.5, 110, 112.5, 115, 117.5,
                     120, 122.5, 125, 127.5, 130, 132.5, 135, 137.5, 140, 142.5, 145, 147.5,
                     150, 152.5, 155, 157.5, 160, 162.5, 165, 167.5, 170, 172.5, 175, 177.5 ])

lat_arr = np.array([ -89.5, -88, -86, -84, -82, -80, -78, -76, -74, -72, -70, -68, -66,
                      -64, -62, -60, -58, -56, -54, -52, -50, -48, -46, -44, -42, -40, -38,
                      -36, -34, -32, -30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10, -8,
                      -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
                      32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66,
                      68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 89.5])


def get_data(filename):
    return np.genfromtxt(filename, delimiter=',')

dd = {}
dd['present'] = {'lon':[-180 ,-140,-180 ,-156,-170 ,-160,145,180],'lat':[-87,-87,-85,-85,-83,-83,-87,-81],'flag':[2,2,2,2]} 
levels = {'present':['']}

cutoff_maps = {}
cutoff_maps['present'] = get_data('../dat/present_VERTICAL_rigidity_2x2.5_individ.csv')
xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)

row = 0
col = 0
for level in levels:
    for az in levels[level]:
        # PLOT STUFF
        fig,axes = plt.subplots(nrows=2,ncols=3)
        fig.set_size_inches(14.5, 8.5)
        cmap = cutoff_maps[level+str(az)]
        im = axes[0,0].pcolormesh(lon_arr,lat_arr,np.log10(cmap),vmin=np.log10(0.06),vmax=np.log10(30.0))
        axes[0,0].set_title(level+str(az))
        for i in range(len(lon_arr)):
            axes[0,1].plot(lat_arr,cmap[:,i]) 
        axes[0,1].set_xlabel('Latitude')
        for i in range(len(lat_arr)):
            axes[0,2].plot(lon_arr,cmap[i,:]) 
        axes[0,2].set_xlabel('Longitude')
#       plt.show()
        x_dd = [(x+180)/2.5 for x in dd[level+str(az)]['lon'] ]
        x_dd[:] = [int(x) for x in x_dd]
        y_dd = [(y+90)/2.0 for y in dd[level+str(az)]['lat'] ]
        y_dd[:] = [int(x) for x in y_dd]
        f_dd = dd[level+str(az)]['flag']
        for idx,flag in enumerate(f_dd):
            if flag == 1: # Set noise to lowest possible value
                cmap[y_dd[2*idx]:y_dd[2*idx+1]+1,x_dd[2*idx]:x_dd[2*idx+1]+1] = 0.06
            if flag == 2:
                cmap[y_dd[2*idx]:y_dd[2*idx+1]+1,x_dd[2*idx]:x_dd[2*idx+1]+1] = np.nan
        if np.isnan(cmap).any():
            Na=len(lat_arr)
            No=len(lon_arr)
            lon_large = np.arange(-540,540,2.5)
            lat_large = np.arange(-270,272,2)
            lat_large = np.append(np.append(lat_arr-180,lat_arr),lat_arr+180)
            cmap_large = np.zeros([3*Na,3*No])
            cmap_large[Na:2*Na,     0:No   ] = cmap
            cmap_large[Na:2*Na,     No:2*No ] = cmap
            cmap_large[Na:2*Na,     2*No:3*No ] = cmap
            cmap_large[0:Na-1,         int(0.5*No):int(1.5*No)   ] = np.flipud(cmap[1:,:])
            cmap_large[0:Na-1,         int(1.5*No):int(2.5*No) ] = np.flipud(cmap[1:,:])
            cmap_large[0:Na-1,         0:int(0.5*No)     ] = np.flipud(cmap[1:,int(0.5*No):No])
            cmap_large[0:Na-1,         int(2.5*No):3*No ] = np.flipud(cmap[1:,0:int(0.5*No)])
            cmap_large[2*Na:3*Na,int(0.5*No):int(1.5*No)   ] = np.flipud(cmap[:,:])
            cmap_large[2*Na:3*Na, int(1.5*No):int(2.5*No) ] = np.flipud(cmap[:,:])
            cmap_large[2*Na:3*Na, 0:int(0.5*No)     ] = np.flipud(cmap[:,int(0.5*No):No])
            cmap_large[2*Na:3*Na, int(2.5*No):3*No ] = np.flipud(cmap[:,0:int(0.5*No)])

            cmap_large = np.ma.masked_invalid(cmap_large)
            lonv,latv = np.meshgrid(lon_large,lat_large ) # THIS IS WHERE WE ARE AT!!!
            lat1 = latv[~cmap_large.mask]
            lon1 = lonv[~cmap_large.mask]
            new_cmap = cmap_large[~cmap_large.mask]
            lon_mesh, lat_mesh = np.meshgrid(lon_arr, lat_arr)
            cmap2 = griddata((lon1, lat1), new_cmap.ravel(),(lon_mesh, lat_mesh), method='cubic')
            cmap = np.reshape(cmap2,np.shape(cmap))
            cmap[cmap<0.06] = 0.06
            n0 = np.nansum(cmap < 0)
            print(n0)
            if (n0>0):
                print("ERROR: Negative rigidity values in interpolation")
            cmap[0,:]  = np.nan 
            cmap[-1,:] = np.nan 
            
        fig2,axes2 = plt.subplots()
        axes2.pcolormesh(lon_large,lat_large,np.log10(cmap_large),vmin=np.log10(0.06),vmax=np.log10(30.0),shading='auto')
        plt.plot()

        fig.set_size_inches(14.5, 8.5)
        im = axes[1,0].pcolormesh(lon_arr,lat_arr,np.log10(cmap),vmin=np.log10(0.06),vmax=np.log10(30.0),shading='auto')
        axes[1,0].set_title(level+str(az)+" denoised")
        for i in range(len(lon_arr)):
            axes[1,1].plot(lat_arr,cmap[:,i]) 
        axes[1,1].set_xlabel('Latitude')
        for i in range(len(lat_arr)):
            axes[1,2].plot(lon_arr,cmap[i,:]) 
        axes[1,2].set_xlabel('Longitude')
        plt.show()
         
        outname = '../dat/denoised/pres_B_'+level+str(az)+'_rigidity_denoised.csv'
        with open(outname,"w+") as my_csv:
            csvWriter = csv.writer(my_csv,delimiter=',')
            csvWriter.writerows(cmap)



import numpy as np
import matplotlib.pyplot as plt
from PLOT_utils import centergrid
from scipy.interpolate import interp2d,RectBivariateSpline,SmoothBivariateSpline,griddata
import csv

lon_arr = np.arange(0,360,5)
lat_arr = np.arange(-90,95,5)

def get_data(filename):
    return np.genfromtxt(filename, delimiter=',')
dd = {}
dd['TOP1'] = {'lon':[165,185],'lat':[-85,-85],'flag':[1]}
dd['TOP2'] = {'lon':[145,150],'lat':[-85,-85],'flag':[1]}
dd['TOP3'] = {'lon':[295,300,150,185,185,270,185,210],'lat':[85,85,-85,-80,-85,-85,-80,-80],'flag':[1,1,2,2]}
dd['MID1'] = {'lon':[145,225,160,165,175,180,195,205,215,220,155,160],'lat':[-85,-85,-85,-85,-85,-85,-85,-85,-85,-85,-75,-75],'flag':[2,2,2,2,2,2]}
dd['MID2'] = {'lon':[145,185],'lat':[-85,-85],'flag':[2]}
dd['MID3'] = {'lon':[150,170,155,160],'lat':[-85,-85,-80,-80],'flag':[2,2]}
dd['MID4'] = {'lon':[130,155,280,285],'lat':[-85,-80,85,85],'flag':[1,1]}
dd['MID5'] = {'lon':[0,360,85,200,140,145,45,50,285,335],'lat':[-85,-85,-80,-80,-75,-75,85,85,85,85],'flag':[2,2,2,2,2]}
dd['MID6'] = {'lon':[140,355,5,15,0,150,180,355,0,5,20,25,40,45,180,185,330,335],'lat':[-85,-70,-85,-85,85,85,85,85,80,80,75,75,80,80,80,80,80,80],'flag':[2,2,2,2,2,2,2,2,2]}
dd['MID7'] = {'lon':[180,180,140,300,310,315,170,260,180,200,215,225,230,240,150,150,215,215,0,355,5,5,25,25,35,35,55,55,70,70,100,100,300,300,315,340],'lat':[75,75,-85,-85,-85,-85,-80,-80,-75,-75,-75,-75,-75,-75,-65,-65,-65,-65,80,85,75,75,75,75,75,75,75,75,75,75,75,75,70,75,75,75],'flag':[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}
dd['MID8'] = {'lon':[0,355,155,260,160,165,205,210,245,245],'lat':[85,85,-85,-85,-80,-80,-80,-80,-80,-80],'flag':[2,2,2,2,2]}
dd['LOW1'] = {'lon':[150,245,160,160,190,190,145,145,260,260],'lat':[-85,-85,-80,-80,-80,-80,-75,-75,85,85],'flag':[2,2,2,2,2]}
dd['LOW2'] = {'lon':[150,215],'lat':[-85,-85],'flag':[2]}
dd['LOW3'] = {'lon':[125,200],'lat':[-85,-80],'flag':[1]}
dd['LOW4'] = {'lon':[130,175],'lat':[-85,-85],'flag':[1]}
dd['LOW5'] = {'lon':[115,170,145,150],'lat':[-85,-85,-80,-80],'flag':[1,1]}
dd['LOW6'] = {'lon':[100,105,125,155,130,135,350,355],'lat':[-85,-85,-85,-85,-75,-75,85,85],'flag':[1,1,1,1]}
dd['LOW7'] = {'lon':[  0,355,145,150,295,355, 0,15],
              'lat':[-85,-85,-80,-80, 80, 85,85,85],
             'flag':[1,1,2,2]}
dd['LOW8'] = {'lon':[ 0,355,25,30,  0,355,   0,115,125,355,   0, 15, 25, 40, 55,355,  0, 20, 40, 75, 90,320,330,355, 315,350,305,310,320,325],
              'lat':[80, 85,75,75,-85,-85, -80,-80,-80,-80, -75,-75,-75,-75,-75,-75,-70,-70,-70,-70,-70,-70,-70,-70, -65,-65, 70, 70, 70, 70],
             'flag':[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}
dd['LOW9'] = {'lon':[0, 355,255,265,285,350,20,50,150,155,195,205,230,235,165,175,100,135, 0, 5,40,85,150,155,290,305,350,355,280,285,345,350,105,110,105,310,315,355,  0, 50, 60, 95,  0, 25, 50, 55],
              'lat':[80, 85, 75, 75, 75, 75,75,75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75,75,75,65,65, 65, 65, 70, 70, 70, 70, 60, 60, 60, 60, 70, 70,-85,-60,-85,-75,-85,-85,-85,-85,-80,-70,-80,-80],
             'flag':[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}
dd['LOW10'] = {'lon':[ 0,355,315,320,110,355,125,130,145,320,140,195,210,250,260,305,140,165,195,210,225,250,270,275,300,305,225,230,280,285,260,265,140,145],
               'lat':[65, 85, 60, 60,-85,-85,-80,-80,-80,-80,-75,-75,-75,-75,-75,-75,-70,-70,-70,-70,-70,-70,-70,-70,-70,-70,-65,-65,-65,-65,-60,-60,-55,-55],'flag':[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}
dd['LOW11'] = {'lon':[ 0,355, 0,140,150,355, 0,25,35,40,45,65,70,175,185,195,205,220,230,235,245,355, 130,305,135,270,195,250,195,200,250,255,220,225],
               'lat':[80, 85,75, 75, 75, 75,70,70,70,70,70,70,70, 70, 70, 70, 70, 70, 70, 70, 70, 70,-85, -85,-80,-80,-75,-75,-70,-70,-70,-70,-65,-65],'flag':[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}
dd['LOW12'] = {'lon':[ 0,355,115,240,250,280,155,170,185,215,230,235,255,260,230,235,150,155,190,195],
               'lat':[85, 85,-85,-85,-85,-85,-80,-80,-80,-80,-80,-80,-80,-80,-75,-75,-70,-70,-70,-70],'flag':[2,2,2,2,2,2,2,2,2,2]}
levels = {'TOP':[1,2,3],'MID':[1,2,3,4,5,6,7,8],'LOW':[1,2,3,4,5,6,7,8,9,10,11,12]}

cutoff_maps = {}

for level in levels:
    for az in levels[level]:
       fname = '../dat/pres_B_'+level+str(az)+'_rigidity.csv'
       cutoff_maps[level+str(az)] = get_data(fname)

xlist = centergrid(lon_arr)
ylist = centergrid(lat_arr)

row = 0
col = 0
for level in levels:
    for az in levels[level]:
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
        x_dd = [(x)/5 for x in dd[level+str(az)]['lon'] ]
        x_dd[:] = [int(x) for x in x_dd]
        y_dd = [(y+90)/5 for y in dd[level+str(az)]['lat'] ]
        y_dd[:] = [int(x) for x in y_dd]
        f_dd = dd[level+str(az)]['flag']
        for idx,flag in enumerate(f_dd):
            if flag == 1: # Set noise to lowest possible value
                cmap[y_dd[2*idx]:y_dd[2*idx+1]+1,x_dd[2*idx]:x_dd[2*idx+1]+1] = 0.06
            if flag == 2:
                cmap[y_dd[2*idx]:y_dd[2*idx+1]+1,x_dd[2*idx]:x_dd[2*idx+1]+1] = np.nan
        if np.isnan(cmap).any():
            N=36
            lon_large = np.arange(-360,720,5)
            lat_large = np.arange(-270,275,5)
            cmap_large = np.zeros([3*N+1,3*2*N])
            cmap_large[N:2*N+1,     0:2*N   ] = cmap
            cmap_large[N:2*N+1,     2*N:4*N ] = cmap
            cmap_large[N:2*N+1,     4*N:6*N ] = cmap
            cmap_large[0:N,         N:3*N   ] = np.flipud(cmap[1:,:])
            cmap_large[0:N,         3*N:5*N ] = np.flipud(cmap[1:,:])
            cmap_large[0:N,         0:N     ] = np.flipud(cmap[1:,N:2*N])
            cmap_large[0:N,         5*N:6*N ] = np.flipud(cmap[1:,0:N])
            cmap_large[2*N+1:3*N+1, N:3*N   ] = np.flipud(cmap[:-1,:])
            cmap_large[2*N+1:3*N+1, 3*N:5*N ] = np.flipud(cmap[:-1,:])
            cmap_large[2*N+1:3*N+1, 0:N     ] = np.flipud(cmap[:-1,N:2*N])
            cmap_large[2*N+1:3*N+1, 5*N:6*N ] = np.flipud(cmap[:-1,0:N])

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
            
    #    fig2,axes2 = plt.subplots()
    #    axes2.pcolormesh(lon_large,lat_large,np.log10(cmap_large),vmin=np.log10(0.06),vmax=np.log10(30.0),shading='auto')
    #    plt.plot()

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



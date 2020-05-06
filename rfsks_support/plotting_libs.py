from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)
import math

DEG2KM = 111.2

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        print("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
 
    #m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)

def plot_topo(map,cmap=plt.cm.jet):
    #20 minute bathymetry/topography data
    etopo = np.loadtxt('topo/etopo20data.gz')
    lons  = np.loadtxt('topo/etopo20lons.gz')
    lats  = np.loadtxt('topo/etopo20lats.gz')
    # shift data so lons go from -180 to 180 instead of 20 to 380.
    etopo,lons = shiftgrid(180.,etopo,lons,start=False)
    lons, lats = np.meshgrid(lons, lats)
    # lons, lats = map(*np.meshgrid(lons,lats))
    cs = map.pcolormesh(lons,lats,etopo,cmap=cmap,latlon=True,shading='gouraud')



def plot_events_loc(eq_map,evlons,evlats, evmgs, evdps,background=True):
    ## Plotting the earthquake map
    # print("------> Plotting the event location map")
    min_marker_size = 0.5
    dp100,dp300,dpelse=[],[],[]
    lt100,lt300,ltelse=[],[],[]
    lo100,lo300,loelse=[],[],[]
    mg100,mg300,mgelse=[],[],[]
    for lon, lat, mag, dp in zip(evlons, evlats, evmgs, evdps):
        if dp < 100.0:
            x,y = eq_map(lon, lat)
            msize = min_marker_size * mag **3
            dp100.append(dp)
            lt100.append(y)
            lo100.append(x)
            mg100.append(msize)
            
        elif dp < 300.0:
            x,y = eq_map(lon, lat)
            msize = min_marker_size * mag **3
            dp300.append(dp)
            lt300.append(y)
            lo300.append(x)
            mg300.append(msize)
            
        else:
            x,y = eq_map(lon, lat)
            msize = min_marker_size * mag **3
            dpelse.append(dp)
            ltelse.append(y)
            loelse.append(x)
            mgelse.append(msize)
    if background:
        eq_map.scatter(lo100, lt100, facecolors='none', marker='o', s=mg100,edgecolors='k',linewidths=0.5, zorder=10)
        eq_map.scatter(lo300, lt300, facecolors='none', marker='o', s=mg300,edgecolors='k',linewidths=0.5, zorder=10)
        eq_map.scatter(loelse, ltelse, facecolors='none', marker='o', s=mgelse,edgecolors='k',linewidths=0.5, zorder=10)
    else:
        eq_map.scatter(lo100, lt100, c='g', marker='o', s=mg100,edgecolors='k',linewidths=0.1, zorder=10)
        eq_map.scatter(lo300, lt300, c='b', marker='o', s=mg300,edgecolors='k',linewidths=0.1, zorder=10)
        eq_map.scatter(loelse, ltelse, c='r', marker='o', s=mgelse,edgecolors='k',linewidths=0.1, zorder=10)


def plot_merc(resolution,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,topo=True):
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
    plt.figure(figsize=(8,8))
    map = Basemap(projection='merc',resolution = resolution, area_thresh = 1000., llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    
    if topo:
        plot_topo(map,cmap=plt.cm.rainbow)

    map.drawcoastlines(color='k',linewidth=0.5)
    # map.fillcontinents()
    map.drawcountries(color='k',linewidth=0.1)
    map.drawstates(color='gray',linewidth=0.05)
    map.drawrivers(color='blue',linewidth=0.05)
    map.drawmapboundary()
    map.drawparallels(np.linspace(llcrnrlat,urcrnrlat,5,dtype='int16').tolist(),labels=[1,0,0,0],linewidth=0)
    map.drawmeridians(np.linspace(llcrnrlon,urcrnrlon,5,dtype='int16').tolist(),labels=[0,0,0,1],linewidth=0)
    return map

def station_map(map, stns_lon, stns_lat,stns_name,figname="", destination="./",figfrmt='png'):
    networkset = set([val.split("_")[0] for val in stns_name])
    stations = [val.split("_")[1] for val in stns_name]
    color=iter(plt.cm.viridis(np.linspace(0,1,len(networkset))))
    netcolor = {net:c for net, c in zip(networkset, color)}
    for lon, lat, sta, net in zip(stns_lon, stns_lat, stations, [val.split("_")[0] for val in stns_name]):
        x,y = map(lon, lat)
        plt.text(x+20000, y, sta, ha="center", va="center", bbox=dict(boxstyle="round",ec=(0., 0., 0.),fc=(1., 1., 1.)),fontsize=5)
        map.plot(x, y,'^',color=netcolor[net], markersize=7,markeredgecolor='k',linewidth=0.1)
    for net in networkset:
        map.plot(np.NaN,np.NaN,'^',color=netcolor[net],label=net, markersize=7,markeredgecolor='k',linewidth=0.1)
    plt.legend(loc=4,fontsize=8)
    plt.savefig(destination+figname+'_map.'+figfrmt,dpi=200,bbox_inches='tight')
    plt.close('all')

def plot_point_on_basemap(map, point, angle, length):
    '''
    point - Tuple (x, y)
    angle - Angle in degrees.
    length - Length of the line to plot.
    '''

    # unpack the point
    x, y = point

    # find the start and end point
    halfleny = length/2 * math.sin(math.radians(float(angle)))
    halflenx = length/2 * math.cos(math.radians(float(angle)))

    endx,endy = map(x+halflenx,y+halfleny)
    startx,starty = map(x-halflenx,y-halfleny)
    map.plot([startx,endx],[starty,endy],color='k',zorder=3)

from cmath import rect, phase
from math import radians, degrees
def mean_angle(deg):
    return degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))


## Fine tuning of SKS
import yaml
with open('Settings/advSKSparam.yaml') as f:
    inpSKSdict = yaml.load(f, Loader=yaml.FullLoader)

plot_params_lev = inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['meas_seg_points']
lev1, lev2, lev3 = int(plot_params_lev['lev1']), int(plot_params_lev['lev2']), int(plot_params_lev['lev3'])
from rfsks_support.rfsks_extras import segregate_measurements
def plot_sks_station_map(sks_meas_all,figname):
    if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements']):
        station_data_0, station_data_14, station_data_4_11, station_data_15 = segregate_measurements(sks_meas_all)
                    
    lonmin, lonmax = sks_meas_all['LON'].min()-1, sks_meas_all['LON'].max()+1
    latmin, latmax = sks_meas_all['LAT'].min()-1, sks_meas_all['LAT'].max()+1
    if np.abs(lonmax-lonmin)<10 or np.abs(latmax-latmin)<2:
        lonmin, lonmax = lonmin - 2, lonmax + 2
        latmin, latmax = latmin - 2, latmax + 2
    else:
        lonmin, lonmax = lonmin - 0.5, lonmax + 0.5
        latmin, latmax = latmin - 0.5, latmax + 0.5

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    map = Basemap(projection='merc',resolution = 'h', area_thresh = 1000., llcrnrlon=lonmin, llcrnrlat=latmin,urcrnrlon=lonmax, urcrnrlat=latmax)
    
    # plot_topo(map,cmap=plt.cm.rainbow)
    map.etopo(zorder=2)
    # map.etopo(scale=2.5, alpha=0.5, zorder=2)

    map.drawcoastlines(color='k',linewidth=0.5)
    map.drawcountries(color='k',linewidth=0.1)
    map.drawstates(color='gray',linewidth=0.05)
    map.drawrivers(color='blue',linewidth=0.05)
    map.drawmapboundary()
    parallelmin = int(latmin)
    parallelmax = int(latmax)+1
    if np.abs(parallelmax - parallelmin)<5:
        parallelmax += 2
        parallelmin -= 2

    meridianmin = int(lonmin)
    meridianmax = int(lonmax)+1
    if np.abs(meridianmax - meridianmin)<5:
        meridianmax += 2
        meridianmin -= 2

    map.drawparallels(np.arange(parallelmin, parallelmax,dtype='int16').tolist(),labels=[1,0,0,0],linewidth=0)
    map.drawmeridians(np.arange(meridianmin, meridianmax,dtype='int16').tolist(),labels=[0,0,0,1],linewidth=0)
    map.drawmapboundary(color='k', linewidth=2, zorder=1)
    legendarray = []
    for a in [1, 2, 3]:
        legendarray.append(map.scatter([], [], c='b', alpha=0.6, s=60*a,label=f"{a}s",edgecolors='k'))



    if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements']):
        ## no measurements
        if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_no_measurement']):
            stlon0s,stlat0s = map(station_data_0['LON'].values,station_data_0['LAT'].values)
            map.scatter(stlon0s, stlat0s, c='red', marker='o', s=60, edgecolors='k',linewidths=0.3, zorder=4)
            legendarray.append(map.scatter([], [], c='r', alpha=0.99, s=60, edgecolors='k'))
        ## plot null measurements
        if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_null_measurements']):
            station_data_null = sks_meas_all.loc[(sks_meas_all['NumNull']>0) & (sks_meas_all['NumMeasurements']==0)]
            stlonnull,stlatnull = map(station_data_null['LON'].values,station_data_null['LAT'].values)
            map.scatter(stlonnull, stlatnull, c='white', marker='o', s=60, edgecolors='k',linewidths=0.3, zorder=4)
            legendarray.append(map.scatter([], [], c='white', alpha=0.99, s=60, edgecolors='k'))

        stlon1s,stlat1s = map(station_data_14['LON'].values,station_data_14['LAT'].values)
        stlon4s,stlat4s = map(station_data_4_11['LON'].values,station_data_4_11['LAT'].values)
        stlon5s,stlat5s = map(station_data_15['LON'].values,station_data_15['LAT'].values)
        
        map.scatter(stlon1s, stlat1s, c='cornflowerblue', marker='o', s=60*station_data_14['AvgLagTime'],edgecolors='k',linewidths=0.1, zorder=4)
        map.scatter(stlon4s, stlat4s, c='navy', marker='o', s=60*station_data_4_11['AvgLagTime'],edgecolors='k',linewidths=0.1, zorder=4)
        map.scatter(stlon5s, stlat5s, c='black', marker='o', s=60*station_data_15['AvgLagTime'],edgecolors='k',linewidths=0.1, zorder=4)


        
        legendarray.append(map.scatter([], [], c='cornflowerblue', alpha=0.99, s=60, edgecolors='k'))
        legendarray.append(map.scatter([], [], c='navy', alpha=0.99, s=60, edgecolors='k'))
        legendarray.append(map.scatter([], [], c='black', alpha=0.99, s=60, edgecolors='k'))



        for jj in range(station_data_14.shape[0]):
            plot_point_on_basemap(map, point=(station_data_14['LON'].values[jj],station_data_14['LAT'].values[jj]), angle = station_data_14['AvgFastDir'].values[jj], length = 1.5)

        for jj in range(station_data_4_11.shape[0]):
            plot_point_on_basemap(map, point=(station_data_4_11['LON'].values[jj],station_data_4_11['LAT'].values[jj]), angle = station_data_4_11['AvgFastDir'].values[jj], length = 1.5)

        for jj in range(station_data_15.shape[0]):
            plot_point_on_basemap(map, point=(station_data_15['LON'].values[jj],station_data_15['LAT'].values[jj]), angle = station_data_15['AvgFastDir'].values[jj], length = 1.5)
    else:
        stlon_all,stlat_all = map(sks_meas_all['LON'].values,sks_meas_all['LAT'].values)
        map.scatter(stlon_all,stlat_all, c='black', marker='o', s=60*sks_meas_all['AvgLagTime'],edgecolors='k',linewidths=0.1, zorder=4)

        legendarray.append(map.scatter([], [], c='black', alpha=0.99, s=60, edgecolors='k'))

        for jj in range(sks_meas_all.shape[0]):
            plot_point_on_basemap(map, point=(sks_meas_all['LON'].values[jj],sks_meas_all['LAT'].values[jj]), angle = sks_meas_all['AvgFastDir'].values[jj], length = 1.5)




    #draw mapscale
    msclon,msclat =lonmax - 0.10*np.abs(lonmax-lonmin),latmin+0.10*np.abs(latmax-latmin)
    len_mapscale = 0.15*np.abs(lonmax-lonmin)*111.1
    msclon0,msclat0 = sks_meas_all['LON'].mean(),sks_meas_all['LAT'].mean()
    map.drawmapscale(msclon,msclat,msclon0,msclat0, len_mapscale, barstyle='fancy', zorder=6)
    if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements']):
        if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_no_measurement']) and bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_null_measurements']):
            leg1 = plt.legend([legendarray[3],legendarray[4],legendarray[5],legendarray[6],legendarray[7]],['No measurement','Null measurements',f'{lev1+1}-{lev2-1} measurements',f'{lev2}-{lev3-1} measurements',f'{lev3}+ measurements'],frameon=False, loc='upper left',labelspacing=1,handletextpad=0.1)
        elif bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_no_measurement']) and not bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_null_measurements']):
            leg1 = plt.legend([legendarray[3],legendarray[4],legendarray[5],legendarray[6]],['No measurement',f'{lev1+1}-{lev2-1} measurements',f'{lev2}-{lev3-1} measurements',f'{lev3}+ measurements'],frameon=False, loc='upper left',labelspacing=1,handletextpad=0.1)
        elif not bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_no_measurement']) and bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_null_measurements']):
            leg1 = plt.legend([legendarray[3],legendarray[4],legendarray[5],legendarray[6]],['Null measurement',f'{lev1+1}-{lev2-1} measurements',f'{lev2}-{lev3-1} measurements',f'{lev3}+ measurements'],frameon=False, loc='upper left',labelspacing=1,handletextpad=0.1)
        else:
            leg1 = plt.legend([legendarray[4],legendarray[5],legendarray[6]],[f'{lev1+1}-{lev2-1} measurements',f'{lev2}-{lev3-1} measurements',f'{lev3}+ measurements'],frameon=False, loc='upper left',labelspacing=1,handletextpad=0.1)
    else:
        leg1 = plt.legend([legendarray[3]],['All measurements'],frameon=False, loc='upper left',labelspacing=1,handletextpad=0.1)
    leg2 = plt.legend(frameon=False, loc='upper right',labelspacing=1,handletextpad=0.1)
    ax.add_artist(leg1)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.close('all')
    
def plot_sks_data_nodata_map(sks_meas_all,all_data_df,figname):

    stations = all_data_df['Station'].values
    laat = all_data_df['Latitude'].values
    loon = all_data_df['Longitude'].values
                    
    lonmin, lonmax = loon.min()-1, loon.max()+1
    latmin, latmax = laat.min()-1, laat.max()+1
    if np.abs(lonmax-lonmin)<10 or np.abs(latmax-latmin)<2:
        lonmin, lonmax = lonmin - 2, lonmax + 2
        latmin, latmax = latmin - 2, latmax + 2
    else:
        lonmin, lonmax = lonmin - 0.5, lonmax + 0.5
        latmin, latmax = latmin - 0.5, latmax + 0.5

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    map = Basemap(projection='merc',resolution = 'h', area_thresh = 1000., llcrnrlon=lonmin, llcrnrlat=latmin,urcrnrlon=lonmax, urcrnrlat=latmax)
    
    # plot_topo(map,cmap=plt.cm.rainbow)
    map.etopo(scale=2.5, alpha=0.5, zorder=2)

    map.drawcoastlines(color='k',linewidth=0.5)
    map.drawcountries(color='k',linewidth=0.1)
    map.drawstates(color='gray',linewidth=0.05)
    map.drawrivers(color='blue',linewidth=0.05)
    map.drawmapboundary()
    parallelmin = int(latmin)
    parallelmax = int(latmax)+1
    if np.abs(parallelmax - parallelmin)<5:
        parallelmax += 2
        parallelmin -= 2

    meridianmin = int(lonmin)
    meridianmax = int(lonmax)+1
    if np.abs(meridianmax - meridianmin)<5:
        meridianmax += 2
        meridianmin -= 2

    map.drawparallels(np.arange(parallelmin, parallelmax,dtype='int16').tolist(),labels=[1,0,0,0],linewidth=0)
    map.drawmeridians(np.arange(meridianmin, meridianmax,dtype='int16').tolist(),labels=[0,0,0,1],linewidth=0)
    map.drawmapboundary(color='k', linewidth=2, zorder=1)
    
    
    allstlons,allstlats = map(loon,laat)
    map.scatter(allstlons,allstlats, c='gold', marker='^', s=60, facecolors='none', edgecolors='b',linewidths=0.3, zorder=4)
        
    stlons,stlats = map(sks_meas_all['LON'].values,sks_meas_all['LAT'].values)
    map.scatter(stlons, stlats, c='red', marker='^', s=60, edgecolors='k',linewidths=0.3, zorder=4)

    #draw mapscale
    msclon,msclat =lonmax - 0.10*np.abs(lonmax-lonmin),latmin+0.10*np.abs(latmax-latmin)
    len_mapscale = 0.15*np.abs(lonmax-lonmin)*111.1
    msclon0,msclat0 = np.mean(loon),np.mean(laat)
    map.drawmapscale(msclon,msclat,msclon0,msclat0, len_mapscale, barstyle='fancy', zorder=6)

    legendarray = []
    for col in ['gold','red']:
        legendarray.append(map.scatter([], [], c=col, marker='^', alpha=0.99, s=60, edgecolors='k'))

    leg2 = plt.legend([legendarray[0],legendarray[1]],['No data','With data'],frameon=False, loc='upper right',labelspacing=1,handletextpad=0.1)
    ax.add_artist(leg2)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.close('all')
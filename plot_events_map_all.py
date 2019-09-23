from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)

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
    # cs = map.contourf(lons,lats,etopo,30,cmap=plt.cm.jet,shading='interp')
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


all_stations_file = "/Users/utpalkumar50/Desktop/RF_SKS/results/InfoRF/all_stations_rf_retrieved.txt"
# all_stations_file = "/Users/utpalkumar50/Desktop/RF_SKS/results/InfoSKS/all_stations_sks_retrieved.txt"
method = all_stations_file.split("/")[-1].split("_")[-2].upper()

info_loc=""
for loc in all_stations_file.split("/")[0:-1]:
    info_loc += loc+"/"
all_stations_df = pd.read_csv(all_stations_file,sep="|")
# print(all_stations_df.head())


for net,sta,stalon,stalat  in zip(all_stations_df['#Network'],all_stations_df['Station'],all_stations_df['Longitude'],all_stations_df['Latitude']):
    stalon = -90
    print(f'Plotting events map for {net}-{sta}')
    event_catalog = info_loc + f"{net}-{sta}-events-info-{method}.txt"
    df_all = pd.read_csv(event_catalog,delimiter="\||,", names=['evtime','evlat','evlon','evdp','evmg','client'],header=None,engine="python")
    evmg_all = [float(val.split()[0]) for val in df_all['evmg']]

    df = pd.read_csv(info_loc + f"{net}-{sta}-events-info-available-{method}.txt",delimiter="\||,", names=['evtime','evlat','evlon','evdp','evmg','client'],header=None,engine="python")
    evmg = [float(val.split()[0]) for val in df['evmg']]

  

    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
    plt.figure(figsize=(16,12))
    lon0 = stalon-360 if int(stalon)>0 else stalon

    eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 20000.0,lon_0=lon0)
    eq_map.drawcoastlines(linewidth=0.5)
    eq_map.drawcountries(linewidth=0.1)
    # eq_map.fillcontinents(color = 'gray')
    eq_map.drawmapboundary()
    eq_map.drawparallels(np.arange(-90, 91,30), color = 'k', linewidth=0.1,labels=[1,1,0,0])
    eq_map.drawmeridians(np.arange(-180,180,30), color = 'k', linewidth=0.1,labels=[0,0,1,1])
    plot_topo(eq_map,cmap=plt.cm.rainbow)

    plot_events_loc(eq_map,df_all['evlon'].values,df_all['evlat'].values, evmg_all, df_all['evdp'].values,background=1)
    plot_events_loc(eq_map,df['evlon'].values,df['evlat'].values, evmg, df['evdp'].values,background=0)


    eq_map.scatter(np.NaN, np.NaN, c='g', marker='o', s=100,edgecolors='k',linewidths=0.1, label="Depth < 100 km", zorder=10)
    eq_map.scatter(np.NaN, np.NaN, c='b', marker='o', s=100,edgecolors='k',linewidths=0.1, label=r"100 km $\leq$ Depth < 300 km", zorder=10)
    eq_map.scatter(np.NaN, np.NaN, c='r', marker='o', s=100,edgecolors='k',linewidths=0.1, label=r"Depth $\geq$ 300 km", zorder=10)




    
    # for net in stns_lon:
    x,y = eq_map(stalon,stalat)
    eq_map.plot(x, y,'^', markersize=10,color='k',markeredgecolor='k',linewidth=0.1)
    plt.legend(loc=3)
    for radius in [30,60,90,120]:
        equi(eq_map, np.mean(stalon), np.mean(stalat), radius*DEG2KM,lw=0.1, color='k')

    plt.savefig(f'{net}-{sta}-all_events.png',dpi=200,bbox_inches='tight')
    plt.close('all')


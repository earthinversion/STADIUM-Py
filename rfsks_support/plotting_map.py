from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.plotting_libs import plot_topo


def plot_topo_simple(map,cmap=plt.cm.jet):
    #20 minute bathymetry/topography data
    etopo = np.loadtxt('topo/etopo20data.gz')
    lons  = np.loadtxt('topo/etopo20lons.gz')
    lats  = np.loadtxt('topo/etopo20lats.gz')
    # shift data so lons go from -180 to 180 instead of 20 to 380.
    etopo,lons = shiftgrid(180.,etopo,lons,start=False)
    lons, lats = np.meshgrid(lons, lats)
    # cs = map.contourf(lons,lats,etopo,30,cmap=plt.cm.jet,shading='interp')
    cs = map.pcolormesh(lons,lats,etopo,cmap=cmap,latlon=True,shading='gouraud')

# [Thomas Lecocq](https://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial-09-drawing-circles/)
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

def latlon_grid(bmap, lon_int, lat_int, labels='lb', **kwargs):
    '''Draws a lat-lon grid in an easy way.

    Some default values are taken from rcParams instead of 'black' (color) and
    1.0 (linewidth) which is the default in Basemap.

    In Basemap, the label pad is computed in projection units. Now you can use
    the keyword argument 'labelpad' to control this separation in points. If
    not specified then this value is taken from rcParams.

    Arguments:

    bmap -- Basemap object.
    lon_int, lat_int -- Difference in degrees from one longitude or latitude to
                        the next.
    labels -- String specifying which margins will be used to write the labels.
              If None, no label will be shown.
              It is assummed that left/right margins (i.e. Y axes) correspond
              to latitudes and top/bottom (X axes) to longitudes. It is valid
              every combination of the characters 't' | 'b' | 'l' | 'r'
              (top|bottom|left|right).
              Ex: 'lrb' means that the longitude values will appear in bottom
              margin and latitudes in left and right.
    **kwargs -- Other arguments to drawparallels, drawmeridians and plt.text.
                labelpad has units of points.
    '''
    # Proccesses arguments and rcParams for defult values
    if 'color' not in kwargs:
        kwargs['color'] = plt.rcParams['grid.color']
    if 'linewidth' not in kwargs:
        kwargs['linewidth'] = plt.rcParams['grid.linewidth']
    if 'labelpad' in kwargs:
        padx = pady = kwargs['labelpad']
        del kwargs['labelpad']
    else:
        pady = plt.rcParams['xtick.major.pad']
        padx = plt.rcParams['ytick.major.pad']
    if 'size' in kwargs:
        xfontsize = yfontsize = kwargs['size']
        del kwargs['size']
    elif 'fontsize' in kwargs:
        xfontsize = yfontsize = kwargs['fontsize']
        del kwargs['fontsize']
    else:
        xfontsize = plt.rcParams['xtick.labelsize']
        yfontsize = plt.rcParams['ytick.labelsize']
    # Vectors of coordinates
    lon0 = bmap.lonmin // lon_int * lon_int
    lat0 = bmap.latmin // lat_int * lat_int
    lon1 = bmap.lonmax // lon_int * lon_int
    lat1 = bmap.latmax // lat_int * lat_int
    nlons = (lon1 - lon0) / lon_int + 1
    nlats = (lat1 - lat0) / lat_int + 1
    assert nlons / int(nlons) == 1, nlons
    assert nlats / int(nlats) == 1, nlats
    lons = np.linspace(lon0, lon1, int(nlons))
    lats = np.linspace(lat0, lat1, int(nlats))
    # If not specified then computes de label offset by 'labelpad'
    xos = yos = None
    if 'xoffset' in kwargs:
        xos = kwargs['xoffset']
    if 'yoffset' in kwargs:
        yos = kwargs['yoffset']
    if xos is None and yos is None:
        # Page size in inches and axes limits
        fig_w, fig_h = plt.gcf().get_size_inches()
        points = plt.gca().get_position().get_points()
        x1, y1 = tuple(points[0])
        x2, y2 = tuple(points[1])
        # Width and height of axes in points
        w = (x2 - x1) * fig_w * 72
        h = (y2 - y1) * fig_h * 72
        # If the aspect relation is fixed then compute the real values
        if bmap.fix_aspect:
            aspect = bmap.aspect * w / h
            if aspect > 1:
                w = h / bmap.aspect
            elif aspect < 1:
                h = w * bmap.aspect
        # Offset in projection units (meters or degrees)
        xos = padx * (bmap.urcrnrx - bmap.llcrnrx) / w
        yos = pady * (bmap.urcrnry - bmap.llcrnry) / h
    # Set the labels
    latlabels = [False] * 4
    lonlabels = [False] * 4
    if labels is not None:
        pst = {'l': 0, 'r': 1, 't': 2, 'b': 3}
        lst = {'l': latlabels, 'r': latlabels, 't': lonlabels, 'b': lonlabels}
        for i in labels.lower():
            lst[i][pst[i]] = True
    # Draws the grid
    bmap.drawparallels(lats, labels=latlabels, fontsize=yfontsize,
                       xoffset=xos, yoffset=yos, **kwargs)
    bmap.drawmeridians(lons, labels=lonlabels, fontsize=xfontsize,
                       xoffset=xos, yoffset=yos, **kwargs)


def plot_merc(resolution,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,topo=True):
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
    plt.figure(figsize=(8,8))
    map = Basemap(projection='merc',resolution = resolution, area_thresh = 1000., llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    
    if topo:
        # plot_topo(map,cmap=plt.cm.rainbow)
        plot_topo(map,lonextent=(llcrnrlon,urcrnrlon),latextent=(llcrnrlat,urcrnrlat))

    map.drawcoastlines(color='k',linewidth=0.5)
    # map.fillcontinents()
    map.drawcountries(color='k',linewidth=0.1)
    map.drawmapboundary()
    # ## to fix overlapping of ticks
    # numlatdiv = np.abs(int(urcrnrlat - llcrnrlat))
    # numlondiv = np.abs(int(urcrnrlon - llcrnrlon))
    # if numlatdiv>5:
    #     numlatdiv=5
    # elif numlatdiv<2:
    #     numlatdiv=3

    # if numlondiv>5:
    #     numlondiv=5
    # elif numlondiv<2:
    #     numlondiv=3
    # latlon_grid(map, numlondiv, numlatdiv, labels='lb',linewidth=0, size=6)
    parallelmin = int(llcrnrlat)
    parallelmax = int(urcrnrlat)+1
    if np.abs(parallelmax - parallelmin)<5:
        parallelmax += 2
        parallelmin -= 2

    meridianmin = int(llcrnrlon)
    meridianmax = int(urcrnrlon)+1
    if np.abs(meridianmax - meridianmin)<5:
        meridianmax += 2
        meridianmin -= 2

    map.drawparallels(np.arange(parallelmin, parallelmax,dtype='int16').tolist(),labels=[1,0,0,0],linewidth=0)
    map.drawmeridians(np.arange(meridianmin, meridianmax,dtype='int16').tolist(),labels=[0,0,0,1],linewidth=0)
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

DEG2KM = 111.2
def events_map(evlons, evlats, evmgs, evdps, stns_lon, stns_lat, destination, outname="all-events_map", topo=True,figfrmt='png', clon = 121):
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
    plt.figure(figsize=(16,12))
    lon0 = clon-360 if int(clon)>0 else clon
    eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 20000.0,lon_0=lon0)
    eq_map.drawcoastlines(linewidth=0.5)
    eq_map.drawcountries(linewidth=0.1)
    # eq_map.fillcontinents(color = 'gray')
    eq_map.drawmapboundary()
    eq_map.drawparallels(np.arange(-90, 91,30), color = 'k', linewidth=0.1,labels=[1,1,0,0])
    eq_map.drawmeridians(np.arange(-180,180,30), color = 'k', linewidth=0.1,labels=[0,0,1,1])
    if topo:
        plot_topo_simple(eq_map,cmap=plt.cm.rainbow)
        

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
            

    eq_map.scatter(lo100, lt100, c='g', marker='o', s=mg100,edgecolors='k',linewidths=0.1, zorder=10)
    eq_map.scatter(np.NaN, np.NaN, c='g', marker='o', s=100,edgecolors='k',linewidths=0.1, label="Depth < 100 km", zorder=10)
    eq_map.scatter(lo300, lt300, c='b', marker='o', s=mg300,edgecolors='k',linewidths=0.1, zorder=10)
    eq_map.scatter(np.NaN, np.NaN, c='b', marker='o', s=100,edgecolors='k',linewidths=0.1, label=r"100 km $\leq$ Depth < 300 km", zorder=10)
    eq_map.scatter(loelse, ltelse, c='r', marker='o', s=mgelse,edgecolors='k',linewidths=0.1, zorder=10)
    eq_map.scatter(np.NaN, np.NaN, c='r', marker='o', s=100,edgecolors='k',linewidths=0.1, label=r"Depth $\geq$ 300 km", zorder=10)
    # for net in stns_lon:
    x,y = eq_map(stns_lon, stns_lat)
    eq_map.plot(x, y,'^', markersize=10,color='k',markeredgecolor='k',linewidth=0.1)
    plt.legend(loc=3)
    for radius in [30,60,90,120]:
        equi(eq_map, np.mean(stns_lon), np.mean(stns_lat), radius*DEG2KM,lw=0.1, color='k')

    plt.savefig(destination+outname+"."+figfrmt,dpi=200,bbox_inches='tight')
    plt.close('all')

def plot_bm_azimuth(map,stlon,stlat,distval_lat=180,distval_lon=180,ndivlat=4,ndivlon=5):    
    xvals_az0 = np.linspace(stlon,stlon+distval_lon/111.2,ndivlon)
    yvals_az0 = [stlat for i in range(ndivlon)]
    xmap_az0,ymap_az0 = map(xvals_az0,yvals_az0)

    xvals_az90 = [stlon for i in range(ndivlat)]
    yvals_az90 = np.linspace(stlat,stlat+distval_lat/111.2,ndivlat)
    xmap_az90,ymap_az90 = map(xvals_az90,yvals_az90)

    map.plot(xmap_az0,ymap_az0,'-|', markersize=5,color='dodgerblue',markeredgewidth=2)
    map.plot(xmap_az90,ymap_az90,'-_', markersize=5,color='dodgerblue',markeredgewidth=2)




from mpl_toolkits.basemap import Basemap
from plotting_libs import plot_topo, plot_merc
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt


def plot_point_on_basemap(map, point, angle, length,color='k',alpha=1):
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
    map.plot([startx,endx],[starty,endy],color=color,zorder=3,alpha=alpha)
    

from cmath import rect, phase
from math import radians, degrees
def mean_angle(deg):
    return degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))

sks_db_GM = pd.read_csv('splittingDB_GM.txt',delimiter='\s+', error_bad_lines=False, skiprows=1, names=['id','Station','Latitude',	'Longitude','phi','dt','refID','Phase','pmax','pmin','dtmax','dtmin','remark'])
sks_db_IRIS = pd.read_csv('splittingDB_IRIS.txt',delimiter='\s+', error_bad_lines=False, skiprows=1, names=['id','Station','Latitude',	'Longitude','phi','dt','refID','Phase','pmax','pmin','dtmax','dtmin','remark'])
lblon = np.amin([sks_db_GM['Longitude'].min(),sks_db_IRIS['Longitude'].min()]) - 0.5
lblat = np.amin([sks_db_GM['Latitude'].min(),sks_db_IRIS['Latitude'].min()]) - 0.5
ublon = np.amax([sks_db_GM['Longitude'].max(),sks_db_IRIS['Longitude'].max()]) + 0.5
ublat = np.amax([sks_db_GM['Latitude'].max(),sks_db_IRIS['Latitude'].max()]) + 0.5


print(sks_db_GM.head())

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
map = Basemap(projection='merc',resolution = 'i', area_thresh = 1000., llcrnrlon=lblon, llcrnrlat=lblat,urcrnrlon=ublon, urcrnrlat=ublat, epsg=4839)
# map.drawmapboundary(color='k', linewidth=2, zorder=1)
map.arcgisimage(service='World_Physical_Map', xpixels = 5000)

# map.etopo(scale=2.5, alpha=0.5, zorder=2) # decrease scale (0-1) to downsample the etopo resolution
#The image has a 1" arc resolution
#map.shadedrelief(scale=1, zorder=1)
map.drawcoastlines(color='k',linewidth=0.5,zorder=2)
# map.fillcontinents()
map.drawcountries(color='k',linewidth=0.5,zorder=3)
map.drawstates(color='gray',linewidth=0.05,zorder=3)
map.drawrivers(color='blue',linewidth=0.05,zorder=3)
map.drawparallels(np.linspace(lblat,ublat,6,dtype='int16').tolist(),labels=[1,0,0,0],linewidth=0,zorder=2)
map.drawmeridians(np.linspace(lblon,ublon,6,dtype='int16').tolist(),labels=[0,0,0,1],linewidth=0,zorder=2)

stlonsGM,stlatsGM = map(sks_db_GM['Longitude'].values,sks_db_GM['Latitude'].values)
stlonsIRIS,stlatsIRIS = map(sks_db_IRIS['Longitude'].values,sks_db_IRIS['Latitude'].values)

map.scatter(stlonsGM,stlatsGM, c='blue', marker='o', s=60,edgecolors='k',linewidths=0.1, zorder=4)
for jj in range(sks_db_GM.shape[0]):
    plot_point_on_basemap(map, point=(sks_db_GM['Longitude'].values[jj],sks_db_GM['Latitude'].values[jj]), angle = sks_db_GM['phi'].values[jj], length = 1.5,color='blue')

map.scatter(stlonsIRIS,stlatsIRIS, c='green', marker='o', s=60,edgecolors='k',linewidths=0.1, zorder=4, alpha=0.5)
for jj in range(sks_db_IRIS.shape[0]):
    plot_point_on_basemap(map, point=(sks_db_IRIS['Longitude'].values[jj],sks_db_IRIS['Latitude'].values[jj]), angle = sks_db_IRIS['phi'].values[jj], length = 1.5,color='green',alpha=0.5)

##legend
map.scatter([], [], c='blue', alpha=0.6, s=60,label="GM",edgecolors='k')
map.scatter([], [], c='green', alpha=0.6, s=60,label="IRIS",edgecolors='k')
plt.legend(loc='upper left')
plt.savefig('SKS_previous_studies.png',bbox_inches='tight',dpi=300)



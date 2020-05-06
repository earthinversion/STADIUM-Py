from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import glob, math
import pandas as pd
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.plotting_libs import shoot, equi, plot_topo, plot_events_loc, plot_merc
import os

# plot_h = 1
# plot_kappa = 1

def plot_h_kappa(h_k_file = "rfsks_support/h-kappa-values.txt",all_stationsfile = "results/InfoRF/all_stations_rf_retrieved.txt",plot_h = 1,plot_kappa = 1):
    fig_loc = all_stationsfile.split("/")[:-1]
    fig_loc_str = "/".join(fig_loc)


    ## read the stations inventory
    h_kappa_df = pd.read_csv(h_k_file,sep=",",header=None,names=['network','station', 'latitude', 'longitude','thickness','kappa'])
    if h_kappa_df.shape[0]:
        all_stations_df = pd.read_csv(all_stationsfile,sep="|")
        outfig_h = fig_loc_str+"/"+'all_stations_thickness_map.png'
        outfig_k = fig_loc_str+"/"+'all_stations_kappa_map.png'
        if plot_h and not os.path.exists(outfig_h):
            fig = plt.figure(figsize=(20,12))
            map = plot_merc(resolution='h',llcrnrlon=int(all_stations_df['Longitude'].min())-1, llcrnrlat=int(all_stations_df['Latitude'].min())-1,urcrnrlon=math.ceil(all_stations_df['Longitude'].max())+1, urcrnrlat=math.ceil(all_stations_df['Latitude'].max())+1,topo=False)


            stalon = h_kappa_df['longitude'].values
            stalat = h_kappa_df['latitude'].values
            x,y = map(stalon, stalat)
            sc = plt.scatter(x, y, c=h_kappa_df['thickness'].values, s=30,edgecolors='k',linewidth=0.1)
            cb = plt.colorbar(sc,fraction=0.025, pad=0.04)
            cb.ax.set_ylabel('Thickness in km',labelpad=15,fontsize=18)

            x1,y1 = map(all_stations_df['Longitude'].values, all_stations_df['Latitude'].values)
            plt.scatter(x1, y1, s=30, facecolors='none', edgecolors='k',linewidth=0.1)
            plt.tight_layout()

                
            # plt.legend(loc=4,fontsize=8)
            plt.savefig(outfig_h,dpi=300,bbox_inches='tight')
            plt.close('all')

        if plot_kappa and not os.path.exists(outfig_k):
            fig = plt.figure(figsize=(20,16))
            map = plot_merc(resolution='h',llcrnrlon=int(all_stations_df['Longitude'].min())-1, llcrnrlat=int(all_stations_df['Latitude'].min())-1,urcrnrlon=math.ceil(all_stations_df['Longitude'].max())+1, urcrnrlat=math.ceil(all_stations_df['Latitude'].max())+1,topo=False)


            stalon = h_kappa_df['longitude'].values
            stalat = h_kappa_df['latitude'].values
            x,y = map(stalon, stalat)
            sc = plt.scatter(x, y, c=h_kappa_df['kappa'].values, s=30,edgecolors='k',linewidth=0.1)
            cb = plt.colorbar(sc,fraction=0.025, pad=0.04)
            cb.set_label(r'$\frac{V_p}{V_s}$',labelpad=15,rotation=0,fontsize=18)

            x1,y1 = map(all_stations_df['Longitude'].values, all_stations_df['Latitude'].values)
            plt.scatter(x1, y1, s=30, facecolors='none', edgecolors='k',linewidth=0.1)
            plt.tight_layout()
            

                
            # plt.legend(loc=4,fontsize=8)
            plt.savefig(outfig_k,dpi=300,bbox_inches='tight')
            plt.close('all')
if __name__=="__main__":
    plot_h_kappa()
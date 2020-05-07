from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import glob,os
import pandas as pd
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.plotting_libs import equi, plot_events_loc
from rfsks_support.plotting_map import plot_topo_simple
DEG2KM = 111.2

def plot_events_map_all(all_stations_file = "results/InfoRF/all_stations_rf_retrieved.txt"):
    # all_stations_file = "/Users/utpalkumar50/Desktop/RF_SKS/results/InfoSKS/all_stations_sks_retrieved.txt"
    method = all_stations_file.split("/")[-1].split("_")[-2].upper()
    fig_loc = all_stations_file.split("/")[:-1]
    fig_loc_str = "/".join(fig_loc)
    

    info_loc=""
    for loc in all_stations_file.split("/")[0:-1]:
        info_loc += loc+"/"
    all_stations_df = pd.read_csv(all_stations_file,sep="|")

    net_sta_stalon_stalat_list = []
    for net,sta,stalon,stalat  in zip(all_stations_df['#Network'],all_stations_df['Station'],all_stations_df['Longitude'],all_stations_df['Latitude']):
        net_sta_stalon_stalat = f"{net}_{sta}_{stalon}_{stalat}"
        if net_sta_stalon_stalat not in net_sta_stalon_stalat_list:
            net_sta_stalon_stalat_list.append(net_sta_stalon_stalat)

    for net_sta_stalon_stalat  in net_sta_stalon_stalat_list:
        net = net_sta_stalon_stalat.split("_")[0]
        sta = net_sta_stalon_stalat.split("_")[1]
        figname = f'{fig_loc_str}/{net}-{sta}-all_events-{method}.png'
        if not os.path.exists(figname):
            stalon = float(net_sta_stalon_stalat.split("_")[2])
            stalat = float(net_sta_stalon_stalat.split("_")[3])
            # stalon = -90
            print(f'Plotting events map for {net}-{sta}')
            event_catalog = info_loc + f"{net}-{sta}-events-info-{method}.txt"
            df_all = pd.read_csv(event_catalog)
            evmg_all = df_all['evmg'].values
            # evmg_all = [float(val.split()[0]) for val in df_all['evmg']]

            df = pd.read_csv(info_loc + f"{net}-{sta}-events-info-available-{method}.txt",delimiter="\||,", names=['evtime','evlat','evlon','evdp','evmg','client'],header=None,engine="python")
            # evmg = df['evmg'].values
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
            plot_topo_simple(eq_map,cmap=plt.cm.rainbow)

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

            plt.savefig(figname,dpi=200,bbox_inches='tight')
            plt.close('all')
if __name__=="__main__":
    plot_events_map_all()

from mpl_toolkits.basemap import Basemap,shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import glob, os
import pandas as pd
import warnings, matplotlib.cbook
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.plotting_libs import shoot, equi, plot_topo, plot_events_loc, plot_merc


def rem_duplicate_stations(retr_stations_df):
    net_sta_stalon_stalat_list = []
    for net,sta,stalon,stalat  in zip(retr_stations_df['#Network'],retr_stations_df['Station'],retr_stations_df['Longitude'],retr_stations_df['Latitude']):
        net_sta_stalon_stalat = f"{net}_{sta}_{stalon}_{stalat}"
        if net_sta_stalon_stalat not in net_sta_stalon_stalat_list:
            net_sta_stalon_stalat_list.append(net_sta_stalon_stalat)
    return net_sta_stalon_stalat_list

def plot_station_map_all(retr_stationsfile = "results/InfoRF/all_stations_rf_retrieved.txt",all_stationsfile = "results/InfoRF/all_stations_RF.txt"):
    # ## all retrieved stations
    # retr_stationsfile = "results/InfoRF/all_stations_rf_retrieved.txt"

    # ## all stations retrieved or not
    # all_stationsfile = "results/InfoRF/all_stations_RF.txt"

    fig_loc = all_stationsfile.split("/")[:-1]
    fig_loc_str = "/".join(fig_loc)

    figname = fig_loc_str+"/"+'all_stations_map.png'
    if not os.path.exists(figname):
        ## selected method automatically from the filename
        method = retr_stationsfile.split("/")[-1].split("_")[-2].upper()

        ## extract the file path
        info_loc=""
        for loc in retr_stationsfile.split("/")[0:-1]:
            info_loc += loc+"/"

        ## read the stations inventory
        retr_stations_df = pd.read_csv(retr_stationsfile,sep="|")
        all_stations_df = pd.read_csv(all_stationsfile,sep="|")

        net_sta_stalon_stalat_list = rem_duplicate_stations(retr_stations_df)
        net_sta_stalon_stalat_list_all = rem_duplicate_stations(all_stations_df)

        new_net_sta_stalon_stalat_list = []
        for net_sta_stalon_stalat_val in net_sta_stalon_stalat_list_all:
            if net_sta_stalon_stalat_val in net_sta_stalon_stalat_list:
                new_net_sta_stalon_stalat = f"{net_sta_stalon_stalat_val}_1"
                new_net_sta_stalon_stalat_list.append(new_net_sta_stalon_stalat)
            else:
                new_net_sta_stalon_stalat = f"{net_sta_stalon_stalat_val}_0"
                new_net_sta_stalon_stalat_list.append(new_net_sta_stalon_stalat)


        all_nets_set =  set(net_sta_stalon_stalat_val.split("_")[0] for net_sta_stalon_stalat_val in new_net_sta_stalon_stalat_list)
        color=iter(plt.cm.viridis(np.linspace(0,1,len(all_nets_set))))
        netcolor = {net:c for net, c in zip(all_nets_set, color)}




        map = plot_merc(resolution='h',llcrnrlon=all_stations_df['Longitude'].min()-1, llcrnrlat=all_stations_df['Latitude'].min()-1,urcrnrlon=all_stations_df['Longitude'].max()+1, urcrnrlat=all_stations_df['Latitude'].max()+1,topo=True)

        for net_sta_stalon_stalat_val in new_net_sta_stalon_stalat_list:
            net = net_sta_stalon_stalat_val.split("_")[0]
            sta = net_sta_stalon_stalat_val.split("_")[1]
            stalon = float(net_sta_stalon_stalat_val.split("_")[2])
            stalat = float(net_sta_stalon_stalat_val.split("_")[3])
            avail = int(net_sta_stalon_stalat_val.split("_")[4])
            # stalon = -90
            print(f'for {net}-{sta}')

            x,y = map(stalon, stalat)
            if avail:
                map.plot(x, y,'^', color=netcolor[net], markersize=7,markeredgecolor=netcolor[net],linewidth=0.1)
            else:
                map.plot(x, y,'^', markerfacecolor="None", markersize=7,markeredgecolor=netcolor[net], markeredgewidth=5)
            
        for net in all_nets_set:
            map.plot(np.NaN,np.NaN,'^',color=netcolor[net],label=net, markersize=7,markeredgecolor='k',linewidth=0.1)
        plt.legend(loc=4,fontsize=8)
        plt.savefig(figname,dpi=300,bbox_inches='tight')
        plt.close('all')

if __name__=="__main__":
    plot_station_map_all()
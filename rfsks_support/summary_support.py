import os, glob
import pandas as pd
import numpy as np
from datetime import datetime

class sum_support:
    def __init__(self,sum_file,res_dir):
        if not os.path.exists(sum_file):
                sum_file_id = open(sum_file,'w')
                sum_file_id.write(f"This is auto-generated summary file for the project: {res_dir}\n")
                sum_file_id.write(f"Created at {datetime.now()}")
                sum_file_id.write("\n")
                mode = 'write'
        else:
            sum_file_id = open(sum_file,'a')
            sum_file_id.write(f"\n########## APPENDING ({datetime.now()}) ##########\n")
            mode = 'append'

        self.sum_file_id = sum_file_id
        self.mode = mode

    def write_initial_summary(self,mnlong,mxlong,mnlat,mxlat,client,network,station,channel):
        if self.mode == 'write':
            self.sum_file_id.write("Requested geograhical coordinates:\nminmax longitude: {:.4f} {:.4f}\nminmax latitude: {:.4f} {:.4f}\n".format(mnlong,mxlong,mnlat,mxlat))
            self.sum_file_id.write("\n")
            self.sum_file_id.write("Requested client: {}\n".format(client))
            if network == "*":
                self.sum_file_id.write("Requested all available networks\n")
            else:
                self.sum_file_id.write("Requested network: {}\n".format(network))

            if station == "*":
                self.sum_file_id.write("Requested all available stations\n")
            else:
                self.sum_file_id.write("Requested station: {}\n".format(station))

            self.sum_file_id.write("Requested channel: {}\n".format(channel))
            self.sum_file_id.write("\n")

    def makeSKSRF(self,makeRF,makeSKS):
        if self.mode == 'write':
            if makeRF and not makeSKS:
                self.sum_file_id.write("Requested Receiver Functions Computations only\n")
            elif not makeRF and makeSKS:
                self.sum_file_id.write("Requested Shear-wave Splitting Computations only\n")
            elif makeRF and makeSKS:
                self.sum_file_id.write("Requested Receiver Functions & Shear-wave Splitting Computations\n")
            self.sum_file_id.write("\n")

    def write_strings(self,string):
        self.sum_file_id.write(string+"\n")

    def write_data_summary(self,SKSsta):
        self.write_strings("----> Data summary")
        self.write_strings("------> Inventory Information:")
        self.SKSsta_path = "/".join(SKSsta.split("/")[0:-1])
        SKSsta = SKSsta.split("/")[-1].split(".")[0]
        full_SKSsta_path = f"{self.SKSsta_path}/{SKSsta}_combined.txt"
        full_SKSsta_path_unfor = f"{self.SKSsta_path}/{SKSsta}.txt"

        ## Station info
        df_avail_station = pd.read_csv(full_SKSsta_path,sep="|")
        self.write_strings("Total available stations in the region: {}".format(df_avail_station.shape[0]))
        self.write_strings("See files: {} [uncombined] and {} [combined events for same station] for details".format(full_SKSsta_path_unfor, full_SKSsta_path))

        ## Corresponding events info
        self.write_strings("")
        all_events_file_list = glob.glob(f"{self.SKSsta_path}/*-*-*-*-events-info-*.txt")
        
        total_events = 0
        minmg = 100
        maxmg = 0
        mindp, maxdp = 99999,-10
        min_evs,max_evs = 999999,-10
        for ff in all_events_file_list:
            dff = pd.read_csv(ff,sep=",")
            total_events+=dff.shape[0]
            if dff.shape[0] < min_evs:
                min_evs = dff.shape[0]
                min_ev_catlog = ff
            if dff.shape[0] > max_evs:
                max_evs = dff.shape[0]
                max_ev_catlog = ff
            ## magnitude
            mag_array = dff['evmg'].values
            minmag = np.min(mag_array)
            maxmag = np.max(mag_array)
            if minmag < minmg:
                minmg = minmag
            if maxmag > maxmg:
                maxmg = maxmag
            ## depth
            mindep = dff['evdp'].min()
            maxdep = dff['evdp'].max()
            if mindep < mindp:
                mindp = mindep
            if maxdep > maxdp:
                maxdp = maxdep


        self.write_strings("Total event catalogs generated: {}".format(len(all_events_file_list)))
        self.write_strings("Total events in total: {}".format(total_events))
        self.write_strings("Catalog with least events: {} ({} events)".format(min_ev_catlog,min_evs))
        self.write_strings("Catalog with most events: {} ({} events)".format(max_ev_catlog,max_evs))
        self.write_strings("Min mag: {}, Max mag: {}".format(minmg,maxmg))
        self.write_strings("Min depth: {} km, Max depth: {} km".format(mindp, maxdp))
    
    def write_data_download_summary(self,datafileloc,retrived_stn_file):
        ## Download info
        self.write_strings("")
        self.write_strings("------> Download Information:")
        total_data_files = len(glob.glob(datafileloc+"*.h5"))

        self.write_strings("Total downloaded data files: {}".format(total_data_files))
        self.write_strings("Stations for which the data has been retrieved successfully: {}".format(retrived_stn_file))
        self.write_strings("Downloaded station map: {}".format(self.SKSsta_path+"/all_stations_map.png"))
        self.write_strings("Downloaded events map (one for each retrieved station): {}".format(self.SKSsta_path+'/net-sta-all_events-method.png'))



    def close_sumfile(self):
        self.sum_file_id.close()


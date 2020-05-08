import os, glob
import pandas as pd
import numpy as np
from datetime import datetime
import yaml

with open('Settings/advSKSparam.yaml') as f:
    inpSKSdict = yaml.load(f, Loader=yaml.FullLoader)

with open('Settings/advRFparam.yaml') as f:
    inpRFdict = yaml.load(f, Loader=yaml.FullLoader)

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
            for i in range(10):
                sum_file_id.write("\n")
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
        if makeRF and not makeSKS:
            self.sum_file_id.write("Requested Receiver Functions Computations only\n")
        elif not makeRF and makeSKS:
            self.sum_file_id.write("Requested Shear-wave Splitting Computations only\n")
        elif makeRF and makeSKS:
            self.sum_file_id.write("Requested Receiver Functions & Shear-wave Splitting Computations\n")
        self.sum_file_id.write("\n")

    def write_strings(self,string):
        self.sum_file_id.write(string+"\n")

    def newline(self):
        self.write_strings("")

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
        self.newline()
        all_events_file_list = glob.glob(f"{self.SKSsta_path}/*-*-*-*-events-info-*.txt")
        
        total_events = 0
        minmg = 100
        maxmg = 0
        mindp, maxdp = 99999,-10
        min_evs,max_evs = 999999,-10
        for ff in all_events_file_list:
            dff = pd.read_csv(ff,sep=",")
            if dff.shape[0]>0:
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
    
    def write_data_download_summary(self,datafileloc,retrived_stn_file,method='SKS'):
        ## Download info
        self.newline()
        self.write_strings("------> Download Information:")
        if method == 'SKS':
            total_data_files = len(glob.glob(datafileloc+"*.h5"))
        elif method == 'RF':
            total_data_files = len(glob.glob(datafileloc+f"*-{str(inpRFdict['filenames']['data_rf_suffix'])}.h5"))


        self.write_strings("Total downloaded data files: {}".format(total_data_files))
        self.write_strings("Stations for which the data has been retrieved successfully: {}".format(retrived_stn_file))
        self.write_strings("Downloaded station map: {}".format(self.SKSsta_path+"/all_stations_map.png"))
        self.write_strings("Downloaded events map (one for each retrieved station): {}".format(self.SKSsta_path+'/net-sta-all_events-method.png'))

        ## finding min and max time of operation of stations
        df_retr = pd.read_csv(retrived_stn_file,sep="|",keep_default_na=False, na_values=[""])
        df_retr['EndTime'].fillna('2599-12-31T23:59:59',inplace=True)
        df_retr['StartTimeNum'] = df_retr['StartTime'].apply(lambda x: int(x.split("-")[0]+x.split("-")[1]+x.split("-")[2][0:2]))
        df_retr['EndTimeNum'] = df_retr['EndTime'].apply(lambda x: int(x.split("-")[0]+x.split("-")[1]+x.split("-")[2][0:2]))
        df_retr['startend_dur'] = df_retr['EndTimeNum'].values-df_retr['StartTimeNum'].values

        startmin_row = df_retr.loc[df_retr['StartTimeNum'].idxmin()]
        endmax_row = df_retr.loc[df_retr['EndTimeNum'].idxmax()]
        maxdur_row = df_retr.loc[df_retr['startend_dur'].idxmax()]

        self.newline()
        self.write_strings("Minimum starttime for the retrieved stations {} ({}-{})".format(startmin_row['StartTime'],startmin_row['#Network'],startmin_row['Station']))
        self.write_strings("Maximum endtime for the retrieved stations {} ({}-{})".format(endmax_row['EndTime'],endmax_row['#Network'],endmax_row['Station']))
        self.write_strings("Longest operating station: {}-{}".format(maxdur_row['#Network'],maxdur_row['Station']))

    def write_sks_meas_sum(self,measure_loc,trace_loc_ENZ,trace_loc_RTZ,trigger_loc):
        self.newline()
        self.write_strings("----> SKS measurement summary:")
        measure_loc_all = "/".join(measure_loc.split("/")[:-2])
        measure_loc_null = "/".join(measure_loc.split("/")[:-1])

        all_measurements_file = measure_loc_all+"/sks_measurements_all.txt"
        self.write_strings("Summary of SKS measurements for all stations stored in: {}".format(all_measurements_file))

        self.write_strings("Summary of SKS measurements for individual stations stored in: {}".format(measure_loc_null+"/net_sta_sks_measurements.txt"))
        self.write_strings("Summary of null measurements for individual stations stored in: {}".format(measure_loc_null+"/net_sta_null_measurements.txt"))

        #read measurement summary file
        self.newline()
        df_meas_sum = pd.read_csv(all_measurements_file,sep='\s+')

        max_num_meas = df_meas_sum.loc[df_meas_sum['NumMeasurements'].idxmax()]
        max_null_meas = df_meas_sum.loc[df_meas_sum['NumNull'].idxmax()]

        self.write_strings("Successful measurements for {} stations".format(df_meas_sum.shape[0]))
        self.write_strings("Max number of measurements for station: {}-{} ({:.4f},{:.4f})".format(max_num_meas['NET'],max_num_meas['STA'],max_num_meas['LAT'],max_num_meas['LON']))
        self.write_strings("Max number of null measurements for station: {}-{} ({:.4f},{:.4f})".format(max_null_meas['NET'],max_null_meas['STA'],max_null_meas['LAT'],max_null_meas['LON']))

        ## plot traces summary
        if trace_loc_ENZ:
            self.newline()
            self.write_strings("ENZ plot for traces stored at {}".format(trace_loc_ENZ))

        if trace_loc_RTZ:
            self.newline()
            self.write_strings("RTZ plot for traces stored at {}".format(trace_loc_RTZ))

        ## Picking
        self.newline()
        self.write_strings("SKS picking algorithm: {} (thresholds: {:.2f},{:.2f})".format(str(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo']),float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1'])))
        if trigger_loc:
            self.write_strings("SKS phase picking plot stored as: {}".format(trigger_loc+'*-eventtime-trigger.png'))

        ## Measurement constrains used
        self.newline()
        self.write_strings("------> Measurement contrains:")
        if str(inpSKSdict['sks_measurement_contrains']['sel_param']) == "lam12":
            self.write_strings("Eigenvalue ratio (lambda1/lambda2), threholds:: fast direction: {:.1f} lag time: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['lam12fast_threh']),float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['lam12lag_threh'])))
        elif str(inpSKSdict['sks_measurement_contrains']['sel_param']) == "snr":
            self.write_strings("signal to noise ratio, threholds:: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['snr_ratio'])))

        self.write_strings("Min allowed lag: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['lag_settings']['minlag'])))
        self.write_strings("Max allowed lag: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['lag_settings']['maxlag'])))
        self.write_strings("Max allowed lag error: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['lag_settings']['maxdlag'])))
        self.write_strings("Max allowed fast direction error: {:.1f}".format(float(inpSKSdict['sks_measurement_contrains']['fast_dir_settings']['maxdfast'])))

        ## Measurement snapshot
        if bool(inpSKSdict['sks_measurement_plot']['measurement_snapshot']):
            self.newline()
            self.write_strings("Measurement snapshot saved as: {}".format(measure_loc_null+'*-evyear_evmonth_evday_evhour_evminute.png'))

        if int(inpSKSdict['error_plot_toggles']['error_plot_indiv']):
            self.newline()
            self.write_strings("Individual error plots saved as: {}".format(measure_loc_null+'errorplot_*-evyear_evmonth_evday_evhour_evminute.png'))

        if bool(inpSKSdict['error_plot_toggles']['error_plot_all']):
            self.newline()
            self.write_strings("Summary error plot for all measurements saved as: {}".format(measure_loc_null+'errorplot_*.png'))

        if bool(inpSKSdict['sks_measurement_plot']['plot_SI']):
            self.newline()
            self.write_strings("Splitting intensity plot w.r.t backazimuth stored at: {}".format(measure_loc_null+'net_sta_BAZ_SI.png'))

        self.newline()
        self.write_strings("SKS measurements map: {}".format(measure_loc_all+'/SKS_station_Map.png'))
        if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements']):
            plot_params_lev = inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['meas_seg_points']
            lev1, lev2, lev3 = int(plot_params_lev['lev1']), int(plot_params_lev['lev2']), int(plot_params_lev['lev3'])
            self.write_strings("Measurement segregated: {},{},{}".format(f'{lev1+1}-{lev2-1} measurements',f'{lev2}-{lev3-1} measurements',f'{lev3}+ measurements'))
            if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_no_measurement']):
                self.write_strings("No measurement shown")
            if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['show_null_measurements']):
                self.write_strings("Null measurement shown")

        self.newline()
        self.write_strings("Map for stations with and without data stored at: {}".format(measure_loc_all+'/data_nodata_map.png'))

    ######## Receiver functions summary ###########
    def write_rf_comp_summary(self,datafileloc,destImg):
        self.newline()
        self.write_strings("----> RF computation summary:")
        total_computed_files = len(glob.glob(datafileloc+f"*-{str(inpRFdict['filenames']['rf_compute_data_suffix'])}.h5"))
        
        self.write_strings("Total number of RF computations: {}".format(total_computed_files))
        self.write_strings("RF computations stored as: {}".format(datafileloc+f"*-{str(inpRFdict['filenames']['rf_compute_data_suffix'])}.h5"))

        self.newline()
        self.write_strings("Bandpass filtered between {:.1f} and {:.1f} for RF computation".format(float(inpRFdict['rf_filter_settings']['minfreq']),float(inpRFdict['rf_filter_settings']['maxfreq'])))

        self.newline()
        self.write_strings("RF plots for L and Q components stored at: {}".format(destImg))

    def write_rf_pp_summary(self,datafileloc,destImg):
        self.newline()
        outimagename = destImg+'piercing_points_map.png'
        outimagename_new = destImg+'piercing_points_map_new.png'
        outputfile = datafileloc+ f"{str(inpRFdict['filenames']['rfprofile_compute_result_prefix'])}AZM_INITDIV_ENDDIV_WIDTH-PROFILE_NUM.h5"
        self.write_strings("RF piercing points calculation for computing RF profile stored as: {}".format(outputfile))
        self.write_strings("RF piercing points map stored at: {}".format(outimagename))
        self.write_strings("RF piercing points map with guides stored at: {}".format(outimagename_new))

        self.newline()
        rf_profile = destImg+"COMP_AZM_NUM_profile.png"
        self.write_strings("RF profile stored at: {}".format(rf_profile))
        self.write_strings("RF profile calculated for the trim range: {:.1f}, {:.1f}".format(int(inpRFdict['rf_display_settings']['trim_min']), int(inpRFdict['rf_display_settings']['trim_max'])))

    def h_kappa_summary(self,figloc,h_k_calc_file):
        self.newline()
        self.write_strings("H-Kappa calculation results: {}".format(figloc+h_k_calc_file))

        outfig_h = figloc+'all_stations_thickness_map.png'
        outfig_k = figloc+'all_stations_kappa_map.png'

        if inpRFdict['filenames']['h_kappa_settings']['plot_h']:
            self.write_strings("Crustal plot: {}".format(outfig_h))
        if inpRFdict['filenames']['h_kappa_settings']['plot_kappa']:
            self.write_strings("Crustal plot: {}".format(outfig_k))



    def close_sumfile(self):
        self.sum_file_id.close()


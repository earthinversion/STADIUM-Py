import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# plt.style.use('ggplot')
plt.style.use('seaborn')
import matplotlib.gridspec as gridspec
import tqdm
import os, glob, yaml
from rf import read_rf, IterMultipleComponents
from obspy import UTCDateTime as UTC
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from obspy.core import read
from obspy.taup import TauPyModel
from rfsks_support.other_support import measure_status, sks_measure_file_start
from rfsks_support.rfsks_extras import plot_trigger, plot_trace, plot_SKS_measure, filter_pick_snr, filter_pick_lam12, errorplot, errorplot_all, auto_null_measure, polar_error_surface, splitting_intensity, segregate_measurements, plot_baz_si_map
from obspy.signal.trigger import recursive_sta_lta,classic_sta_lta,z_detect,carl_sta_trig,delayed_sta_lta, trigger_onset
import splitwavepy as sw
import logging
import math
from mpl_toolkits.basemap import Basemap
from rfsks_support.plotting_libs import plot_topo, plot_merc, plot_point_on_basemap, mean_angle, plot_sks_station_map, plot_sks_data_nodata_map


## Fine tuning of SKS
with open('Settings/advSKSparam.yaml') as f:
    inpSKSdict = yaml.load(f, Loader=yaml.FullLoader)

class sks_measurements:

    def __init__(self,plot_measure_loc=None):
        self.logger = logging.getLogger(__name__)
        self.plot_measure_loc= plot_measure_loc
        # print(self.plot_measure_loc)
        # pass

    ## Pre-processing
    def SKScalc(self, dataSKSfileloc,trace_loc_ENZ=None,trace_loc_RTZ=None,trigger_loc=None,method = 'None'):
        
        self.logger.info("Cut the traces around the SKS arrival")
        sksfiles = glob.glob(dataSKSfileloc+f"*-{str(inpSKSdict['filenames']['data_sks_suffix'])}.h5")
        # self.logger.info(sksfiles)
        
        # all_measurements = open(self.plot_measure_loc+"../"+"sks_measurements_all.txt",'w')
        # all_measurements.write("NET STA LON LAT AvgFastDir AvgLagTime NumMeasurements NumNull\n")
        all_meas_start,all_meas_close = True, False

        meas_file = self.plot_measure_loc+'done_measurements.txt'
        f, finished_file, finished_events = measure_status(meas_file) #track the measurements
        
        for i,sksfile in enumerate(sksfiles):
            count=0
            data = read_rf(sksfile, 'H5')
            self.logger.info(f"Calculating SKS arrival times for {sksfile}\n")
            net_name = os.path.basename(sksfile).split("-")[0]
            stn_name = os.path.basename(sksfile).split("-")[1]

            stn_meas_close = False
            # if stn_meas_start:
            sks_measurements_stn = self.plot_measure_loc+f"{net_name}_{stn_name}_{str(inpSKSdict['filenames']['sks_meas_indiv'])}"
            null_measurements_stn = self.plot_measure_loc+f"{net_name}_{stn_name}_null_measurements.txt"
            if not os.path.exists(sks_measurements_stn):
                sks_meas_file = sks_measure_file_start(sks_measurements_stn,data[0].stats.station_longitude,data[0].stats.station_latitude,"EventTime EvLong EvLat Evdp Baz FastDirection(degs) deltaFastDir(degs) LagTime(s) deltaLagTime(s) SI\n")
                
                sks_meas_file_null = sks_measure_file_start(null_measurements_stn,data[0].stats.station_longitude,data[0].stats.station_latitude,"EventTime EvLong EvLat Evdp Baz\n")
                stn_meas_close = True

            plt_id=f"{net_name}-{stn_name}"
            measure_list,squashfast_list,squashlag_list=[],[],[]
            fast_dir_all, lag_time_all = [], []
            num_measurements, num_null = 0, 0
            for stream3c in IterMultipleComponents(data, 'onset', 3):
                count+=1
                ## check if the length of all three traces are equal
                for tr in stream3c:
                    lentr=tr.stats.npts
                    lengt= tr.stats.sampling_rate * 100
                    if lentr != lengt:
                        # print('Wrong trace length ', lentr,lengt)
                        continue
                    
                if sksfile in finished_file and str(stream3c[0].stats.event_time) in finished_events:
                    continue
                else:
                    if all_meas_start:
                        all_measurements = open(self.plot_measure_loc+"../"+"sks_measurements_all.txt",'w')
                        all_measurements.write("NET STA LON LAT AvgFastDir AvgLagTime NumMeasurements NumNull\n")
                        all_meas_start = False
                        all_meas_close = True
                    f.write("{},{}\n".format(sksfile,stream3c[0].stats.event_time))

                ## check if the length of all three traces are equal
                len_tr_list=list()
                for tr in stream3c:
                    len_tr_list.append(len(tr))
                if len(set(len_tr_list))!=1:
                    self.logger.warning(f"{count}/{int(len(data)/3)} Bad trace: {stream3c[0].stats.event_time}")
                    continue

                ## filter the trace
                st = stream3c.filter('bandpass', freqmin=float(inpSKSdict['sks_filter_settings']['minfreq']), freqmax=float(inpSKSdict['sks_filter_settings']['maxfreq']))
                st.detrend('linear')
                # st.taper(max_percentage=0.05, type="hann")
                sps = st[0].stats.sampling_rate
                t = st[0].stats.starttime
                ## trim the trace
        
                trace1 = st.trim(t+int(inpSKSdict['sks_picking']['trimstart']), t+int(inpSKSdict['sks_picking']['trimend']))



                ## plot the ENZ
                if trace_loc_ENZ:
                    plot_trace(trace1,trace_loc_ENZ)
                
                ## Rotate to RTZ
                ## trace2[0]->BHT; trace2[1]->BHR; trace2[2]->BHZ;
                trace1.rotate('NE->RT')


                evyear = trace1[0].stats.event_time.year
                evmonth = trace1[0].stats.event_time.month
                evday = trace1[0].stats.event_time.day
                evhour = trace1[0].stats.event_time.hour
                evminute = trace1[0].stats.event_time.minute
                
                # ## plot all three traces RTZ
                if trace_loc_RTZ:
                    plot_trace(trace1,trace_loc_RTZ)

                ######################
                #  Different picker methods
                ######################
                ### operating on transverse component
                if method=="recursive_sta_lta":
                    # self.logger.info(f"Method is {method}")
                    cft = recursive_sta_lta(trace1[1].data, int(1 * sps), int(5 * sps))
                    threshold = (float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1']))#(2.5,0.65)
                    on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                    
                    if trigger_loc and on_off.shape[0]==1:
                        outfile = trigger_loc+f'{plt_id}-{trace1[0].stats.event_time}-trigger.png'
                        plot_trigger(trace1[1], cft, on_off, threshold[0], threshold[1], outfile=outfile)
                    

                elif method=="classic_sta_lta":
                    cft = classic_sta_lta(trace1[1].data, int(5 * sps), int(10 * sps))
                    threshold = (float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1']))#(1.5, 0.5)
                    on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                elif method=="z_detect":
                    cft = z_detect(trace1[1].data, int(10 * sps))
                    threshold = (float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1']))#(-0.4, -0.3)
                    on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                elif method=="carl_sta_trig":
                    cft = carl_sta_trig(trace1[1].data, int(5 * sps), int(10 * sps), 0.8, 0.8)
                    threshold = (float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1']))#(20.0, -20.0)
                    on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                elif method=="delayed_sta_lta":
                    cft = delayed_sta_lta(trace1[1].data, int(5 * sps), int(10 * sps))
                    threshold = (float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr0']), float(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo_thr1']))#(5, 10)
                    on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                else:
                    self.logger.info("No valid method specified")
                    pass

                if on_off.shape[0]==1:
                    trace1.rotate('RT->NE')
                    trace2 = trace1
                    realdata = sw.Pair(trace2[1].data,trace2[0].data, delta=1/sps) #creates Pair from two traces, delta: sample interval
                    try:
                        measure = sw.EigenM(realdata, lags=(float(inpSKSdict['sks_measurement_contrains']['lag_settings']['minlag']), float(inpSKSdict['sks_measurement_contrains']['lag_settings']['maxlag']), 40))
                        
                    except Exception as e:
                        self.logger.error(e)
                        continue
                    d = measure.srcpoldata_corr().chop()
                    snr = sw.core.snrRH(d.x,d.y) #Restivo and Helffrich (1999) signal to noise ratio
                    # print(d.x,d.y)
                    # print("splitting intensity",splitting_intensity(d))

                    ##sum the error surfaces along each of the axes, to "squash" the surface into two profiles, one for fast and one for lag
                    ## the result is best defined for the lam1/lam2 surface than the lam1 surface, or the lam2 surface
                    #- Jack Walpole
                    squashfast = np.sum(measure.lam1/measure.lam2, axis=0)
                    squashlag = np.sum(measure.lam1/measure.lam2, axis=1)

                    mean_max_lam12_fast = np.max(squashfast)/np.mean(squashfast)
                    mean_max_lam12_lag = np.max(squashlag)/np.mean(squashlag)

                    ## Null test
                    ## The measurements that fail the constrain of the maximum allowed error in delay time and the maximum delay time can be associated with null measurements because this happens due little energy on the transverse component to constrain delay time (Evans et al., 2006). 

                    diff_mult = auto_null_measure(measure,squashfast,squashlag,plot_null=False)
                    null_thresh = 0.05 #below this value, the measurement is classified as null
                    if diff_mult<null_thresh:
                        if stn_meas_close:
                            sks_meas_file_null.write("{} {:8.4f} {:8.4f} {:4.1f}\n".format(trace1[0].stats.event_time,trace1[0].stats.event_longitude,trace1[0].stats.event_latitude,trace1[0].stats.event_depth,trace1[0].stats.back_azimuth))
                        self.logger.info("{}/{} Null measurement {}".format(count,int(len(data)/3),trace1[0].stats.event_time))
                        num_null+=1
                    else:
                        if str(inpSKSdict['sks_measurement_contrains']['sel_param']) == "snr":
                            filtres = filter_pick_snr(measure,inpSKSdict,snr)
                        elif str(inpSKSdict['sks_measurement_contrains']['sel_param']) == "lam12":
                            filtres = filter_pick_lam12(measure,inpSKSdict,mean_max_lam12_fast,mean_max_lam12_lag)

                        ## 
                        if filtres:
                            num_measurements+=1
                            if stn_meas_close:
                                sks_meas_file.write("{} {:8.4f} {:8.4f} {:4.1f} {:6.1f} {:6.1f} {:.1f} {:.1f} {:.2f} {:.2f}\n".format(trace1[0].stats.event_time,trace1[0].stats.event_longitude,trace1[0].stats.event_latitude,trace1[0].stats.event_depth,trace1[0].stats.back_azimuth,measure.fast,measure.dfast,measure.lag,measure.dlag,splitting_intensity(d)))

                            if self.plot_measure_loc and bool(inpSKSdict['sks_measurement_plot']['measurement_snapshot']):
                                plot_SKS_measure(measure)
                                plt.savefig(self.plot_measure_loc+f'{plt_id}-{evyear}_{evmonth}_{evday}_{evhour}_{evminute}.png')
                                plt.close('all')  
                                self.logger.info("{}/{} Good measurement: {}; fast = {:.2f}+-{:.2f}, lag = {:.2f}+-{:.2f}".format(count,int(len(data)/3),trace1[0].stats.event_time,measure.fast,measure.dfast,measure.lag,measure.dlag))
                                

                            if int(inpSKSdict['error_plot_toggles']['error_plot_indiv']):
                                errorplot(measure,squashfast,squashlag,figname=self.plot_measure_loc+f'errorplot_{plt_id}-{evyear}_{evmonth}_{evday}_{evhour}_{evminute}.png')
                                polar_error_surface(measure,figname=self.plot_measure_loc+f'errorplot_polar_{plt_id}-{evyear}_{evmonth}_{evday}_{evhour}_{evminute}.png')

                            if int(inpSKSdict['error_plot_toggles']['error_plot_all']):
                                measure_list.append(measure)
                                squashfast_list.append(squashfast)
                                squashlag_list.append(squashlag)

                            fast_dir = measure.degs[0,np.argmax(squashfast)]

                            #to be sure the measurements are on the same half of projection
                            if fast_dir<-45 and fast_dir>-91:
                                fast_dir = fast_dir+180
                            else:
                                fast_dir = fast_dir

                            fast_dir_all.append(fast_dir)
                            lag_time_all.append(measure.lags[np.argmax(squashlag),0])
                        else:
                            self.logger.info("{}/{} Bad measurement: {}! dfast = {:.1f}, dlag = {:.1f}, snr: {:.1f}".format(count,int(len(data)/3),stream3c[0].stats.event_time,measure.dfast,measure.dlag,snr))#; Consider changing the trim window
                else:
                    self.logger.info(f"{count}/{int(len(data)/3)} Bad phase pick: {stream3c[0].stats.event_time}")
            if stn_meas_close:
                sks_meas_file.close()
                sks_meas_file_null.close()

            if bool(inpSKSdict['error_plot_toggles']['error_plot_all']) and count>0:
                errorplot_all(measure_list,squashfast_list,squashlag_list,np.array(fast_dir_all),np.array(lag_time_all),figname=self.plot_measure_loc+f'errorplot_{plt_id}.png')

            ## Splitting intensity vs backazimuth
            if bool(inpSKSdict['sks_measurement_plot']['plot_SI']):
                sks_meas_file = self.plot_measure_loc+f"{net_name}_{stn_name}_{str(inpSKSdict['filenames']['sks_meas_indiv'])}"
                outfig = self.plot_measure_loc+f"{net_name}_{stn_name}_BAZ_SI.png"
                if os.path.exists(sks_meas_file) and not os.path.exists(outfig):
                    plot_baz_si_map(sks_meas_file = sks_meas_file, outfig = outfig)
            
            if all_meas_close:
                mean_fast_dir_all = mean_angle(fast_dir_all) if len(fast_dir_all) else 0
                
                all_measurements.write("{} {} {:.4f} {:.4f} {:.2f} {:.1f} {} {}\n".format(net_name,stn_name,data[0].stats.station_longitude,data[0].stats.station_latitude,mean_fast_dir_all,np.mean(lag_time_all),num_measurements, num_null))

        f.close()
        if all_meas_close:
            all_measurements.close()


    ## plotting the measurement
    def plot_sks_map(self):
        figname = self.plot_measure_loc+'../SKS_station_Map.png'
        if not os.path.exists(figname):
            self.logger.info("##Plotting SKS map")
            sks_meas_all = pd.read_csv(self.plot_measure_loc+"../"+"sks_measurements_all.txt",delimiter="\s+")
            

            ## Segregate data based on num of measurements
            if bool(inpSKSdict['sks_measurement_plot']['segregate_measurements_options']['segregate_measurements_tofile']):
                segregate_measurements(sks_meas_all,toTextFile=True,outloc=self.plot_measure_loc+"../")

            plot_sks_station_map(sks_meas_all,figname)
            self.logger.info(f"SKS measurement figure: SKS_station_Map.png")
    
    def plot_data_nodata_map(self,sks_stations_infofile):
        figname = self.plot_measure_loc+'../data_nodata_map.png'
        all_data_df = pd.read_csv(sks_stations_infofile,delimiter='|')
        if not os.path.exists(figname):
            self.logger.info("##Plotting data-nodata map")
            sks_meas_all = pd.read_csv(self.plot_measure_loc+"../"+"sks_measurements_all.txt",delimiter="\s+")
            plot_sks_data_nodata_map(sks_meas_all,all_data_df,figname)
        
            
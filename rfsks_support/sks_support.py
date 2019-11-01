import matplotlib.pyplot as plt
import numpy as np
# plt.style.use('ggplot')
plt.style.use('seaborn')
import matplotlib.gridspec as gridspec
import tqdm
import os, glob
from rf import read_rf, IterMultipleComponents
from obspy import UTCDateTime as UTC
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from obspy.core import read
from obspy.taup import TauPyModel
from rfsks_support.rfsks_extras import plot_trigger, plot_trace, plot_SKS_measure
from obspy.signal.trigger import recursive_sta_lta,classic_sta_lta,z_detect,carl_sta_trig,delayed_sta_lta, trigger_onset
import splitwavepy as sw
import logging




## Pre-processing
def SKScalc(dataSKSfileloc,trace_loc_ENZ=None,trace_loc_RTZ=None,trigger_loc=None,plot_measure_loc=None,method = 'None'):
    logger = logging.getLogger(__name__)
    logger.info("Cut the traces around the SKS arrival")
    sksfiles = glob.glob(dataSKSfileloc+'*-sks_profile_data.h5')
    logger.info(sksfiles)
    count=0
    for i,sksfile in enumerate(sksfiles):
        data = read_rf(sksfile, 'H5')
        logger.info(f"Calculating SKS arrival times for {sksfile}")
        net_name = os.path.basename(sksfile).split("-")[0]
        stn_name = os.path.basename(sksfile).split("-")[1]
        # print("file name",net_name+stn_name)
        sks_meas_file = open(plot_measure_loc+f"{net_name}_{stn_name}_sks_measurements.txt",'w')
        sks_meas_file.write(f"Stlon: {data[0].stats.station_longitude}; Stlat: {data[0].stats.station_latitude}; Stbaz: {data[0].stats.back_azimuth}\n")
        sks_meas_file.write("EventTime EvLong EvLat FastDirection(degs) deltaFastDir(degs) LagTime(s) deltaLagTime(s)\n")
        
        for stream3c in IterMultipleComponents(data, 'onset', 3):
            count+=1
            logger.info(f"Working on {count}/{int(len(data)/3)}")

            ## check if the length of all three traces are equal
            len_tr_list=list()
            for tr in stream3c:
                len_tr_list.append(len(tr))
            if len(set(len_tr_list))!=1:
                continue


            ## filter the trace
            st = stream3c.filter('bandpass', freqmin=0.05, freqmax=0.6)
            st.detrend('linear')
            # st.taper(max_percentage=0.05, type="hann")
            sps = st[0].stats.sampling_rate
            t = st[0].stats.starttime
            ## trim the trace
            # print('st end',st[0].stats.starttime,st[0].stats.endtime)
            ev_sttime = st[0].stats.starttime
            ev_endtime = st[0].stats.endtime
            # print('t',t,t+30, t + 110)
            trace1 = st.trim(t+30, t + 110)

            ## plot the ENZ
            if trace_loc_ENZ:
                plot_trace(trace1,trace_loc_ENZ)
            
            ## Rotate to RTZ
            ## trace2[0]->BHT; trace2[1]->BHR; trace2[2]->BHZ;
            trace1.rotate('NE->RT')
            # print('all stats',dir(trace1[0].stats))
            # print()
            snr_rt = sw.core.snrRH(trace1[1].data,trace1[0].data)
            if snr_rt>3:
                pass
            else:
                continue

            
            plt_id = f"{trace1[0].stats.network}-{trace1[0].stats.station}"
            # print("latlon",dir(trace1[0].stats))
            # print(dir(trace1[0].stats.event_time))
            evyear = trace1[0].stats.event_time.year
            evmonth = trace1[0].stats.event_time.month
            evday = trace1[0].stats.event_time.day
            evhour = trace1[0].stats.event_time.hour
            evminute = trace1[0].stats.event_time.minute
            
            # ## plot all three traces RTZ
            if trace_loc_RTZ:
                plot_trace(trace1,trace_loc_RTZ)

            ######################
            #  Different picker methods. User's choice?
            ######################
            # method = 'None'
            ### operating on transverse component
            if method=="recursive_sta_lta":
                # logger.info(f"Method is {method}")
                cft = recursive_sta_lta(trace1[1].data, int(1 * sps), int(5 * sps))
                threshold = (2.5, 0.65)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                
                if trigger_loc and on_off.shape[0]==1:
                    outfile = trigger_loc+f'{plt_id}-{trace1[0].stats.event_time}-trigger.png'
                    plot_trigger(trace1[1], cft, on_off, threshold[0], threshold[1], outfile=outfile)
                    tav = int( (on_off[:, 0] /sps + on_off[:, 1] /(sps)))
                    # print('t, tav, sps =', t, tav,sps)
                    # print('trim =', t+tav-30,t+tav+30)
                    # print('edtimediff',(t+tav+30),ev_endtime,(t+tav+30) - ev_endtime)
                    # print('sttimediff',(t+tav-30),ev_sttime, (t+tav-30) - ev_sttime)

                    begttrim = t+tav-30
                    endttrim = t+tav+30

                    if (t+tav+30) - ev_endtime < 0:
                        endttrim = ev_endtime
                    
                    if (t+tav-30) - ev_sttime < 0:
                        begttrim = ev_sttime
                    trace1 = st.trim(begttrim,endttrim)

            elif method=="classic_sta_lta":
                cft = classic_sta_lta(trace1[1].data, int(5 * sps), int(10 * sps))
                threshold = (1.5, 0.5)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
            elif method=="z_detect":
                cft = z_detect(trace1[1].data, int(10 * sps))
                threshold = (-0.4, -0.3)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
            elif method=="carl_sta_trig":
                cft = carl_sta_trig(trace1[1].data, int(5 * sps), int(10 * sps), 0.8, 0.8)
                threshold = (20.0, -20.0)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
            elif method=="delayed_sta_lta":
                cft = delayed_sta_lta(trace1[1].data, int(5 * sps), int(10 * sps))
                threshold = (5, 10)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
            else:
                logger.info("No valid method specified")
                pass

            if on_off.shape[0]==1:
                # print(snr_rt,'\n')
                # t = trace1[0].stats.starttime

                
                # print(trace1)
                # trace1 = trace1.trim()
                # logger.info("------> Measure the splitting")
                # print(trace1[0].stats)
                trace1.rotate('RT->NE')
                trace2 = trace1
                logger.info(f"Measure splitting for {plt_id}-{trace1[0].stats.event_time}: {trace2[1].stats.channel},{trace2[0].stats.channel}")
                realdata = sw.Pair(trace2[1].data,trace2[0].data, delta=1/sps)
                try:
                    measure = sw.EigenM(realdata, lags=(3,))
                    
                except Exception as e:
                    logger.error(e)
                    continue
                if measure.dfast < 5 and measure.dlag < 1.5:
                    '''
                    Number of degrees of freedom is less than 3 may lead to a spurios measurement.
                    '''
                    sks_meas_file.write("{} {:8.4f} {:8.4f} {:6.1f} {:6.1f} {:.1f} {:.1f}\n".format(trace1[0].stats.event_time,trace1[0].stats.event_longitude,trace1[0].stats.event_latitude,measure.fast,measure.dfast,measure.lag,measure.dlag))
                    # sks_meas_file.write("{} {:8.4f} {:8.4f} {:6.1f} {:6.1f} {:.1f} {:.1f}\n".format(trace1[0].stats.event_time,trace1[0].stats.station_longitude,trace1[0].stats.station_latitude,measure.fast,measure.dfast,measure.lag,measure.dlag))
                    if plot_measure_loc:
                        plot_SKS_measure(measure)
                        plt.savefig(plot_measure_loc+f'{plt_id}-{evyear}_{evmonth}_{evday}_{evhour}_{evminute}.png')
                        plt.close('all')  
                else:
                    logger.warning("Measurement rejected! Consider changing the trim window")

        sks_meas_file.close()
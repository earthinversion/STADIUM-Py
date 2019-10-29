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
            st = stream3c.filter('bandpass', freqmin=0.01, freqmax=0.6)
            st.detrend('linear')
            st.taper(max_percentage=0.05, type="hann")
            sps = st[0].stats.sampling_rate
            t = st[0].stats.starttime
            ## trim the trace
            trace1 = st.trim(t, t + 70)

            ## plot the ENZ
            if trace_loc_ENZ:
                plot_trace(trace1,trace_loc_ENZ)
            
            ## Rotate to RTZ
            ## trace2[0]->BHT; trace2[1]->BHR; trace2[2]->BHZ;
            trace1.rotate('NE->RT')
            
            plt_id = f"{trace1[0].stats.network}-{trace1[0].stats.station}"
            
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
                trig_on = on_off[:,0][0]-10
                trig_off = on_off[:,1][0]+30
                t = trace1[0].stats.starttime
                # logger.info(trig_on,trig_off,UTC(trig_on),UTC(trig_off),t)
                # logger.info("------> Measure the splitting")
                trace1.rotate('RT->NE')
                trace2 = trace1
                # trace2 = trace1.trim(UTC(trig_on), UTC(trig_off))
                # logger.info('\nTrace1',trace1[0].stats.channel,trace1[1].stats.channel,trace1[2].stats.channel)
                logger.info(f"Measure splitting for {plt_id}-{trace1[0].stats.event_time}: {trace2[1].stats.channel},{trace2[0].stats.channel}")
                realdata = sw.Pair(trace2[1].data,trace2[0].data, delta=1/sps)
                try:
                    measure = sw.EigenM(realdata, lags=(3,))
                    # print(dir(measure))
                    print(measure.fast,measure.dfast,measure.lag,measure.dlag)
                    # print(measure.lam1)
                    
                except Exception as e:
                    logger.error(e)
                    continue
                
                if plot_measure_loc:
                    plot_SKS_measure(measure)
                    # plt.savefig(f"{plt_id}-{trace1[0].stats.event_time}.png")
                    plt.savefig(plot_measure_loc+f'{plt_id}-{trace1[0].stats.event_time}.png')
                    plt.close('all')              
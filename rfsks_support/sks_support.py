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
from rfsks_support.rfsks_extras import plot_trigger, plot_trace, rotateNE_RT
from obspy.signal.trigger import recursive_sta_lta,classic_sta_lta,z_detect,carl_sta_trig,delayed_sta_lta, trigger_onset
import splitwavepy as sw



## Pre-processing
def SKScalc(dataSKSfileloc,trace_loc_ENZ=None,trace_loc_RTZ=None,trigger_loc=None,plot_measure_loc=1,method = 'None'):
    print("--> Cut the traces around the SKS arrival")
    sksfiles = glob.glob(dataSKSfileloc+'*-sks_profile_data.h5')
    for i,sksfile in enumerate(sksfiles):
        data = read_rf(sksfile, 'H5')
        print(f"----> Calculating SKS arrival times for {sksfile}")
        for stream3c in tqdm.tqdm(IterMultipleComponents(data, 'onset', 3)):
            ## filter the trace
            st = stream3c.filter('bandpass', freqmin=0.01, freqmax=0.6)
            sps = st[0].stats.sampling_rate
            t = st[0].stats.starttime
            ## trim the trace
            trace1 = st.trim(t + 30, t + 90)

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
                # print(f"Method is {method}")
                cft = recursive_sta_lta(trace1[1].data, int(1 * sps), int(5 * sps))
                threshold = (2.5, 0.65)
                on_off = np.array(trigger_onset(cft, threshold[0], threshold[1]))
                
                if trigger_loc and on_off.shape[0]==1:
                    outfile = trigger_loc+f'{plt_id}-{trace1[0].stats.event_time}-trigger.png'
                    plot_trigger(trace1[0], cft, on_off, threshold[0], threshold[1], outfile=outfile)

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
                print("No valid method specified")
                pass

            if on_off.shape[0]==1:
                trig_on = on_off[:,0][0]-5
                trig_off = on_off[:,1][0]+5
                t = trace1[0].stats.starttime
                print(trig_on,trig_off,UTC(trig_on),UTC(trig_off),t)
                # print("------> Measure the splitting")
                trace1.rotate('RT->NE')
                trace2 = trace1
                # trace2 = trace1.trim(UTC(trig_on), UTC(trig_off))
                # print('\nTrace1',trace1[0].stats.channel,trace1[1].stats.channel,trace1[2].stats.channel)
                realdata = sw.Pair(trace2[1].data,trace2[0].data, delta=1/sps)
                measure = sw.EigenM(realdata)
                
                if plot_measure_loc:
                    # setup figure and subplots
                    fig = plt.figure(figsize=(12,6)) 
                    gs = gridspec.GridSpec(2, 3,width_ratios=[1,1,2])
                    ax0 = plt.subplot(gs[0,0])
                    ax1 = plt.subplot(gs[0,1])
                    ax2 = plt.subplot(gs[1,0])
                    ax3 = plt.subplot(gs[1,1])
                    ax4 = plt.subplot(gs[:,2])
                    # data to plot
                    d1 = measure.data.chop()
                    d1f = measure.srcpoldata().chop()
                    d2 = measure.data_corr().chop()
                    d2s = measure.srcpoldata_corr().chop()
                    # get axis scaling
                    lim = np.abs(d2s.data()).max() * 1.1
                    ylim = [-lim,lim]

                    # original
                    d1f._ptr(ax0,ylim=ylim)
                    d1._ppm(ax1,lims=ylim)
                    # corrected
                    d2s._ptr(ax2,ylim=ylim)
                    d2._ppm(ax3,lims=ylim)
                    kwargs={}
                    kwargs['vals'] = measure.lam1 / measure.lam2
                    kwargs['title'] = r'$\lambda_1 / \lambda_2$'
                    
                    # add marker and info box by default
                    kwargs['marker'] = True
                    kwargs['info'] = True
                    kwargs['conf95'] = True

                    measure._psurf(ax4,**kwargs)
                    plt.tight_layout()
                    plt.savefig(plot_measure_loc+f'{plt_id}-{trace1[0].stats.event_time}.png')
                    plt.close('all')
                    

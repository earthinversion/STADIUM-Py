import itertools
import numpy as np
from rf import RFStream
import tqdm
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy import UTCDateTime as UTC
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# plt.style.use('ggplot')
plt.style.use('seaborn')
import logging
from rfsks_support.other_support import Timeout
import matplotlib.gridspec as gridspec




def minendtime(alledtimes):
        minval0 = alledtimes[0]
        for val in alledtimes[1:]:
            if val<minval0:
                minval0 = val
        return minval0

def maxstarttime(allsttimes):
    maxval0 = allsttimes[0]
    for val in allsttimes[1:]:
        if val>maxval0:
            maxval0 = val
    return maxval0

def filter_traces_rf(stream,pharr = None):
    logger = logging.getLogger(__name__)
    #Filtering traces
    
    allsttimes = []
    alledtimes = []
    for tr in stream:
        sttime = tr.stats.starttime
        edtime = tr.stats.endtime
        allsttimes.append(sttime)
        alledtimes.append(edtime)
        
    
    if minendtime(alledtimes)-pharr>=75 and pharr - maxstarttime(allsttimes)>=25:
        
        stream = stream.trim(pharr - 25, pharr+ 75)
        for tr in stream:
        
            if tr.stats.sampling_rate < 20:
                logger.warning(f"Sampling rate too low: {tr.stats.sampling_rate}, required >= 20Hz")
                stream.remove(tr)
                continue
            elif tr.stats.sampling_rate >= 20:
                if tr.stats.sampling_rate % 20 == 0:
                    factor = int(tr.stats.sampling_rate / 20)
                    # logger.warning(f"Downsampling to 20 Hz, current sr: {tr.stats.sampling_rate}, factor: {factor}")
                    tr.decimate(factor, strict_length=False, no_filter=True) 
                    continue 
                    # logger.warning(f"After Downsampling to 20 Hz, current sr: {tr.stats.sampling_rate}")
                else:
                    tr.resample(20.0)
                    logger.warning(f"Resampling traces; New sampling rate: {tr.stats.sampling_rate}")
                    # stream.remove(tr)
                    continue
            else:
                pass
        
        

DEG2KM = 111.2  #: Conversion factor from degrees epicentral distance to km

def filter_traces_sks(stream,pharr = None):
    logger = logging.getLogger(__name__)
    #Filtering traces
    # obtaining the start and end time for all traces in stream
    allsttimes = []
    alledtimes = []
    for tr in stream:
        sttime = tr.stats.starttime
        edtime = tr.stats.endtime
        allsttimes.append(sttime)
        alledtimes.append(edtime)
    
    # filtering for defined window
    if minendtime(alledtimes)-pharr>=60 and pharr - maxstarttime(allsttimes)>=60:
        mincut = np.min([minendtime(alledtimes)-pharr,pharr - maxstarttime(allsttimes)]) #min dist from arrival time
        mincut = 60 if mincut > 60 else None
        
        stream = stream.trim(pharr - mincut, pharr+ mincut)
        for tr in stream:            
            if tr.stats.sampling_rate < 20:
                logger.warning(f"Sampling rate too low: {tr.stats.sampling_rate}, required >= 20Hz")
                stream.remove(tr)
                continue
            elif tr.stats.sampling_rate >= 20:
                if tr.stats.sampling_rate % 20 == 0:
                    factor = int(tr.stats.sampling_rate / 20)
                    # logger.warning(f"Downsampling to 20 Hz, current sr: {tr.stats.sampling_rate}, factor: {factor}")
                    tr.decimate(factor, strict_length=False, no_filter=True) 
                    continue 
                    # logger.warning(f"After Downsampling to 20 Hz, current sr: {tr.stats.sampling_rate}")
                else:
                    tr.resample(20.0)
                    logger.warning(f"Resampling traces; New sampling rate: {tr.stats.sampling_rate}")
                    # stream.remove(tr)
                    continue
            else:
                pass
        
        


def _get_stations(inventory):
    channels = inventory.get_contents()['channels']
    stations = {ch[:-1] + '?': ch[-1] for ch in channels}
    return stations


def iter_event_data(events, inventory, get_waveforms, phase='P',
                    request_window=None, pad=10, pbar=None, **kwargs):
    """
    Return iterator yielding three component streams per station and event.

    :param events: list of events or `~obspy.core.event.Catalog` instance
    :param inventory: `~obspy.core.inventory.inventory.Inventory` instance
        with station and channel information
    :param get_waveforms: Function returning the data. It has to take the
        arguments network, station, location, channel, starttime, endtime.
    :param phase: Considered phase, e.g. 'P', 'S', 'PP'
    :type request_window: tuple (start, end)
    :param request_window: requested time window around the onset of the phase
    :param float pad: add specified time in seconds to request window and
       trim afterwards again
    :param pbar: tqdm_ instance for displaying a progressbar
    :param kwargs: all other kwargs are passed to `~rf.rfstream.rfstats()`

    :return: three component streams with raw data

    Example usage with progressbar::

        from tqdm import tqdm
        from rf.util import iter_event_data
        with tqdm() as t:
            for stream3c in iter_event_data(*args, pbar=t):
                do_something(stream3c)

    .. _tqdm: https://pypi.python.org/pypi/tqdm
    """
    logger = logging.getLogger(__name__)
    from rf.rfstream import rfstats, RFStream
    method = phase[-1].upper()
    if request_window is None:
        request_window = (-50, 150) if method == 'P' else (-100, 50)
    # else:
        # print(f"\nRequested window is {request_window}")
    stations = _get_stations(inventory)
    # print(f"Available stations are {stations}")
    if pbar is not None:
        pbar.total = len(events) * len(stations)
    for event, seedid in itertools.product(events, stations):
        if pbar is not None:
            pbar.update(1)
        origin_time = (event.preferred_origin() or event.origins[0])['time']
        # print(f"Origin time is {origin_time}")
        try:
            args = (seedid[:-1] + stations[seedid], origin_time)
            coords = inventory.get_coordinates(*args)
            # print(f"Station coords are {coords}")
        except Exception as exception:  # station not available at that time
            logger.error("Station not available at given time", exc_info=True)
            continue
        try:
            stats = rfstats(station=coords, event=event, phase=phase, **kwargs)
        except:
            logger.warning(f"Unable to read for: {event}")
            pass
            
        if not stats:
            continue
        net, sta, loc, cha = seedid.split('.')
        starttime = stats.onset + request_window[0]
        endtime = stats.onset + request_window[1]
        kws = {'network': net, 'station': sta, 'location': loc,
               'channel': cha, 'starttime': starttime - pad,
               'endtime': endtime + pad}
        try:
            stream = get_waveforms(**kws)
            # print(f"stream obtained {len(stream)}")
        except:  # no data available
            # logger.warning(f"No data available for {event}")
            continue
        stream.trim(starttime, endtime)
        stream.merge()
        # print(f"length of the stream is {len(stream)}")
        if len(stream) != 3:
            from warnings import warn
            warn(f'Need 3 component seismograms. {len(stream)} components '
                 'detected for event {event.resource_id}, station {seedid}.')
            logger.warning(f'Need 3 component seismograms. {len(stream)} components '
                 'detected for event {event.resource_id}, station {seedid}.')
            continue
        if any(isinstance(tr.data, np.ma.masked_array) for tr in stream):
            from warnings import warn
            logger.warning(f'Gaps or overlaps detected for event {event.resource_id}, station {seedid}.')
            continue
        for tr in stream:
            tr.stats.update(stats)
        yield RFStream(stream)

def retrieve_waveform(client,net,stn,t1,t2,stats_dict=None,cha="BHE,BHN,BHZ",attach_response=False,loc="",pharr=None, phasenm = 'P'):  
    try:  
        st = client.get_waveforms(net, stn, loc, cha, t1, t2,attach_response=attach_response)
    except:
        return False
    # print("Retrieving")
    if phasenm == 'P':
        # filter_traces(st,lenphase=int(t2-t1))
        filter_traces_rf(st,pharr = pharr)
    elif phasenm == 'SKS':
        filter_traces_sks(st,pharr = pharr)

    if len(st) != 3:
        # print(f"All three components not available: {len(st)}")
        return False
    if stats_dict:
        dist, baz, _ = gps2dist_azimuth(stats_dict['station_latitude'],stats_dict['station_longitude'],stats_dict['event_latitude'],stats_dict['event_longitude'])
        for tr in st:
            for key, value in stats_dict.items():
                tr.stats[key] = value
            tr.stats['distance'] = dist / 1000 / DEG2KM
            tr.stats['back_azimuth'] = baz
    
    st.merge()
    # if all 3 components present and no gap or overlap in data
    if len(st) == 3 and not any(isinstance(tr.data, np.ma.masked_array) for tr in st):
        # print(f"Stream obtained {len(st)}")
        return RFStream(st)
    elif not any(isinstance(tr.data, np.ma.masked_array) for tr in st):
        # print("--------> There's a gap/overlap in the data")
        return False

def multi_download(client,inv,net,stn,slat,slon,elat,elon,evdp,evtime,em,emt,fcat,stalons,stalats,staNetNames,phase='P',locations=[""]):
    logger = logging.getLogger(__name__)
    strm = None
    j=0
    msg = None
    model = TauPyModel('iasp91')
    arrivals = model.get_travel_times_geo(float(evdp),slat,slon,float(elat),float(elon),phase_list=[phase])
    if phase=='P':
        t1 = UTC(str(evtime)) + int(arrivals[0].time - 50)
        t2 = UTC(str(evtime)) + int(arrivals[0].time + 110)
        # t1 = UTC(str(evtime)) + int(arrivals[0].time - 25)
        # t2 = UTC(str(evtime)) + int(arrivals[0].time + 75)
    elif phase=='SKS':
        t1 = UTC(str(evtime)) + int(arrivals[0].time - 80)
        t2 = UTC(str(evtime)) + int(arrivals[0].time + 80)
    # sel_inv = inv.select(network=net).select(station=stn)[0][0]
    # if not sel_inv.is_active(starttime=t1, endtime=t2):
    #     # logger.warning(f"------> Station not active during {evtime}")
    #     msg = f"Station not active during {evtime}"
    #     return strm, 0, msg
    # process_id = os.getpid()
    pharr = UTC(str(evtime)) + arrivals[0].time
    while not strm:
        client_local = Client(client[j])
        stats_args = {"_format":'H5', "onset" : UTC(str(evtime)) + arrivals[0].time, "event_latitude": elat, "event_longitude": elon,"event_depth":evdp, "event_magnitude":em,"event_time":UTC(str(evtime)),"phase":phase,"station_latitude":slat,"station_longitude":slon,"inclination":arrivals[0].incident_angle,"slowness":arrivals[0].ray_param_sec_degree}
        if phase=='P':
            for loc in locations:
                with Timeout(5):
                    strm = retrieve_waveform(client_local,net,stn,t1,t2,stats_dict=stats_args,cha="BHE,BHN,BHZ",loc=loc,pharr = pharr, phasenm = phase)
                    if strm:
                        break
                
        elif phase=='SKS':
            for loc in locations:
                # print(f"Location: {loc}")
                with Timeout(5):
                    
                    strm = retrieve_waveform(client_local,net,stn,t1,t2,stats_dict=stats_args,cha="BHE,BHN,BHZ",attach_response=True,loc=loc,pharr = pharr, phasenm = phase)
                    if strm:
                        break #break the locations loop
                
        if strm:
            fcat.write('{} | {:9.4f}, {:9.4f} | {:5.1f} | {:5.1f} {:4s} | {}\n'.format(evtime,elat,elon,evdp,em,emt,client[j]))
            stalons.append(slon)
            stalats.append(slat)
            staNetNames.append(f"{net}_{stn}")
            # print("stream obtained\n")
            msg = f"Data {evtime}"
            res = 1
            break
        elif j == len(client)-1:
            res = 0
            msg = f"No data {evtime}"
            break
        j+=1
    return strm, res, msg


def plot_trigger(trace, cft, on_off, thr_on, thr_off,outfile):
    """
    Plot characteristic function of trigger along with waveform data and
    trigger On/Off from given thresholds.

    :type trace: :class:`~obspy.core.trace.Trace`
    :param trace: waveform data
    :type cft: :class:`numpy.ndarray`
    :param cft: characteristic function as returned by a trigger in
        :mod:`obspy.signal.trigger`
    :type thr_on: float
    :param thr_on: threshold for switching trigger on
    :type thr_off: float
    :param thr_off: threshold for switching trigger off
    :type show: bool
    :param show: Do not call `plt.show()` at end of routine. That way,
        further modifications can be done to the figure before showing it.
    """
    logger = logging.getLogger(__name__)
    from obspy.signal.trigger import trigger_onset

    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / df
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, trace.data, 'k',lw=0.5)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.plot(t, cft, 'k',lw=0.5)
    # on_off = np.array(trigger_onset(cft, thr_on, thr_off))
    i, j = ax1.get_ylim()
    try:
        ax1.vlines(on_off[:, 0] / df, i, j, color='r', lw=2,
                   label="Trigger On")
        ax1.vlines(on_off[:, 1] / df, i, j, color='b', lw=2,
                   label="Trigger Off")
        ax1.legend()
    except IndexError:
        logger.error('IndexError', exc_info=True)
        pass
    ax2.axhline(thr_on, color='red', lw=1, ls='--')
    ax2.axhline(thr_off, color='blue', lw=1, ls='--')
    ax2.set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())
    fig.suptitle(trace.id)
    fig.canvas.draw()
    plt.savefig(outfile,dpi=150,bbox_inches='tight')
    plt.close('all')
    # return on_off

def plot_trace(trace2,trace_loc):
    logger = logging.getLogger(__name__)
    try:
        fig, ax = plt.subplots(3,1,figsize=(10,6),sharex=True)
        ax[0].plot(trace2[0].times("matplotlib"), trace2[0].data, "r-",lw=0.5,label=trace2[0].stats.channel)
        ax[0].legend(loc='best')
        ax[1].plot(trace2[1].times("matplotlib"), trace2[1].data, "g-",lw=0.5,label=trace2[1].stats.channel)
        ax[1].legend(loc='best')
        ax[2].plot(trace2[2].times("matplotlib"), trace2[2].data, "b-",lw=0.5,label=trace2[2].stats.channel)
        ax[2].legend(loc='best')
        ax[2].xaxis_date()
        fig.autofmt_xdate()
        plt_id = f"{trace2[0].stats.network}-{trace2[0].stats.station}"
        ax[0].set_title(plt_id+f" Onset: {trace2[0].stats.onset}")
        
        plt.savefig(trace_loc+f'{plt_id}-{trace2[0].stats.event_time}-RTZ.png',dpi=150,bbox_inches='tight')
        plt.close('all')
    except Exception as exception:
        logger.error("Plot trace error", exc_info=True)
def plot_SKS_measure(measure):
    fig = plt.figure(figsize=(12,6)) 
    gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,2])    
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[:,2])
    d1 = measure.data
    # print(measure.data._w0(), measure.data._w1())
    # print(dir(measure))
    # print(d1.t())
    d1f = measure.srcpoldata()
    d2 = measure.data_corr()
    d2s = measure.srcpoldata_corr()
    # # get axis scaling
    lim = np.abs(d2s.data()).max() * 1.1
    ylim = [-lim,lim]
    # print(d1f.wbeg(), d1f.wend())
    # print(d1f._centresamp(), d1f.window.offset, d1f.window.width/2)
    d1f._ptr(ax0,ylim=ylim, cmplabels=['N','E'])
    d1._ppm(ax1,lims=ylim, cmplabels=['',''])
    # # corrected
    d2s._ptr(ax2,ylim=ylim, cmplabels=['R','T'])
    d2._ppm(ax3,lims=ylim, cmplabels=['',''])
    # add marker and info box by default
    kwargs={}
    kwargs['marker'] = True
    kwargs['info'] = True
    kwargs['conf95'] = True
    kwargs['cmap'] = 'rainbow'
    measure._psurf(ax4,**kwargs)
    plt.tight_layout()


def filter_pick_snr(measure,inpSKSdict,snr):
    if measure.dfast < int(inpSKSdict['sks_measurement_contrains']['fast_dir_settings']['maxdfast']) and measure.dlag < float(inpSKSdict['sks_measurement_contrains']['lag_settings']['maxdlag','VALUES']) and snr > float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['snr_ratio']):
        '''
        Uses the one sigma error in fast direction and lag time. Calculated by taking a quarter of the width of 95% confidence region (found using F-test) of lambda2. And signal to noise ratio of the trace
        '''
        return True
    else:
        return False

def filter_pick_lam12(measure,inpSKSdict,mean_max_lam12_fast,mean_max_lam12_lag):
    if measure.dfast < int(inpSKSdict['sks_measurement_contrains']['fast_dir_settings']['maxdfast']) and measure.dlag < float(inpSKSdict['sks_measurement_contrains']['lag_settings']['maxdlag']) and mean_max_lam12_fast > float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['lam12fast_threh']) and mean_max_lam12_lag > float(inpSKSdict['sks_measurement_contrains']['sel_param_settings']['lam12lag_threh']):
        '''
        Uses the one sigma error in fast direction and lag time. Calculated by taking a quarter of the width of 95% confidence region (found using F-test) of lambda2. And ratio of the lambda1/lambda2 for fast direction and lag time. Empirically found as more robust than snr.
        '''
        return True
    else:
        return False

def errorplot(measure,squashfast,squashlag,figname):
    plt.close('all')
    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(measure.degs[0,:],measure.fastprofile(),'b')
    ax[0,0].axvline(measure.fast,color='r')
    ax[0,0].axvline(measure.fast-2*measure.dfast,alpha=0.5,color='r')
    ax[0,0].axvline(measure.fast+2*measure.dfast,alpha=0.5,color='r')
    ax[0,0].set_title('fast direction')

    ax[0,1].plot(measure.lags[:,0],measure.lagprofile(),'b')
    ax[0,1].axvline(measure.lag,color='r')
    ax[0,1].axvline(measure.lag-2*measure.dlag,alpha=0.5,color='r')
    ax[0,1].axvline(measure.lag+2*measure.dlag,alpha=0.5,color='r')
    ax[0,1].set_title('lag time')

    ax[1,0].plot(measure.degs[0,:],squashfast)
    ax[1,0].axvline(x=measure.degs[0,np.argmax(squashfast)],color='r')
    ax[1,0].set_title(f'L1/L2 Fast: {measure.degs[0,np.argmax(squashfast)]}')

    ax[1,1].plot(measure.lags[:,0],squashlag)
    ax[1,1].axvline(x=measure.lags[np.argmax(squashlag),0],color='r')
    ax[1,1].set_title(f'L1/L2 Lag: {measure.lags[np.argmax(squashlag),0]}')
    plt.savefig(figname,bbox_inches='tight')

def errorplot_all(measure_list,squashfast_list,squashlag_list,fast_dir_all,lag_time_all,figname):
    plt.close('all')
    fig,ax = plt.subplots(1,2,figsize=(10,6))
    for meas,sqfast,sqlag in zip(measure_list,squashfast_list,squashlag_list):
        ax[0].plot(meas.degs[0,:],sqfast,lw=0.5)
        ax[1].plot(meas.lags[:,0],sqlag,lw=0.5)
    
    #plot avg fast dir and lag time
    # fast_dir_pos = fast_dir_all[fast_dir_all>0]
    # ax[0].axvline(x=np.mean(fast_dir_pos),color='r',lw=1)
    # ax[0].axvline(x=np.mean(fast_dir_pos)+np.std(fast_dir_pos),color='r',alpha=0.5,lw=1)
    # ax[0].axvline(x=np.mean(fast_dir_pos)-np.std(fast_dir_pos),color='r',alpha=0.5,lw=1)
    # ax[0].set_title('Avg. fast direction (abs): {:.2f}'.format(np.mean(fast_dir_pos)))

    # fast_dir_neg = fast_dir_all[fast_dir_all<0]
    # ax[0].axvline(x=-np.mean(fast_dir_neg),color='r',lw=1)
    # ax[0].axvline(x=-np.mean(fast_dir_neg)+np.std(fast_dir_neg),color='r',alpha=0.5,lw=1)
    # ax[0].axvline(x=-np.mean(fast_dir_neg)-np.std(fast_dir_neg),color='r',alpha=0.5,lw=1)

    ax[1].axvline(x=np.mean(lag_time_all),color='r',lw=1)
    min_lag_err = np.mean(lag_time_all)-np.std(lag_time_all)
    if min_lag_err<0:
        min_lag_err = 0
    ax[1].axvline(x=np.mean(lag_time_all)+np.std(lag_time_all),color='r',alpha=0.5,lw=1)
    ax[1].axvline(x=min_lag_err,color='r',alpha=0.5,lw=1)
    ax[1].set_title('Avg. lag time: {:.2f}'.format(np.mean(lag_time_all)))

    plt.savefig(figname,bbox_inches='tight',dpi=300)

def null_intensity(diff,mult):
    sumdiffsq = np.sum(diff**2) #should be maximized for null measurement
    summult = np.sum(mult) #should be minimized
    return sumdiffsq/summult

def auto_null_measure(measure,squashfast,squashlag,plot_null=False):
    '''
    - Jack Walpole (Modified by Utpal Kumar)
    We can use the property of self-similarity in null error surfaces at 90 degrees rotation to automatically detect null measurements.
    The approach developed here is to shift the squashed fast profile by 90 degrees. If this is a perfect null measurement this shifted trace will be very similar to the original trace. The result of multiplying, sample by sample, the two traces will be maximised; the result of subtracting one trace from the other will be minimised. We call these two traces mult and diff. The ratio diff / mult will be small and tend to zero for null measurements.
    '''

    diff = squashfast - np.roll(squashfast,45)
    mult = squashfast * np.roll(squashfast,45)

    ratio_diff_mult = null_intensity(diff,mult)    
    return ratio_diff_mult

    if plot_null:
        fig, ax = plt.subplots(figsize=(10,6))
        ax.plot(measure.degs[0,:],squashfast)
        ax.plot(measure.degs[0,:],np.roll(squashfast,45))

def polar_error_surface(measure,figname):
    plt.close('all')
    rads = np.deg2rad(np.column_stack((measure.degs,measure.degs+180,measure.degs[:,0]+360)))
    lags = np.column_stack((measure.lags,measure.lags,measure.lags[:,0]))
    lam1 = np.column_stack((measure.lam1,measure.lam1,measure.lam1[:,0]))
    lam2 = np.column_stack((measure.lam2,measure.lam2,measure.lam2[:,0]))                  


    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'),nrows=1,ncols=3,figsize=(20,8))
    ax[0].contourf(rads,lags,lam1,50,cmap='magma')
    ax[0].set_title('Lam1')
    ax[1].contourf(rads,lags,lam2,50,cmap='magma')
    ax[1].set_title('Lam2')
    ax[2].contourf(rads,lags,lam1/lam2,50,cmap='viridis')
    ax[2].set_title('Lam1/Lam2')
    plt.savefig(figname,dpi=300,bbox_inches='tight')

def splitting_intensity(d):
    '''
    Calculate splitting intensity.
    - Utpal Kumar (modified after Walpole et al.)
    #Chevrot (2000)
    '''
    rdiff = np.gradient(d.x)
    trans = d.y
    si = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
    return si

## Fine tuning of SKS
import yaml
with open('Settings/advSKSparam.yaml') as f:
    inpSKSdict = yaml.load(f, Loader=yaml.FullLoader)

def segregate_measurements(sks_meas_all,toTextFile=False,outloc="./"):
    plot_params_lev = inpSKSdict['sks_measurement_plot']['meas_seg_points']
    lev1, lev2, lev3 = int(plot_params_lev['lev1']), int(plot_params_lev['lev2']), int(plot_params_lev['lev3'])
    station_data_0 = sks_meas_all.loc[sks_meas_all['NumMeasurements']==lev1]
    station_data_14 = sks_meas_all.loc[(sks_meas_all['NumMeasurements']>=lev1+1) & (sks_meas_all['NumMeasurements']<=lev2)]
    station_data_4_11 = sks_meas_all.loc[(sks_meas_all['NumMeasurements']>=lev2+1) & (sks_meas_all['NumMeasurements']<=lev3-1)]
    station_data_15 = sks_meas_all.loc[sks_meas_all['NumMeasurements']>=lev3]
    if toTextFile:
        station_data_0.to_csv(outloc+f"{lev1}_sks_measure.txt",index=None, header=True,sep=' ', float_format='%.4f')
        station_data_14.to_csv(outloc+f"{lev2}_sks_measure.txt",index=None, header=True,sep=' ', float_format='%.4f')
        station_data_4_11.to_csv(outloc+f"{lev2}_{lev3}_sks_measure.txt",index=None, header=True,sep=' ', float_format='%.4f')
        station_data_15.to_csv(outloc+f"{lev3}_sks_measure.txt",index=None, header=True,sep=' ', float_format='%.4f')

    return station_data_0, station_data_14, station_data_4_11, station_data_15


    
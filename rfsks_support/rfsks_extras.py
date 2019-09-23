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
# plt.style.use('ggplot')
plt.style.use('seaborn')





DEG2KM = 111.2  #: Conversion factor from degrees epicentral distance to km


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
        except:  # station not available at that time
            # print("Station not available at given time")
            continue
        try:
            stats = rfstats(station=coords, event=event, phase=phase, **kwargs)
        except:
            pass
            # print(f"Error: {event}")
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
        except:  # no data available
            # print(f"No data available for {event}")
            continue
        stream.trim(starttime, endtime)
        stream.merge()
        # print(f"length of the stream is {len(stream)}")
        if len(stream) != 3:
            from warnings import warn
            warn(f'Need 3 component seismograms. {len(stream)} components '
                 'detected for event {event.resource_id}, station {seedid}.')
            continue
        if any(isinstance(tr.data, np.ma.masked_array) for tr in stream):
            from warnings import warn
            warn(f'Gaps or overlaps detected for event {event.resource_id}, station {seedid}.')
            continue
        for tr in stream:
            tr.stats.update(stats)
        yield RFStream(stream)

def retrieve_waveform(client,net,stn,t1,t2,stats_dict=None,cha="BHE,BHN,BHZ",attach_response=False,loc=""):    
    st = client.get_waveforms(net, stn, loc, cha, t1, t2,attach_response=attach_response)
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
        return RFStream(st)
    elif len(st) != 3:
        # print("--------> All requested channels not present!")
        return False
    elif not any(isinstance(tr.data, np.ma.masked_array) for tr in st):
        # print("--------> There's a gap/overlap in the data")
        return False

def multi_download(client,inv,net,stn,slat,slon,elat,elon,evdp,evtime,em,emt,fcat,stalons,stalats,staNetNames,phase='P',locations=[""]):
    strm = None
    j=0
    
    model = TauPyModel('iasp91')
    arrivals = model.get_travel_times_geo(float(evdp),slat,slon,float(elat),float(elon),phase_list=[phase])
    if phase=='P':
        t1 = UTC(str(evtime)) + int(arrivals[0].time - 25)
        t2 = UTC(str(evtime)) + int(arrivals[0].time + 75)
    elif phase=='SKS':
        t1 = UTC(str(evtime)) + int(arrivals[0].time - 60)
        t2 = UTC(str(evtime)) + int(arrivals[0].time + 60)
    sel_inv = inv.select(network=net).select(station=stn)[0][0]
    if not sel_inv.is_active(starttime=t1, endtime=t2):
        print(f"------> Station not active during {evtime}")
        strm, 0
    # process_id = os.getpid()
    while not strm:
        client_local = Client(client[j])
        stats_args = {"_format":'H5', "onset" : UTC(str(evtime)) + arrivals[0].time, "event_latitude": elat, "event_longitude": elon,"event_depth":evdp, "event_magnitude":em,"event_time":UTC(str(evtime)),"phase":phase,"station_latitude":slat,"station_longitude":slon,"inclination":arrivals[0].incident_angle,"slowness":arrivals[0].ray_param_sec_degree}
        if phase=='P':
            for loc in locations:
                # print(f"loc is {loc}.")
                try:
                    strm = retrieve_waveform(client_local,net,stn,t1,t2,stats_dict=stats_args,cha="BHE,BHN,BHZ",loc=loc)
                    if strm:
                        break
                except:
                    pass
        elif phase=='SKS':
            for loc in locations:
                try:
                    strm = retrieve_waveform(client_local,net,stn,t1,t2,stats_dict=stats_args,cha="BHE,BHN,BHZ",attach_response=True,loc=loc)
                    if strm:
                        break
                except:
                    pass
        if strm:
            fcat.write('{} | {:9.4f}, {:9.4f} | {:5.1f} | {:5.1f} {:4s} | {}\n'.format(evtime,elat,elon,evdp,em,emt,client[j]))
            stalons.append(slon)
            stalats.append(slat)
            staNetNames.append(f"{net}_{stn}")
            res = 1
            break
        elif j == len(client)-1:
            res = 0
            break
        j+=1
    return strm, res


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



def rotateNE_RT(trace,back_azimuth=None):
    import obspy.signal.rotate as rotate
    north_data = trace[1].data
    east_data = trace[0].data
    # Figure out back-azimuth.
    if back_azimuth is None:
        try:
            back_azimuth = trace[1].stats.back_azimuth
        except Exception:
            msg = "No back-azimuth specified."
            raise TypeError(msg)

        # Do one of the two-component rotations.
        input_1 = trace.select(component='N')
        input_2 = trace.select(component='E')
        for i_1, i_2 in zip(input_1, input_2):
            dt = 0.5 * i_1.stats.delta
            if (len(i_1) != len(i_2)) or \
                    (abs(i_1.stats.starttime - i_2.stats.starttime) > dt) \
                    or (i_1.stats.sampling_rate !=
                        i_2.stats.sampling_rate):
                msg = "All components need to have the same time span."
                raise ValueError(msg)
        for i_1, i_2 in zip(input_1, input_2):
            output_1, output_2 = rotate.rotate_ne_rt(i_1.data, i_2.data, back_azimuth)
            i_1.data = output_1
            i_2.data = output_2
            # Rename the components.
            i_1.stats.channel = i_1.stats.channel[:-1] + \
                'R'
            i_2.stats.channel = i_2.stats.channel[:-1] + \
                'T'
            # Add the azimuth and inclination to the stats object.
            for comp in (i_1, i_2):
                comp.stats.back_azimuth = back_azimuth
    return input_1,input_2
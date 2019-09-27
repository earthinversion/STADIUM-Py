import matplotlib.pyplot as plt
import tqdm, sys
import numpy as np
import os, glob
from rf import RFStream, read_rf, IterMultipleComponents, get_profile_boxes
from rfsks_support.plotting_map import plot_merc, plot_bm_azimuth
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.other_support import avg
from rfsks_support.profile import profile
from rfsks_support.rfsks_extras import filter_traces
import logging
from numba import jit



### Compute RF
def compute_rf(dataRFfileloc):
    logger = logging.getLogger(__name__)
    all_rfdatafile = glob.glob(dataRFfileloc+'*-rf_profile_data.h5')
    for rfdatafile in all_rfdatafile:
        network = rfdatafile.split("-")[0]
        station = rfdatafile.split("-")[1]
        rffile = f'{network}-{station}-rf_profile_rfs.h5'
        datatmp = read_rf(rfdatafile, 'H5')
        if not os.path.exists(rffile):
            logger.info(f"--> Computing RF for {rfdatafile}")
            data = read_rf(rfdatafile, 'H5')
            stream = RFStream()
            for stream3c in tqdm.tqdm(IterMultipleComponents(data, 'onset', 3)):
                if len(stream3c) != 3:
                    continue
                
                ## check if the length of all three traces are equal
                for tr in stream3c:
                    lentr=tr.stats.npts
                    lengt= tr.stats.sampling_rate * 100
                    if lentr != lengt:
                        print('Wrong trace length ', lentr,lengt)
                        continue
                
                stream3c.filter('bandpass', freqmin=0.5, freqmax=2)

                try:
                    stream3c.rf()
                except Exception as e:
                    logger.error("Problem applying rf method", exc_info=True)
                stream3c.moveout()
                stream.extend(stream3c)
            stream.write(rffile, 'H5')
        else:
            logger.info(f"--> {rffile} already exists!")

def plot_RF(dataRFfileloc,destImg,fig_frmt="png"):
    logger = logging.getLogger(__name__)
    logger.info("--> Plotting the receiver functions")
    rffiles = glob.glob(dataRFfileloc+'*-rf_profile_rfs.h5')
    for i,rffile in enumerate(rffiles):
        stream = read_rf(rffile, 'H5')
    
        kw = {'trim': (-5, 20), 'fillcolors': ('black', 'gray'), 'trace_height': 0.1}
        num_trace=len(stream.select(component='L', station=stream[0].stats.station).sort(['back_azimuth']))
        if num_trace > 0:
            try:
                stream.select(component='L', station=stream[0].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                plt.savefig(destImg + f"{stream[0].stats.station}"+'_L.'+fig_frmt)
                stream.select(component='Q', station=stream[0].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                plt.savefig(destImg + f"{stream[0].stats.station}"+'_Q.'+fig_frmt)
                logger.info("----> Working on {}/{}, {}-{} total traces: {}".format(i+1,len(rffiles),stream[0].stats.network, stream[0].stats.station, num_trace))
            except Exception as e:
                logger.error("Unexpected error", exc_info=True)
        else:
            logger.info("----> {} traces for {}-{}".format(num_trace,stream[0].stats.network, stream[0].stats.station))

def plot_pp_profile_map(dataRFfileloc,profilefileloc,catalogtxtloc,topo=True,destination="./",depth=70,fig_frmt="png",ndivlat = 2, ndivlon=3):
    logger = logging.getLogger(__name__)
    ppoints_df = pd.DataFrame()
    list_of_dfs = []
    list_of_streams = []

    logger.info("--> Plotting the piercing points on map")

    stlons,stlats=[],[]
        
    rffiles = glob.glob(dataRFfileloc+'*-rf_profile_rfs.h5')
    stream = read_rf(rffiles[0], 'H5')
    filter_traces(stream,lenphase=100)
    ppoints_tmp = stream.ppoints(depth)
    ppdf_tmp = pd.DataFrame({"pplon":ppoints_tmp[:,1],"pplat":ppoints_tmp[:,0]})
    list_of_dfs.append(ppdf_tmp)

    stlons.append(stream[0].stats.station_longitude)
    stlats.append(stream[0].stats.station_latitude)

    
    for rffile in rffiles[1:]:
        logger.info(f"----> reading profile {rffile}")
        try:
            st_tmp = read_rf(rffile, 'H5')
            filter_traces(st_tmp,lenphase=100)
            list_of_streams.append(st_tmp)
            ppoints_tmp = st_tmp.ppoints(depth)
            ppdf_tmp = pd.DataFrame({"pplon":ppoints_tmp[:,1],"pplat":ppoints_tmp[:,0]})
            list_of_dfs.append(ppdf_tmp)
            stlons.append(st_tmp[0].stats.station_longitude)
            stlats.append(st_tmp[0].stats.station_latitude)
        except:
            logger.error("Error", exc_info=True)
    logger.info("Calculating profile")
    ppoints_df = ppoints_df.append(list_of_dfs , ignore_index=True)


    for st in list_of_streams:
        for tr in st:
            stream.append(tr)


    abs_latvals = np.absolute(ppoints_df["pplat"].values)
    abs_lonvals = np.absolute(ppoints_df["pplon"].values)
    degkmfac = 111.2
    width_lat = (np.amax(abs_latvals)-np.amin(abs_latvals)) * degkmfac
    # width_lat +=0.2*width_lat
    width_lon = (np.amax(abs_lonvals)-np.amin(abs_lonvals)) * degkmfac
    # width_lon +=0.2*width_lon
    # print(width_lat,width_lon)
        
    ndivlat += 1
    ndivlon += 1
    stlat = ppoints_df["pplat"].min()
    stlon = ppoints_df["pplon"].min()
    for azimuth in [0,90]:
        if azimuth == 90:
            mxbin = width_lat
            ndiv=ndivlat
            divisions = np.linspace(stlon,stlon+mxbin/degkmfac,ndivlon)
        elif azimuth == 0:
            mxbin = width_lon
            ndiv=ndivlon
            divisions = np.linspace(stlat,stlat+mxbin/degkmfac,ndivlat)

        for n in range(len(divisions)-1):
            initdiv = divisions[n]
            enddiv = divisions[n+1]
            widthprof = int(np.abs(enddiv-initdiv)*degkmfac)
            # print(initdiv,enddiv,widthprof)
            # print(initdiv,enddiv)
            outputfile = profilefileloc+ f"rf_profile_profile{azimuth}_{int(initdiv)}_{int(enddiv)}_{widthprof}_{n}.h5"
            if not os.path.exists(outputfile):
                # logger.info("Calculating the boxes")
                if azimuth == 90:
                    logger.info(f'For az: {azimuth} startlat: {stlat}, startlon: {initdiv}, endlon: {enddiv}, length: {widthprof}')
                    boxes = get_profile_boxes((stlat, initdiv), azimuth, np.linspace(0, widthprof, int((widthprof)/5)), width=mxbin)
                elif azimuth == 0:
                    logger.info(f'For az: {azimuth} startlat: {initdiv}, startlon: {stlon}, endlat: {enddiv}, length: {widthprof}')
                    boxes = get_profile_boxes((initdiv, stlon), azimuth, np.linspace(0, widthprof, int((widthprof)/5)), width=mxbin)
                    
                pstream = profile(tqdm.tqdm(stream), boxes)

                if len(pstream):
                    logger.info(f"------> Calculated profile for azimuth {azimuth}: {outputfile}; Number of traces in the box: {len(pstream)}\n")
                    pstream.write(outputfile, 'H5')
                else:
                    logger.warning(f"------> No output file written for {outputfile}; Number of traces in the box: {len(pstream)}\n")
            else:
                logger.info(f"----> {outputfile} already exists!")


    logger.info("Plotting the piercing point maps")

    map = plot_merc(resolution='h', llcrnrlon=np.amin(stlons)-1, llcrnrlat=np.amin(stlats)-1,urcrnrlon=np.amax(stlons)+1, urcrnrlat=np.amax(stlats)+1,topo=True)
    
    for lat,lon in zip(stlats,stlons):
        x,y = map(lon, lat)
        map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
    
    xpp,ypp = map(ppoints_df["pplon"].values, ppoints_df["pplat"].values)
    map.plot(xpp, ypp,'x', markersize=5,color='b',markeredgewidth=0.3, zorder=1)
    plt.savefig(destination+'piercing_points_map.'+fig_frmt,dpi=200,bbox_inches='tight')

    plot_bm_azimuth(map,stlon=stlon,stlat=stlat,distval_lat=width_lat,distval_lon=width_lon,ndivlat=ndivlat,ndivlon=ndivlon)
    outimagename = destination+'piercing_points_map_new.'+fig_frmt
    plt.savefig(outimagename,dpi=200,bbox_inches='tight')
    logger.info(f"----> Output image is {outimagename}")
    plt.close('all')


def plot_RF_profile(profilefileloc,destination="./",trimrange=(-5,20)):
    logger = logging.getLogger(__name__)
    logger.info("--> Plotting the RF profile")
    plt.style.use('classic')
    for azimuth in [0,90]:
        inpfiles = glob.glob(profilefileloc+ f"rf_profile_profile{azimuth}_*.h5")

        for inpfile in inpfiles:
            logger.info(f"----> Working on {inpfile}")
            pstream = read_rf(inpfile)
            divparam = inpfile.split("_")[-4:-1]
            divsuffix = inpfile.split("_")[-1].split(".")[0]
            # logger.info(pstream)
            pstream.trim2(trimrange[0], trimrange[1], 'onset')
            for chn in ['L','Q']:
                plt.figure()
                pstream.select(channel='??'+chn).normalize().plot_profile(scale=1.5, top='hist', fillcolors=('r', 'b'))
                plt.gcf().set_size_inches(15, 10)
                plt.title(f'Channel: {chn} Azimuth {azimuth}')
                outputimage = destination+f"{chn}_{azimuth}_{divparam[0]}_{divparam[1]}_{divparam[2]}_{divsuffix}_profile_plot.png"
                plt.savefig(outputimage,dpi=200,bbox_inches='tight')
                logger.info(f"------> Output image is {outputimage}")



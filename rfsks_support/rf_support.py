import matplotlib.pyplot as plt
import tqdm
import numpy as np
import os, glob
from rf import RFStream, read_rf, IterMultipleComponents, get_profile_boxes
from rfsks_support.plotting_map import plot_merc, plot_bm_azimuth
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.other_support import avg
from rfsks_support.profile import profile
import logging



### Compute RF
def compute_rf(dataRFfileloc):
    logger = logging.getLogger(__name__)
    all_rfdatafile = glob.glob(dataRFfileloc+'*-rf_profile_data.h5')
    for rfdatafile in all_rfdatafile:
        network = rfdatafile.split("-")[0]
        station = rfdatafile.split("-")[1]
        rffile = f'{network}-{station}-rf_profile_rfs.h5'
        if not os.path.exists(rffile):
            logger.info(f"--> Computing RF for {rfdatafile}")
            data = read_rf(rfdatafile, 'H5')
            stream = RFStream()
            for stream3c in tqdm.tqdm(IterMultipleComponents(data, 'onset', 3)):
                stream3c.filter('bandpass', freqmin=0.5, freqmax=2)
                if len(stream3c) != 3:
                    continue
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
        # logger.info(stream[i].stats.network,stream[i].stats.station, stream[i].stats.channel)

        kw = {'trim': (-5, 20), 'fillcolors': ('black', 'gray'), 'trace_height': 0.1}

        num_trace=len(stream.select(component='L', station=stream[i].stats.station).sort(['back_azimuth']))
        if num_trace > 0:
            try:
                stream.select(component='L', station=stream[i].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                plt.savefig(destImg + f"{stream[i].stats.station}"+'_L.'+fig_frmt)
                stream.select(component='Q', station=stream[i].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                plt.savefig(destImg + f"{stream[i].stats.station}"+'_Q.'+fig_frmt)
                logger.info("----> Working on {}/{}, network: {} station name: {} total traces: {}".format(i+1,len(rffiles),stream[i].stats.network, stream[i].stats.station, num_trace))
            except Exception as e:
                logger.error("Unexpected error", exc_info=True)
        else:
            logger.info("----> {} traces for network: {} station: {}".format(num_trace,stream[i].stats.network, stream[i].stats.station))


def plot_pp_profile_map(dataRFfileloc,profilefileloc,catalogtxtloc,topo=True,destination="./",depth=70,fig_frmt="png",ndiv = 2):
    logger = logging.getLogger(__name__)
    df = pd.read_csv(catalogtxtloc+'all_stations_rf_retrieved.txt',sep="|")

    logger.info("--> Plotting the piercing points on map")
    clon = df['Longitude'].mean()
    clat = df['Latitude'].mean()
    mnlon,mxlon,mnlat,mxlat = df['Longitude'].min(),df['Longitude'].max(),df['Latitude'].min(),df['Latitude'].max()
        
    rffiles = glob.glob(dataRFfileloc+'*-rf_profile_rfs.h5')
    stream = read_rf(rffiles[0], 'H5')
    map = plot_merc(resolution='h', llcrnrlon=mnlon-1, llcrnrlat=mnlat-1,urcrnrlon=mxlon+1, urcrnrlat=mxlat+1,topo=True)
    x,y = map(stream[0].stats.station_longitude, stream[0].stats.station_latitude)
    map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
    for rffile in rffiles[1:]:
        logger.info(f"----> for {rffile}")
        st_tmp = read_rf(rffile, 'H5')
        stream += read_rf(rffile, 'H5')
        x,y = map(st_tmp[0].stats.station_longitude, st_tmp[0].stats.station_latitude)
        map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
    ppoints = stream.ppoints(depth)

    xpp,ypp = map(ppoints[:,1], ppoints[:,0])
    map.plot(xpp, ypp,'x', markersize=5,color='b',markeredgewidth=0.3, zorder=1)
    plt.savefig(destination+'piercing_points_map.'+fig_frmt,dpi=200,bbox_inches='tight')
    plt.close('all')

    
    width_lat = (np.amax(ppoints[:,0])-np.amin(ppoints[:,0])) * 111.2
    # width_lat +=0.3*width_lat
    width_lon = (np.amax(ppoints[:,1])-np.amin(ppoints[:,1])) * 111.2
    # width_lon +=0.3*width_lon
        

    for azimuth in [0,90]:
        if azimuth == 0:
            mxbin = width_lat
        elif azimuth == 90:
            mxbin = width_lon

        divlocs = 0

        
        for dd in range(ndiv):
            initdiv = divlocs
            enddiv = divlocs+(mxbin/ndiv)   

            outputfile = profilefileloc+ f"rf_profile_profile{azimuth}_{int(initdiv)}-{int(enddiv)}.h5"
            if not os.path.exists(outputfile):
                boxes = get_profile_boxes((clat, clon), azimuth, np.linspace(initdiv, enddiv, int(np.abs(enddiv-initdiv)*5)), width=mxbin)
                pstream = profile(tqdm.tqdm(stream), boxes)
                if len(pstream):
                    logger.info(f"------> Calculated profile for azimuth {azimuth}: {outputfile}\n")
                    pstream.write(outputfile, 'H5')
                else:
                    logger.warning(f"------> No output file written for {outputfile}; Number of traces in the box: {len(pstream)}\n")
            else:
                logger.info(f"----> {outputfile} already exists!")
            divlocs = enddiv

    ## Profile 2
    rffiles = glob.glob(dataRFfileloc+'*-rf_profile_rfs.h5')
    stream = read_rf(rffiles[0], 'H5')
    map = plot_merc(resolution='h', llcrnrlon=mnlon-1, llcrnrlat=mnlat-1,urcrnrlon=mxlon+1, urcrnrlat=mxlat+1,topo=topo)
    x,y = map(stream[0].stats.station_longitude, stream[0].stats.station_latitude)
    map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
    for rffile in rffiles[1:]:
        logger.info(f"----> for {rffile}")
        st_tmp = read_rf(rffile, 'H5')
        stream += read_rf(rffile, 'H5')
        x,y = map(st_tmp[0].stats.station_longitude, st_tmp[0].stats.station_latitude)
        map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
    
    ppoints = stream.ppoints(depth)

    xpp,ypp = map(ppoints[:,1], ppoints[:,0])
    map.plot(xpp, ypp,'x', markersize=5,color='b',markeredgewidth=0.3, zorder=1)
    plot_bm_azimuth(map,clon=clon,clat=clat,distval_lat=width_lat,distval_lon=width_lon,ndiv=ndiv)

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
        # logger.info(inpfiles)
        # inpfile = profilefileloc+ 'rf_profile_profile'+str(azimuth)+'.h5'
        for inpfile in inpfiles:
            logger.info(f"----> Working on {inpfile}")
            pstream = read_rf(inpfile)
            divparam = inpfile.split("_")[-1].split(".")[0]
            # logger.info(pstream)
            pstream.trim2(trimrange[0], trimrange[1], 'onset')
            for chn in ['L','Q']:
                plt.figure()
                pstream.select(channel='??'+chn).normalize().plot_profile(scale=1.5, top='hist', fillcolors=('r', 'b'))
                plt.gcf().set_size_inches(15, 10)
                plt.title(f'Channel: {chn} Azimuth {azimuth}')
                outputimage = destination+f"{chn}_{azimuth}_{divparam}_profile_plot.png"
                plt.savefig(outputimage,dpi=200,bbox_inches='tight')
                logger.info(f"------> Output image is {outputimage}")



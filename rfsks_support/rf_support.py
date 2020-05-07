import matplotlib.pyplot as plt
import tqdm, sys
import numpy as np
import os, glob
from rf import RFStream, read_rf, IterMultipleComponents , get_profile_boxes
from rfsks_support.plotting_map import plot_merc, plot_bm_azimuth
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from rfsks_support.other_support import avg
from rfsks_support.profile import profile
# from rfsks_support.rfsks_extras import get_profile_boxes
import logging, yaml



# advinputRF = "Settings/advRFparam.txt"
# inpRF = pd.read_csv(advinputRF,sep="|",index_col ='PARAMETERS')
with open('Settings/advRFparam.yaml') as f:
    inpRFdict = yaml.load(f, Loader=yaml.FullLoader)

### Compute RF
def compute_rf(dataRFfileloc):
    logger = logging.getLogger(__name__)
    all_rfdatafile = glob.glob(dataRFfileloc+f"*-{str(inpRFdict['filenames']['data_rf_suffix'])}.h5")
    for jj,rfdatafile in enumerate(all_rfdatafile):
        network = rfdatafile.split("-")[0]
        station = rfdatafile.split("-")[1]
        rffile = f"{network}-{station}-{str(inpRFdict['filenames']['rf_compute_data_suffix'])}.h5"
        datatmp = read_rf(rfdatafile, 'H5')
        if not os.path.exists(rffile):
            logger.info(f"--> Computing RF for {rfdatafile}, {jj+1}/{len(all_rfdatafile)}")
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
                        # print('Wrong trace length ', lentr,lengt)
                        continue
                
                stream3c.filter('bandpass', freqmin=float(inpRFdict['rf_filter_settings']['minfreq']), freqmax=float(inpRFdict['rf_filter_settings']['maxfreq']))

                try:
                    stream3c.rf()
                except Exception as e:
                    logger.error("Problem applying rf method", exc_info=True)
                stream3c.moveout()
                stream.extend(stream3c)
            stream.write(rffile, 'H5')
        else:
            # logger.info(f"--> {rffile} already exists!, {jj}/{len(all_rfdatafile)}")
            logger.info(f"--> Verifying RF computation {jj+1}/{len(all_rfdatafile)}")

def plot_RF(dataRFfileloc,destImg,fig_frmt="png"):
    logger = logging.getLogger(__name__)
    logger.info("--> Plotting the receiver functions")
    rffiles = glob.glob(dataRFfileloc+f"*-{str(inpRFdict['filenames']['rf_compute_data_suffix'])}.h5")
    for i,rffile in enumerate(rffiles):
        stream = read_rf(rffile, 'H5')
        outfigname1 = destImg + f"{stream[0].stats.station}"+'_L.'+fig_frmt
        outfigname2 = destImg + f"{stream[0].stats.station}"+'_Q.'+fig_frmt

        if not os.path.exists(outfigname1) and not os.path.exists(outfigname2):
            kw = {'trim': (int(inpRFdict['rf_display_settings']['trim_min']), int(inpRFdict['rf_display_settings']['trim_max'])), 'fillcolors': ('black', 'gray'), 'trace_height': float(inpRFdict['rf_display_settings']['trace_height'])}
            if str(inpRFdict['rf_display_settings']['rf_info']) == "default":
                kw['info'] = (('back_azimuth', u'baz (°)', 'C0'),('distance', u'dist (°)', 'C3'))
            else:
                kw['info'] = None

            num_trace=len(stream.select(component='L', station=stream[0].stats.station).sort(['back_azimuth']))
            if num_trace > 0:
                try:
                    stream.select(component='L', station=stream[0].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                    plt.savefig(destImg + f"{stream[0].stats.station}"+'_L.'+fig_frmt)
                    plt.close('all')
                    stream.select(component='Q', station=stream[0].stats.station).sort(['back_azimuth']).plot_rf(**kw)
                    plt.savefig(destImg + f"{stream[0].stats.station}"+'_Q.'+fig_frmt)
                    plt.close('all')
                    logger.info("----> Plotting RF {}/{}, {}-{} Traces: {}".format(i+1,len(rffiles),stream[0].stats.network, stream[0].stats.station, num_trace))
                except Exception as e:
                    logger.error("Unexpected error", exc_info=True)
            else:
                logger.info("----> {} traces for {}-{}".format(num_trace,stream[0].stats.network, stream[0].stats.station))

def write_profile_boxes(outputfile,stream,azimuth,stlat,stlon,initdiv,enddiv,widthprof,mxbin,dbff,done_box_list):
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(outputfile) and outputfile not in done_box_list:
        logger.info("Calculating the boxes")
        if azimuth == 90:
            logger.info('For az: {} startlat: {:.4f}, startlon: {:.4f}, endlon: {:.4f}, length: {}'.format(azimuth,stlat,initdiv,enddiv,widthprof))
            boxes = get_profile_boxes((stlat, initdiv), azimuth, np.linspace(0, widthprof, int((widthprof)/5)), width=mxbin)
        elif azimuth == 0:
            # logger.info(f'For az: {azimuth} startlat: {initdiv}, startlon: {stlon}, endlat: {enddiv}, length: {widthprof}')
            logger.info('For az: {} startlat: {:.4f}, startlon: {:.4f}, endlat: {:.4f}, length: {}'.format(azimuth,initdiv,stlon,enddiv,widthprof))
            boxes = get_profile_boxes((initdiv, stlon), azimuth, np.linspace(0, widthprof, int((widthprof)/5)), width=mxbin)
            
        pstream = profile(tqdm.tqdm(stream), boxes)
        if len(pstream):
            logger.info("------> Calculated profile for azimuth {}: {}; Number of traces in the box: {}\n".format(azimuth,outputfile,len(pstream)))
            pstream.write(outputfile, 'H5')
            dbff.write("{}\n".format(outputfile))
        else:
            logger.warning(f"------> No output file written for {outputfile}; Number of traces in the box: {len(pstream)}\n")
            dbff.write("{}\n".format(outputfile))
    else:
        logger.info(f"----> {outputfile} already exists!")


def plot_pp_profile_map(dataRFfileloc,profilefileloc,catalogtxtloc,topo=True,destination="./",depth=int(inpRFdict['rf_profile_settings']['ppdepth']),fig_frmt="png",ndivlat = 2, ndivlon=3):
    logger = logging.getLogger(__name__)
    list_of_dfs = []
    list_of_streams = []

    logger.info("--> Plotting the piercing points on map")

    stlons,stlats=[],[]
        
    rffiles = glob.glob(dataRFfileloc+f"*-{str(inpRFdict['filenames']['rf_compute_data_suffix'])}.h5")
    if len(rffiles):
        stream = read_rf(rffiles[0], 'H5')

        if not os.path.exists(profilefileloc+"ppoints_df.pkl"):
            ppoints_tmp = stream.ppoints(depth)
            ppdf_tmp = pd.DataFrame({"pplon":ppoints_tmp[:,1],"pplat":ppoints_tmp[:,0]})
            list_of_dfs.append(ppdf_tmp)

            stlons.append(stream[0].stats.station_longitude)
            stlats.append(stream[0].stats.station_latitude)
            for ii,rffile in enumerate(rffiles[1:]):
                logger.info(f"----> reading profile {rffile}, {ii+1}/{len(rffiles[1:])}")
                try:
                    st_tmp = read_rf(rffile, 'H5')
                    list_of_streams.append(st_tmp)
                    ppoints_tmp = st_tmp.ppoints(depth)
                    ppdf_tmp = pd.DataFrame({"pplon":ppoints_tmp[:,1],"pplat":ppoints_tmp[:,0]})
                    list_of_dfs.append(ppdf_tmp)
                    stlons.append(st_tmp[0].stats.station_longitude)
                    stlats.append(st_tmp[0].stats.station_latitude)
                except:
                    logger.error("Error", exc_info=True)

            logger.info("Calculating profile")
            ppoints_df = pd.DataFrame({"list_of_dfs":list_of_dfs,"stlons": stlons, "stlats": stlats})
            ppoints_lst_df = pd.DataFrame({ "list_of_streams":list_of_streams})
            ppoints_df.to_pickle(profilefileloc+"ppoints_df.pkl")
            ppoints_lst_df.to_pickle(profilefileloc+"ppoints_lst_df.pkl")

            list_of_dfs = ppoints_df["list_of_dfs"].values
            pp_lon_lat = pd.concat(list_of_dfs)
            list_of_streams = ppoints_lst_df["list_of_streams"].values
            stlons = ppoints_df["stlons"].values
            stlats = ppoints_df["stlats"].values
        else:
            ppoints_df = pd.read_pickle(profilefileloc+"ppoints_df.pkl")
            ppoints_lst_df = pd.read_pickle(profilefileloc+"ppoints_lst_df.pkl")

            list_of_dfs = ppoints_df["list_of_dfs"].values
            pp_lon_lat = pd.concat(list_of_dfs)
            list_of_streams = ppoints_lst_df["list_of_streams"].values
            stlons = ppoints_df["stlons"].values
            stlats = ppoints_df["stlats"].values



        for st in list_of_streams:
            for tr in st:
                stream.append(tr)


        abs_latvals = np.absolute(pp_lon_lat["pplat"].values)
        abs_lonvals = np.absolute(pp_lon_lat["pplon"].values)
        degkmfac = 111.2
        width_lat = (np.amax(abs_latvals)-np.amin(abs_latvals)) * degkmfac
        width_lon = (np.amax(abs_lonvals)-np.amin(abs_lonvals)) * degkmfac
            
        ndivlat += 1
        ndivlon += 1
        stlat = pp_lon_lat["pplat"].min()
        stlon = pp_lon_lat["pplon"].min()
        for azimuth in [0,90]:
            if azimuth == 90:
                mxbin = width_lat
                ndiv=ndivlat
                divisions = np.linspace(stlon,stlon+mxbin/degkmfac,ndivlon)
            elif azimuth == 0:
                mxbin = width_lon
                ndiv=ndivlon
                divisions = np.linspace(stlat,stlat+mxbin/degkmfac,ndivlat)
            
            ## storing done boxing info
            done_boxing = profilefileloc+'done_boxing.txt'
            if not os.path.exists(done_boxing):
                dbff = open(done_boxing,'w')
                dbff.write("outputfile\n")
                done_box_list = []
            else:
                dbff = open(done_boxing,'a')
                done_box_array_df = pd.read_csv(done_boxing)
                done_box_array_slice = done_box_array_df.iloc[:,0]
                done_box_list = done_box_array_slice.values.tolist()

            for n in range(len(divisions)-1):
                initdiv = divisions[n]
                enddiv = divisions[n+1] #length of profile
                widthprof = int(np.abs(enddiv-initdiv)*degkmfac) #width of profile
                
                outputfile = profilefileloc+ f"{str(inpRFdict['filenames']['rfprofile_compute_result_prefix'])}{azimuth}_{int(initdiv)}_{int(enddiv)}_{widthprof}_{n}.h5"
                write_profile_boxes(outputfile,stream,azimuth,stlat,stlon,initdiv,enddiv,widthprof,mxbin,dbff,done_box_list)



        logger.info("Plotting the piercing points map")
        outimagename = destination+'piercing_points_map_new.'+fig_frmt
        if not os.path.exists(outimagename):
            map = plot_merc(resolution='h', llcrnrlon=np.amin(stlons)-1, llcrnrlat=np.amin(stlats)-1,urcrnrlon=np.amax(stlons)+1, urcrnrlat=np.amax(stlats)+1,topo=True)
            
            for lat,lon in zip(stlats,stlons):
                x,y = map(lon, lat)
                map.plot(x, y,'^', markersize=10,color='r',markeredgecolor='k',markeredgewidth=0.3, zorder=2)
            
            xpp,ypp = map(pp_lon_lat["pplon"].values, pp_lon_lat["pplat"].values)
            map.plot(xpp, ypp,'x', markersize=5,color='b',markeredgewidth=0.3, zorder=1)
            plt.savefig(destination+'piercing_points_map.'+fig_frmt,dpi=200,bbox_inches='tight')

            plot_bm_azimuth(map,stlon=stlon,stlat=stlat,distval_lat=width_lat,distval_lon=width_lon,ndivlat=ndivlat,ndivlon=ndivlon)
            plt.savefig(outimagename,dpi=200,bbox_inches='tight')
            logger.info(f"----> PP map: {outimagename}")
            plt.close('all')
    

def plot_RF_profile(profilefileloc,destination="./",trimrange=(int(inpRFdict['rf_display_settings']['trim_min']), int(inpRFdict['rf_display_settings']['trim_max']))):
    logger = logging.getLogger(__name__)
    logger.info("--> Plotting the RF profile")
    # plt.style.use('classic')
    for azimuth in [0,90]:
        inpfiles = glob.glob(profilefileloc+ f"{str(inpRFdict['filenames']['rfprofile_compute_result_prefix'])}{azimuth}_*.h5")
        if len(inpfiles):
            for inpfile in inpfiles:
                logger.info(f"----> RF profile {inpfile}")
                pstream = read_rf(inpfile)
                divparam = inpfile.split("_")[-4:-1]
                divsuffix = inpfile.split("_")[-1].split(".")[0]
                
                pstream.trim2(trimrange[0], trimrange[1], 'onset')
                for chn in ['L','Q']:
                    outputimage = destination+f"{chn}_{azimuth}_{divsuffix}_profile.png"
                    if not os.path.exists(outputimage):
                        plt.figure()
                        pstream.select(channel='??'+chn).normalize().plot_profile(scale=1.5, top='hist', fillcolors=('r', 'b'))
                        plt.gcf().set_size_inches(15, 10)
                        plt.title(f'Channel: {chn} Azimuth {azimuth}')
                        plt.savefig(outputimage,dpi=200,bbox_inches='tight')
                        logger.info(f"------> Output image: {outputimage}")
    plt.close('all')


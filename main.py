### Some doc at: https://github.com/trichter/rf @cplegendre
import os, sys, glob
import obspy
import pandas as pd
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from rfsks_support.other_support import create_dir, rem_dir, read_directories, setup_logging
import rfsks_support.rf_support as rfs
import rfsks_support.sks_support as skss
from rfsks_support.download_large_data import downloadDataclass
#######################
## Logging
import logging, requests
logfiles = glob.glob("*.log")
for log in logfiles:
    if os.path.exists(log):
        os.remove(log)
setup_logging()
logger = logging.getLogger(__name__)

############################


inp = pd.read_csv("input_parameters.txt",sep="|",index_col ='PARAMETERS')
res_dir = 'results/'
dirs = read_directories(res_dir)


## Input parameters  ## General
fresh_start=int(inp.loc['fresh_start','VALUES'])       #0/1
mnlong,mxlong=float(inp.loc['mnlong','VALUES']),float(inp.loc['mxlong','VALUES'])   #min and max longitude 
mnlat,mxlat=float(inp.loc['mnlat','VALUES']),float(inp.loc['mxlat','VALUES'])   #min and max latitude 
client=inp.loc['client','VALUES'].split(",")   #client name to retrieve the data
network=str(inp.loc['network','VALUES'])
station=str(inp.loc['station','VALUES'])
fig_frmt="png"

## Input parameters  ## User's choice
makeRF=int(inp.loc['makeRF','VALUES'])                ######   0/1
makeSKS=int(inp.loc['makeSKS','VALUES'])               ######   0/1

## Input parameters  ## Plotting
plot_stations=int(inp.loc['plot_stations','VALUES'])
plot_events=int(inp.loc['plot_events','VALUES'])
plot_RF = int(inp.loc['plot_RF','VALUES']) #Plotting the receiver functions
plot_ppoints=int(inp.loc['plot_ppoints','VALUES'])
plot_RF_profile = int(inp.loc['plot_RF_profile','VALUES'])
plot_SKS_stations=int(inp.loc['plot_SKS_stations','VALUES'])
plot_SKS_events=int(inp.loc['plot_SKS_events','VALUES'])
plot_SKS = int(inp.loc['plot_SKS','VALUES']) #Plotting the receiver functions
picking_SKS=int(inp.loc['picking_SKS','VALUES'])




## Input parameters  ## RF
minradiusRF=float(inp.loc['minradiusRF','VALUES'])
maxradiusRF=float(inp.loc['maxradiusRF','VALUES'])
minmagnitudeRF=float(inp.loc['minmagnitudeRF','VALUES'])
maxmagnitudeRF=float(inp.loc['maxmagnitudeRF','VALUES'])

## Input parameters  ## SKS
minradiusSKS=float(inp.loc['minradiusSKS','VALUES'])
maxradiusSKS=float(inp.loc['maxradiusSKS','VALUES'])
minmagnitudeSKS=float(inp.loc['minmagnitudeSKS','VALUES'])
maxmagnitudeSKS=float(inp.loc['maxmagnitudeSKS','VALUES'])

## Download data
download_data_RF = int(inp.loc['download_data_RF','VALUES'])
download_data_SKS = int(inp.loc['download_data_SKS','VALUES'])

## Log
logger.info(f"Running the program for makeRF: {makeRF}; makeSKS: {makeSKS}")



invRFfile = str(dirs.loc['RFinfoloc','DIR_NAME']) + 'rf_stations.xml'
RFsta = str(dirs.loc['RFinfoloc','DIR_NAME']) + 'all_stations_RF.txt'


## Defining paths for SKS
invSKSfile = str(dirs.loc['SKSinfoloc','DIR_NAME']) + 'sks_stations.xml'
SKSsta = str(dirs.loc['SKSinfoloc','DIR_NAME']) + 'stations_SKS.txt'


#############################################################
#############################################################
#####################                 #######################
##################### Initialisation  #######################
#####################                 #######################
#############################################################
#############################################################
if fresh_start:
    rem_dir(res_dir)

## Creating directories
for direc in list(dirs['DIR_NAME']):
    if not os.path.exists(direc):
        create_dir(direc)

if inp.loc['locations','VALUES'] is np.nan:
    locations=[""]
else:
    locations = inp.loc['locations','VALUES'].split(",")
if not len(locations):
    locations=[""]
#############################################################
#############################################################
#####################                 #######################
#####################       RF        #######################
#####################                 #######################
#############################################################
#############################################################
if makeRF:
    logger.info("WORKING ON RF")
    logger.info("WORKING ON RF")
    logger.info("# Initializing the downloadDataclass")
    logger.info("# Initializing the downloadDataclass")

    rf_data=downloadDataclass(inventoryfile=invRFfile,inventorytxtfile=RFsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='RF')
    ## Obtain inventory and events info
    if int(inp.loc['obtain_inventory','VALUES']):
        if not os.path.exists(invRFfile):
                try:
                    logger.info("Trying to operate the get_stnxml method")
                    logger.info("\n")
                    logger.info("## Operating get_stnxml method")
                    rf_data.get_stnxml(network=network, station=station)
                except Exception as e:
                    logger.error('Error occurred ' + str(e))
                    logger.info("Timeout while requesting...Please try again after some time")
                    sys.exit()
        
        rf_data.obtain_events(catalogxmlloc=str(dirs.loc['RFinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['RFinfoloc','DIR_NAME']),minmagnitude=minmagnitudeRF,maxmagnitude=maxmagnitudeRF)

    ## Download waveforms
    if download_data_RF:
        
        retrived_stn_file = str(dirs.loc['RFinfoloc','DIR_NAME'])+'all_stations_rf_retrieved.txt'
        if not os.path.exists(retrived_stn_file):
            if not os.path.exists(invRFfile):
                try:
                    logger.info("Trying get_stnxml method for RF")
                    logger.info("\n")
                    logger.info("## Operating get_stnxml method")
                    rf_data.get_stnxml(network=network, station=station)
                except requests.Timeout as err:
                    logger.error({"message": err.message})
                    logger.info("Timeout while requesting...Please try again after some time")
                    sys.exit()

            all_stations_df = pd.read_csv(RFsta, sep="|")
            nets = all_stations_df['#Network'].values
            stns = all_stations_df['Station'].values
            

            for net, sta in zip(nets,stns):
                catfile = str(dirs.loc['RFinfoloc','DIR_NAME'])+f"{net}-{sta}-events-info-RF.txt"
                if not os.path.exists(catfile):
                    rf_data.obtain_events(catalogxmlloc=str(dirs.loc['RFinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['RFinfoloc','DIR_NAME']),minmagnitude=minmagnitudeRF,maxmagnitude=maxmagnitudeRF)
                    break

            total_events=0
            for net, sta in zip(nets,stns):
                catfile = str(dirs.loc['RFinfoloc','DIR_NAME'])+f"{net}-{sta}-events-info-RF.txt"
                total_events += int(pd.read_csv(catfile,sep="|",header=None).shape[0])

            if total_events:
                logger.info("\n")
                logger.info("## Operating download method")
                rf_data.download_data(catalogxmlloc=str(dirs.loc['RFinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['RFinfoloc','DIR_NAME']),datafileloc=str(dirs.loc['dataRFfileloc','DIR_NAME']),tot_evnt_stns=total_events, plot_stations=plot_stations, plot_events=plot_events,dest_map=str(dirs.loc['RFstaevnloc','DIR_NAME']),locations=locations)
            else:
                logger.info("No events found!")

    

    if plot_RF:
        try:
            logger.info("\n## Computing RF")
            rfs.compute_rf(str(dirs.loc['dataRFfileloc','DIR_NAME']))
            logger.info("\n")
            logger.info("## Operating plot_RF method")
            rfs.plot_RF(str(dirs.loc['dataRFfileloc','DIR_NAME']),destImg=str(dirs.loc['plotRFloc','DIR_NAME']))
        except Exception as e:
            logger.info(e)

    if plot_ppoints:
        try:
            logger.info("\n")
            logger.info("## Operating plot_priercingpoints_RF method")
            rfs.plot_pp_profile_map(str(dirs.loc['dataRFfileloc','DIR_NAME']),str(dirs.loc['dataRFfileloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['RFinfoloc','DIR_NAME']),destination=str(dirs.loc['RFprofilemaploc','DIR_NAME']),depth=70,ndiv = int(inp.loc['num_profile_divs','VALUES']))
            # rfs.plot_pp_profile_map_new(str(dirs.loc['dataRFfileloc','DIR_NAME']),str(dirs.loc['dataRFfileloc','DIR_NAME']),catalogtxtloc=catalogtxtloc,destination=str(dirs.loc['RFprofilemaploc','DIR_NAME']),depth=70)

            if plot_RF_profile:
                logger.info("\n")
                logger.info("## Operating plot_RF_profile method")
                rfs.plot_RF_profile(str(dirs.loc['dataRFfileloc','DIR_NAME']),destination=str(dirs.loc['RFprofilemaploc','DIR_NAME']))
        except Exception as e:
            logger.info(e)



#############################################################
#############################################################
#####################                 #######################
#####################      SKS        #######################
#####################                 #######################
#############################################################
#############################################################
if makeSKS:
    # logger.info("\nWORKING ON SKS")
    # logger.info("\n# Initializing the downloadDataclass")
    logger.info("\nWORKING ON SKS")
    logger.info("\n# Initializing the downloadDataclass")
    sks_data=downloadDataclass(inventoryfile=invSKSfile,inventorytxtfile=SKSsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='SKS')

    ## Obtain inventory and events info
    if int(inp.loc['obtain_inventory_SKS','VALUES']):
        if not os.path.exists(invSKSfile):
            try:
                logger.info("Trying to operate the get_stnxml method")
                logger.info("\n")
                logger.info("## Operating get_stnxml method")
                sks_data.get_stnxml(network=network, station=station)
            except requests.Timeout as err:
                logger.error({"message": err.message})
                logger.info("Timeout while requesting...Please try again after some time")
                sys.exit()
        
        sks_data.obtain_events(catalogxmlloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),minmagnitude=minmagnitudeSKS,maxmagnitude=maxmagnitudeSKS)
        
    

    ## Download waveforms
    if download_data_SKS:
        if not os.path.exists(invSKSfile):
            try:
                logger.info("\n")
                logger.info("## Operating get_stnxml method")
                sks_data.get_stnxml(network=network, station=station)
            except:
                logger.info("Timeout while requesting...Please try again after some time")
                sys.exit()

        all_stations_df = pd.read_csv(SKSsta, sep="|")
        nets = all_stations_df['#Network'].values
        stns = all_stations_df['Station'].values
        

        for net, sta in zip(nets,stns):
            catfile = str(dirs.loc['SKSinfoloc','DIR_NAME'])+f"{net}-{sta}-events-info-SKS.txt"
            if not os.path.exists(catfile):
                logger.info(f"{catfile} does not exist!\nObtaining {catfile}")
                sks_data.obtain_events(catalogxmlloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),minmagnitude=minmagnitudeSKS,maxmagnitude=maxmagnitudeSKS)
                break

        total_events=0
        for net, sta in zip(nets,stns):
            catfile = str(dirs.loc['SKSinfoloc','DIR_NAME'])+f"{net}-{sta}-events-info-SKS.txt"
            total_events += int(pd.read_csv(catfile,sep="|",header=None).shape[0])

        if total_events:
            for net, sta in zip(nets,stns):
                if not os.path.exists(str(dirs.loc['dataSKSfileloc','DIR_NAME'])+f"{net}-{sta}-sks_profile_data.h5"):
                    logger.info(f"Downloading {str(dirs.loc['dataSKSfileloc','DIR_NAME'])}+{net}-{sta}-sks_profile_data.h5")
                    logger.info("\n")
                    logger.info("## Operating download method")
                    sks_data.download_data(catalogxmlloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['SKSinfoloc','DIR_NAME']),datafileloc=str(dirs.loc['dataSKSfileloc','DIR_NAME']),tot_evnt_stns=total_events, plot_stations=plot_SKS_stations, plot_events=plot_SKS_events,dest_map=str(dirs.loc['SKSstaevnloc','DIR_NAME']))
                else:
                    logger.info(f"{dirs.loc['dataSKSfileloc','DIR_NAME']} {net}-{sta}-sks_profile_data.h5 already exists!")

        else:
            logger.info("No events found!")


    if picking_SKS:
        logger.info("\n")
        logger.info("## Pre-processing")
        plot_traces_ENZ=int(inp.loc['plot_traces_ENZ','VALUES'])
        plot_traces_RTZ=int(inp.loc['plot_traces_RTZ','VALUES'])
        plot_trigger=int(inp.loc['plot_trigger','VALUES'])
        plot_SKS_measure=int(inp.loc['plot_SKS_measure','VALUES'])

        trace_loc_ENZ = str(dirs.loc['SKStracesloc_ENZ','DIR_NAME']) if plot_traces_ENZ else None
        trace_loc_RTZ = str(dirs.loc['SKStracesloc_RTZ','DIR_NAME']) if plot_traces_RTZ else None
        trigger_loc = str(dirs.loc['SKS_trigger_loc','DIR_NAME']) if plot_trigger else None

        plot_measure_loc = str(dirs.loc['plot_measure_loc','DIR_NAME']) if plot_SKS_measure else None

        skss.SKScalc(str(dirs.loc['dataSKSfileloc','DIR_NAME']),trace_loc_ENZ,trace_loc_RTZ,trigger_loc,plot_measure_loc,method = 'recursive_sta_lta')




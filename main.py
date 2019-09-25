### Some doc at: https://github.com/trichter/rf @cplegendre
import os, sys, glob
import obspy
import pandas as pd
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
import rfsks_support.other_support as oss
# from rfsks_support.other_support import create_dir, oss.rem_dir, oss.read_directories, oss.setup_logging, oss.obtain_inventory_events, oss.select_to_download_events
import rfsks_support.rf_support as rfs
import rfsks_support.sks_support as skss
from rfsks_support.download_large_data import downloadDataclass
import warnings
warnings.filterwarnings("ignore")
import time





inp = pd.read_csv("input_parameters.txt",sep="|",index_col ='PARAMETERS')
res_dir = 'results/'
dirs = oss.read_directories(res_dir)


#######################
## Logging
import logging, requests
logfiles = glob.glob(res_dir+"tmp/*.log")
for log in logfiles:
    if os.path.exists(log):
        os.remove(log)
oss.setup_logging()
logger = logging.getLogger(__name__)
print(f"\nCheck file {res_dir+'tmp/info.log'} for details\n")
time.sleep(5)
############################

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
    oss.rem_dir(res_dir)

## Creating directories
for direc in list(dirs['DIR_NAME']):
    if not os.path.exists(direc):
        oss.create_dir(direc)

## List the given station locations
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
    logger.info("# Initializing the downloadDataclass")

    rf_data=downloadDataclass(inventoryfile=invRFfile,inventorytxtfile=RFsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='RF')
    catalogxmlloc = str(dirs.loc['RFinfoloc','DIR_NAME'])
    ## Obtain inventory and events info
    if int(inp.loc['obtain_inventory','VALUES']):
        logger.info("Obtaining Inventory")
        oss.obtain_inventory_events(rf_data,invRFfile,catalogxmlloc,network,station,dirs,minmagnitudeRF,maxmagnitudeRF)
        logger.info(f"Catalog xml/txt files saved at {dirs.loc['RFinfoloc','DIR_NAME']}")

    ## Download waveforms
    if download_data_RF:
        logger.info("Downloading the RF data")
        try:
            if not os.path.exists(RFsta):
                logger.info(f"{RFsta.split('/')[-1]} does not exist...obtaining")
                oss.obtain_inventory_events(rf_data,invRFfile,catalogxmlloc,network,station,dirs,minmagnitudeRF,maxmagnitudeRF)
            else:
                logger.info(f"{RFsta.split('/')[-1]} exists!")
        except:
            sys.exit()
    

        retrived_stn_file = str(dirs.loc['RFinfoloc','DIR_NAME'])+'all_stations_rf_retrieved.txt'
        if not os.path.exists(retrived_stn_file):
            logger.info(f"{retrived_stn_file} does not exist...obtaining events catalog!")
            catalogloc = str(dirs.loc['RFinfoloc','DIR_NAME'])
            datafileloc=str(dirs.loc['dataRFfileloc','DIR_NAME'])
            dest_map=str(dirs.loc['RFstaevnloc','DIR_NAME'])
            ## The stations list can be edited
            oss.select_to_download_events(catalogloc,datafileloc,dest_map,RFsta,rf_data,minmagnitudeRF,maxmagnitudeRF,plot_stations,plot_events,locations,method='RF')

    

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
   
    logger.info("\nWORKING ON SKS")
    logger.info("\n# Initializing the downloadDataclass")
    sks_data=downloadDataclass(inventoryfile=invSKSfile,inventorytxtfile=SKSsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='SKS')
    catalogxmlloc=str(dirs.loc['SKSinfoloc','DIR_NAME'])



    ## Obtain inventory and events info
    if int(inp.loc['obtain_inventory_SKS','VALUES']):
        
        logger.info("Obtaining Inventory")
        oss.obtain_inventory_events(sks_data,invSKSfile,catalogxmlloc,network,station,dirs,minmagnitudeSKS,maxmagnitudeSKS)
        logger.info(f"Catalog xml/txt files saved at {catalogxmlloc}")
 
        
    

    ## Download waveforms
    if download_data_SKS:
        logger.info("Downloading the SKS data")
        if not os.path.exists(SKSsta):
            logger.info(f"{SKSsta} does not exist...obtaining")
            oss.obtain_inventory_events(sks_data,invSKSfile,catalogxmlloc,network,station,dirs,minmagnitudeSKS,maxmagnitudeSKS)
    

        retrived_stn_file = str(dirs.loc['SKSinfoloc','DIR_NAME'])+'all_stations_rf_retrieved.txt'
        if not os.path.exists(retrived_stn_file):
            logger.info(f"{retrived_stn_file} does not exist...obtaining inventory!")
            catalogloc = str(dirs.loc['SKSinfoloc','DIR_NAME'])
            datafileloc=str(dirs.loc['dataSKSfileloc','DIR_NAME'])
            dest_map=str(dirs.loc['SKSstaevnloc','DIR_NAME'])
            ## The stations list can be edited
            oss.select_to_download_events(catalogloc,datafileloc,dest_map,SKSsta,sks_data,minmagnitudeSKS,maxmagnitudeSKS,plot_SKS_stations,plot_SKS_events,locations,method='SKS')


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




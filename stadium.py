### Some doc at: https://github.com/trichter/rf @cplegendre
import os, sys, glob, yaml
import obspy
import pandas as pd
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
import rfsks_support.other_support as oss
import rfsks_support.rf_support as rfs
import rfsks_support.sks_support as skss
import rfsks_support.summary_support as sumsup
from rfsks_support.download_large_data import downloadDataclass
from rfsks_support.plot_events_map_all import plot_events_map_all
from rfsks_support.plot_station_map_all import plot_station_map_all
from rfsks_support.plotting_h_k import plot_h_kappa
from rfsks_support.calculate_h_k import calc_h_kappa
import warnings
warnings.filterwarnings("ignore")
import time


def main():
    with open('input_file.yaml') as f:
        inp = yaml.load(f, Loader=yaml.FullLoader)

    res_dir = str(inp['project_name']) #'results/'
    dirs,rfdirs,sksdirs,otherdirs = oss.read_directories(res_dir)

    ## Step wise mode
    with open('Settings/stepwise.yaml') as f:
        inp_step = yaml.load(f, Loader=yaml.FullLoader)


    ## Fine tuning of RF
    with open('Settings/advRFparam.yaml') as f:
        inpRFdict = yaml.load(f, Loader=yaml.FullLoader)

    ## Fine tuning of SKS
    with open('Settings/advSKSparam.yaml') as f:
        inpSKSdict = yaml.load(f, Loader=yaml.FullLoader)



    ## Input parameters  ## General
    fresh_start=int(inp['fresh_start'])       #0/1
    mnlong,mxlong=float(inp['mnlong']),float(inp['mxlong'])   #min and max longitude 
    mnlat,mxlat=float(inp['mnlat']),float(inp['mxlat'])   #min and max latitude 
    client=inp_step['data_settings']['client'].split(",")   #client name to retrieve the data
    network=str(inp_step['data_settings']['network'])

    station=str(inp_step['data_settings']['station'])
    
    channel=str(inp_step['data_settings']['channel'])
    


    fig_frmt="png"

    ## Input parameters  ## User's choice
    makeRF=int(inp['makeRF'])                ######   0/1
    makeSKS=int(inp['makeSKS'])               ######   0/1


    ## Input parameters  ## Plotting
    plot_stations=int(inp_step['plot_settings']['plot_stations'])
    plot_events=int(inp_step['plot_settings']['plot_events'])
    plot_all_retrieved = int(inp_step['plot_settings']['plot_all_retrieved_events_stations']) 


    compute_plot_RF = int(inp_step['rf_stepwise']['compute_plot_RF']) #Plotting the receiver functions
    plot_ppoints=int(inp_step['rf_stepwise']['plot_ppoints'])
    plot_RF_profile = int(inp_step['rf_stepwise']['plot_RF_profile'])

    # # plot_SKS = int(inp_step['sks_stepwise']['plot_SKS']) #Plotting the receiver functions
    picking_SKS=int(inp_step['sks_stepwise']['picking_SKS'])


    ## Input parameters  ## RF
    minmagnitudeRF=float(inpRFdict['rf_event_search_settings']['minmagnitudeRF'])
    maxmagnitudeRF=float(inpRFdict['rf_event_search_settings']['maxmagnitudeRF'])


    ## Input parameters  ## SKS
    minmagnitudeSKS=float(inpSKSdict['sks_event_search_settings']['minmagnitudeSKS'])
    maxmagnitudeSKS=float(inpSKSdict['sks_event_search_settings']['maxmagnitudeSKS'])

    ## Download data
    download_data_RF = int(inp_step['rf_stepwise']['download_data_RF'])
    download_data_SKS = int(inp_step['sks_stepwise']['download_data_SKS'])


    ## Station inventory files
    invRFfile = str(dirs.loc['RFinfoloc','DIR_NAME']) + str(inpRFdict['filenames']['invRFfile'])
    RFsta = str(dirs.loc['RFinfoloc','DIR_NAME']) + str(inpRFdict['filenames']['RFsta'])


    ## Defining paths for SKS
    invSKSfile = str(dirs.loc['SKSinfoloc','DIR_NAME'])+ str(inpSKSdict['filenames']['invSKSfile'])
    SKSsta = str(dirs.loc['SKSinfoloc','DIR_NAME']) + str(inpSKSdict['filenames']['SKSsta'])



    #############################################################
    #############################################################
    #####################                 #######################
    ##################### Initialisation  #######################
    #####################                 #######################
    #############################################################
    #############################################################
    if fresh_start:
        if os.path.exists(res_dir):
            response = input("Are you sure you want to start fresh? (Input 'yes' to continue): ")
            if response == "yes":
                oss.rem_dir(res_dir)
            else:
                print(f"Response: {response}, Exiting!")
                sys.exit()

    ## Creating directories
    if makeRF:
        for direc in rfdirs:
            if not os.path.exists(direc):
                oss.create_dir(direc)
    if makeSKS:
        for direc in sksdirs:
            if not os.path.exists(direc):
                oss.create_dir(direc)

    for direc in otherdirs:
        if not os.path.exists(direc):
            oss.create_dir(direc)

    ## List the given station locations
    if inp_step['data_settings']['locations'] is np.nan:
        locations=[""]
    else:
        locations = list(inp_step['data_settings']['locations'])
    if not len(locations):
        locations=[""]

    ##################
    ## Initializing Summary file
    ## Summary File
    sum_file = res_dir+str(inp['summary_file'])
    sum_sup_class = sumsup.sum_support(sum_file,res_dir)
    sum_sup_class.write_initial_summary(mnlong,mxlong,mnlat,mxlat,client,network,station,channel)
    sum_sup_class.makeSKSRF(makeRF,makeSKS)


    

    

    #######################
    ## Logging
    import logging
    logfiles = glob.glob(res_dir+"tmp/*.log")
    for log in logfiles:
        if os.path.exists(log):
            os.remove(log)
    oss.setup_logging(default_level=logging.INFO,dirname=res_dir)
    logger = logging.getLogger(__name__)

    print(f"\nCheck file {res_dir+'tmp/info.log'} for details\n")
    time.sleep(3)
    ## Log
    logger.info(f"Running the program for makeRF: {makeRF}; makeSKS: {makeSKS}")

    ############################
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
        sum_sup_class.write_strings("--> RECEIVER FUNCTIONS PART:")

        rf_data=downloadDataclass(inventoryfile=invRFfile,inventorytxtfile=RFsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='RF',channel=channel)
        catalogxmlloc = str(dirs.loc['RFinfoloc','DIR_NAME'])
        ## Obtain inventory and events info
        if int(inp_step['rf_stepwise']['obtain_inventory_RF']):
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
        

            retrived_stn_file = str(dirs.loc['RFinfoloc','DIR_NAME'])+str(inpRFdict['filenames']['retr_stations'])
            if not os.path.exists(retrived_stn_file):
                logger.info(f"{retrived_stn_file} does not exist...obtaining events catalog..")
                catalogloc = str(dirs.loc['RFinfoloc','DIR_NAME'])
                datafileloc=str(dirs.loc['RFdatafileloc','DIR_NAME'])
                dest_map=str(dirs.loc['RFstaevnloc','DIR_NAME'])
                ## The stations list can be edited
                oss.select_to_download_events(catalogloc,datafileloc,dest_map,RFsta,rf_data,minmagnitudeRF,maxmagnitudeRF,plot_stations,plot_events,locations,method='RF')

        
        if plot_all_retrieved and os.path.exists(str(dirs.loc['RFinfoloc','DIR_NAME'])+str(inpRFdict['filenames']['retr_stations'])):
            logger.info("\n")
            logger.info("## Plotting retrieved stations")
            RFsta_path = "/".join(RFsta.split("/")[0:-1])
            RFsta_prex = RFsta.split("/")[-1].split(".")[0]
            full_RFsta_path = f"{RFsta_path}/{RFsta_prex}_combined.txt"
            plot_events_map_all(all_stations_file = str(dirs.loc['RFinfoloc','DIR_NAME'])+str(inpRFdict['filenames']['retr_stations']))
            plot_station_map_all(retr_stationsfile = str(dirs.loc['RFinfoloc','DIR_NAME'])+str(inpRFdict['filenames']['retr_stations']),all_stationsfile=full_RFsta_path)

        if compute_plot_RF:
            dataRFfileloc = str(dirs.loc['RFdatafileloc','DIR_NAME'])
            all_rfdatafile = glob.glob(dataRFfileloc+f"*-{str(inpRFdict['filenames']['data_rf_suffix'])}.h5")
            if len(all_rfdatafile)>=1:
                try:
                    logger.info("\n## Computing RF")
                    rfs.compute_rf(dataRFfileloc)
                    logger.info("\n")
                    logger.info("## Operating plot_RF method")
                    rfs.plot_RF(dataRFfileloc,destImg=str(dirs.loc['RFplotloc','DIR_NAME']))
                except Exception as e:
                    logger.info(e)
            else:
                logger.error("No RF data files present...download the data")
                sys.exit()

        if plot_ppoints:
            # try:
            logger.info("\n")
            logger.info("## Operating plot_priercingpoints_RF method")
            rfs.plot_pp_profile_map(str(dirs.loc['RFdatafileloc','DIR_NAME']),str(dirs.loc['RFdatafileloc','DIR_NAME']),catalogtxtloc=str(dirs.loc['RFinfoloc','DIR_NAME']),destination=str(dirs.loc['RFprofilemaploc','DIR_NAME']), ndivlat = int(inpRFdict['rf_profile_settings']['num_profile_divs_lat']), ndivlon=int(inpRFdict['rf_profile_settings']['num_profile_divs_lon']))

            if plot_RF_profile:
                logger.info("\n")
                logger.info("## Operating plot_RF_profile method")
                rfs.plot_RF_profile(str(dirs.loc['RFdatafileloc','DIR_NAME']),destination=str(dirs.loc['RFprofilemaploc','DIR_NAME']))
            # except Exception as e:
            #     logger.info(e)

        ## H-kappa calculation
        if int(inpRFdict['filenames']['h_kappa_settings']['plot_h']) or int(inpRFdict['filenames']['h_kappa_settings']['plot_kappa']):
            logger.info("\n")
            logger.info("## H-kappa implementation")
            outloc=str(dirs.loc['RFinfoloc','DIR_NAME'])
            outfile = str(inpRFdict['filenames']['h_kappa_settings']['h_kappa_res_file'])
            calc_h_kappa(outfile = outfile,data_dir_loc = str(dirs.loc['RFdatafileloc','DIR_NAME']), outloc=outloc)
            if os.path.exists(outloc+outfile):
                plot_h_kappa(h_k_file = outloc+outfile,all_stationsfile = str(inpRFdict['filenames']['retr_stations']),plot_h = int(inpRFdict['filenames']['h_kappa_settings']['plot_h']),plot_kappa = int(inpRFdict['filenames']['h_kappa_settings']['plot_kappa']))



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
        sum_sup_class.write_strings("--> SHEAR-WAVE SPLITTING PART:")

        sks_data=downloadDataclass(inventoryfile=invSKSfile,inventorytxtfile=SKSsta,client=client,minlongitude=mnlong,maxlongitude=mxlong,minlatitude=mnlat,maxlatitude=mxlat,fig_frmt=fig_frmt,method='SKS',channel=channel)

        catalogxmlloc=str(dirs.loc['SKSinfoloc','DIR_NAME'])


        ## Obtain inventory and events info
        if int(inp_step['sks_stepwise']['obtain_inventory_SKS']):
            
            logger.info("Obtaining Inventory")
            oss.obtain_inventory_events(sks_data,invSKSfile,catalogxmlloc,network,station,dirs,minmagnitudeSKS,maxmagnitudeSKS)
            logger.info(f"Catalog xml/txt files saved at {catalogxmlloc}")
            sum_sup_class.write_data_summary(SKSsta)

    

        ## Download waveforms
        datafileloc=str(dirs.loc['SKSdatafileloc','DIR_NAME'])
        if download_data_SKS:
            logger.info("Downloading the SKS data")
            if not os.path.exists(SKSsta):
                logger.info(f"{SKSsta} does not exist...obtaining")
                oss.obtain_inventory_events(sks_data,invSKSfile,catalogxmlloc,network,station,dirs,minmagnitudeSKS,maxmagnitudeSKS)
                sum_sup_class.write_data_summary(SKSsta)

            retrived_stn_file = str(dirs.loc['SKSinfoloc','DIR_NAME'])+str(inpSKSdict['filenames']['retr_stations'])
            if not os.path.exists(retrived_stn_file):
                logger.info(f"{retrived_stn_file} does not exist...obtaining inventory!")
                catalogloc = str(dirs.loc['SKSinfoloc','DIR_NAME'])
                dest_map=str(dirs.loc['SKSstaevnloc','DIR_NAME'])
                ## The stations list can be edited
                oss.select_to_download_events(catalogloc,datafileloc,dest_map,SKSsta,sks_data,minmagnitudeSKS,maxmagnitudeSKS,plot_stations,plot_events,locations,method='SKS')
            sum_sup_class.write_data_download_summary(datafileloc,retrived_stn_file)

        if plot_all_retrieved and os.path.exists(str(dirs.loc['SKSinfoloc','DIR_NAME'])+str(inpSKSdict['filenames']['retr_stations'])):
            logger.info("\n")
            logger.info("## Plotting retrieved stations")
            SKSsta_path = "/".join(SKSsta.split("/")[0:-1])
            SKSsta_prex = SKSsta.split("/")[-1].split(".")[0]
            full_SKSsta_path = f"{SKSsta_path}/{SKSsta_prex}_combined.txt"
            # all_station_file = str(dirs.loc['SKSinfoloc','DIR_NAME'])+str(inpSKSdict['filenames']['SKSsta'])
            plot_events_map_all(all_stations_file = str(dirs.loc['SKSinfoloc','DIR_NAME'])+str(inpSKSdict['filenames']['retr_stations']))
            plot_station_map_all(retr_stationsfile = str(dirs.loc['SKSinfoloc','DIR_NAME'])+str(inpSKSdict['filenames']['retr_stations']),all_stationsfile=full_SKSsta_path)
        if len(glob.glob(datafileloc+"*.h5"))>0:
            if picking_SKS:
                logger.info("\n")
                logger.info("## Pre-processing")
                plot_traces_ENZ=int(inp_step['sks_stepwise']['plot_traces_ENZ'])
                plot_traces_RTZ=int(inp_step['sks_stepwise']['plot_traces_RTZ'])
                plot_trigger=int(inp_step['sks_stepwise']['plot_trigger'])
                plot_SKS_measure=int(inp_step['sks_stepwise']['plot_SKS_measure'])

                trace_loc_ENZ = str(dirs.loc['SKStracesloc_ENZ','DIR_NAME']) if plot_traces_ENZ else None
                trace_loc_RTZ = str(dirs.loc['SKStracesloc_RTZ','DIR_NAME']) if plot_traces_RTZ else None
                trigger_loc = str(dirs.loc['SKS_trigger_loc','DIR_NAME']) if plot_trigger else None

                plot_measure_loc = str(dirs.loc['SKSplot_measure_loc','DIR_NAME']) if plot_SKS_measure else None
                sksMeasure = skss.sks_measurements(plot_measure_loc=plot_measure_loc)
                sksMeasure.SKScalc(str(dirs.loc['SKSdatafileloc','DIR_NAME']),trace_loc_ENZ,trace_loc_RTZ,trigger_loc,method = str(inpSKSdict['sks_picking']['picking_algo']['sks_picking_algo']))
                
            # sksMeasure.plot_sks_map()
            sks_measurement_file = plot_measure_loc+"../"+"sks_measurements_all.txt"
            if os.path.exists(SKSsta) and os.path.exists(sks_measurement_file):
                if bool(inp_step['sks_stepwise']['plot_data_nodata_map']):
                    sksMeasure.plot_data_nodata_map(sks_stations_infofile=SKSsta)
        
    sum_sup_class.close_sumfile()
if __name__ == '__main__':
    main()


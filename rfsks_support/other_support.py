import os, tqdm
import shutil
import numpy as np
from obspy import UTCDateTime as UTC
import pandas as pd
import logging
import logging.config
logger = logging.getLogger(__name__)

import yaml

def setup_logging(
    default_path='rfsks_support/logging.yaml',
    default_level=logging.INFO,
    env_key='LOG_CFG'
):
    """Setup logging configuration

    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


def create_dir(direc):
    '''
    Create a directory
    '''
    try:
        os.makedirs(direc, exist_ok=True)
    except OSError:
        logger.info("--> Creation of the directory {} failed".format(direc))
    else:
        logger.info("--> Successfully created the directory {}".format(direc))

def rem_dir(direc):
    '''
    Delete a directory
    '''
    if os.path.exists(direc):
        shutil.rmtree(direc)

def read_directories(res_dir):
    dirs = pd.read_csv("rfsks_support/directories_names.txt",sep="|",index_col ='DIR_VAR')
    dirs['DIR_NAME'] = np.array([res_dir+val for val in dirs['DIR_NAME'].values])
    newdirname=[]
    for direc in dirs['DIR_NAME']:
        if direc[-1]!="/":
            newdirname.append(f"{direc}/")
        else:
            newdirname.append(f"{direc}")
    dirs['DIR_NAME']= np.array(newdirname)
    return dirs


avg = lambda num1,num2: (int(num1)+int(num2))/2.0


def date2time(sta_sdate,sta_edate):
    smonth = f'0{sta_sdate.month}' if sta_sdate.month < 10 else f'{sta_sdate.month}'
    emonth = f'0{sta_edate.month}' if sta_edate.month < 10 else f'{sta_edate.month}'
    sday = f'0{sta_sdate.day}' if sta_sdate.day < 10 else f'{sta_sdate.day}'
    eday = f'0{sta_edate.day}' if sta_edate.day < 10 else f'{sta_edate.day}'
    stime = f'{sta_sdate.year}-{smonth}-{sday}'
    etime = f'{sta_edate.year}-{emonth}-{eday}'

    return UTC(stime), UTC(etime)

def write_station_file(inventorytxtfile,rf_staNetNames,outfile):
    df_stations = pd.read_csv(inventorytxtfile,sep="|")
    allnetsrf = [netsta.split("_")[0] for netsta in set(rf_staNetNames)]
    allstnsrf = [netsta.split("_")[1] for netsta in set(rf_staNetNames)]
    df_stations_new = pd.DataFrame()
    for net, stn in zip(allnetsrf,allstnsrf):
        df_stations_extract = df_stations[(df_stations['#Network']==net) & (df_stations['Station']==stn)]
        df_stations_new = df_stations_new.append(df_stations_extract, ignore_index=True)
    df_stations_new.to_csv(outfile, index = None, header=True, sep = "|")
import os, tqdm
import shutil
import numpy as np
from obspy import UTCDateTime as UTC
import pandas as pd




def create_dir(direc):
    '''
    Create a directory
    '''
    try:
        os.makedirs(direc, exist_ok=True)
    except OSError:
        print ("--> Creation of the directory {} failed".format(direc))
    else:
        print ("--> Successfully created the directory {}".format(direc))

def rem_dir(direc):
    '''
    Delete a directory
    '''
    if os.path.exists(direc):
        shutil.rmtree(direc)


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
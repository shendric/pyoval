# -*- coding: utf-8 -*-

from helpers import EmptyObject
import numpy as np
import datetime
import pytz


def read_embird(filename):
    """
    Returns an object with time/longitude/latitude/thickness
    """
    from netCDF4 import Dataset
    data = EmptyObject()
    f = Dataset(filename)
    for key in f.variables.keys():
        setattr(data, key.lower(), np.array(f.variables[key][:]))
    f.close()
    # Create datetime object
    n_records = len(data.thickness)
    data.dt_time = np.ndarray(shape=(n_records), dtype=object)

    secs = np.floor(data.time)
    msecs = 1e6 * (data.time - np.floor(secs))
    for i in np.arange(n_records):
        data.dt_time[i] = datetime.datetime(data.year[i], data.month[i], data.day[i], tzinfo=pytz.utc) + \
                          datetime.timedelta(seconds=int(secs[i]), microseconds=int(msecs[i]))
    return data

def read_asiras(filename):
    """
    Returns an object with time/longitude/latitude/thickness
    """
    from netCDF4 import Dataset
    data = EmptyObject()
    f = Dataset(filename)
    for key in f.variables.keys():
        setattr(data, key.lower(), np.array(f.variables[key][:]))
        print key.lower()
    f.close()
    return data
    
def read_embird_dat(filename):
    """
    Returns an object with time/longitude/latitude/thickness from an EMBird .dat
    file.
    """
    data=EmptyObject()
    names = ('year','month','day','time','fid','latitude','longitude',
             'distance','thickness','height')
    formats = ('i4','i4','i4','f8','i4','f8','f8','f8','f8','f8')
    content=np.loadtxt(filename,dtype={'names' : names,'formats' : formats})
    for key in content.dtype.names:
        setattr(data,key,content[key])
    n_records = len(data.latitude)
    data.dt_time = np.ndarray(shape=(n_records), dtype=object)
    secs = np.floor(data.time)
    msecs = 1e6 * (data.time - np.floor(secs))
    for i in np.arange(n_records):
        data.dt_time[i] = datetime.datetime(data.year[i], data.month[i], data.day[i], tzinfo=pytz.utc) + \
                          datetime.timedelta(seconds=int(secs[i]), microseconds=int(msecs[i]))
    return data

def read_OIB(filename):
    '''
    Returns an object with time/longitude/latitude/fb/etc
    '''
    data=EmptyObject()
    names = ('latitude','longitude','thickness','thickness_unc','mean_fb','freeboard',
            'fb_unc','snow_depth','snow_depth_unc','n_atm','pcnt_ow',
            'pcnt_thin_ice','pcnt_grey_ice','corr_elev','elev','date','elapsed',
            'atmos_corr','geoid_corr','ellip_corr','tidal_corr',
            'ocean_tide_corr_part','load_tide_corr_part','earth_tide_corr_part',
            'ssh','n_ssh','ssh_sd','ssh_diff','ssh_elapsed','ssh_tp_dist',
            'surface_roughness','ATM_file_name','Tx','Rx','KT19_surf',
            'KT19_int','low_en_corr','sa_int_elev','si_int_elev','my_ice_flag',
            'KT19_unc','empty1','empty2','empty3','empty4','empty5','empty6',
            'empty7','empty8','empty9')
            #Note that there is a difference between IDCSI2 and IDCSI2 QL Files:
            #in IDCSI2 mean_fb == freeboard (ATM_fb) and KT19_unc is Empty0 in 
            #IDCSI2_QL
    formats = ('f8','f8','f8','f8','f8','f8','f8','f8','f8','i4','f8','f8','f8',
               'f8','f8','a8','f8','f8','f8','f8','f8','f8','f8','f8','f8','i4',
               'f8','f8','f8','f8','f8','a27','f8','f8','f8','f8','f8','f8',
               'f8','i4','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')         
    content = np.loadtxt(filename,delimiter=',',
                            dtype = {'names': names,'formats' : formats},
                            skiprows=1)
                            
    
    for key in content.dtype.names:
        setattr(data,key,content[key])
    n_records = len(data.latitude)
    data.dt_time = np.ndarray(shape=(n_records), dtype=object)
    secs = np.floor(data.elapsed)
    msecs = np.floor(1e6 * (data.elapsed - np.floor(secs)))
    for i in np.arange(n_records):
        data.dt_time[i] = datetime.datetime.strptime(data.date[i],"%Y%m%d").replace(tzinfo=pytz.utc) + \
                          datetime.timedelta(seconds=int(secs[i]), microseconds=int(msecs[i]))
        if data.longitude[i] >= 180.0:
            data.longitude[i]=data.longitude[i]-360.0
        if data.freeboard[i]==-99999:
            data.freeboard[i]=np.nan
  
    return data
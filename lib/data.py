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
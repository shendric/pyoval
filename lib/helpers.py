# -*- coding: utf-8 -*-

import datetime

class EmptyObject():
    """ Container for arbitrary python object"""
    pass


def timedelta_str(td):
    """
    Solve an issue with negative timedelta strings
    e.g. '-0:50:00' instead of '-1 day, 23:10:00'
    returns unchanged str for positive timedeltas
    """
    if td.days < 0:
        return '-' + str(datetime.timedelta() - td)
    return str(td)
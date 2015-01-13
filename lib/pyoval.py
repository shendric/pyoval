# -*- coding: utf-8 -*-
#
# License: GNU General Public License (GPL) v3
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Author: Stefan Hendricks
#
# Institution: Alfred-Wegener-Institut, Helmholtz Zentrum für
#              Polar und Meeresforschung
#
# Created for ESA CryoVal-SI under contract AO/1-­/7340/12/NL/CT

__version__ = 1.2
__author__ = "Stefan Hendricks, AWI"


import numpy as np
import os
import sys

# Timing modules
import datetime
import pytz

# Projections
import pyproj

# Filtering
from scipy.ndimage.filters import maximum_filter1d

# Plotting routines -> matplotlib (version 1.3.1 and higher)
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.dates import  DateFormatter, date2num, num2date
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import MultipleLocator
from matplotlib.collections import PatchCollection
from matplotlib.spines import Spine
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# user library
from helpers import timedelta_str


class CS2RefOrbit(object):
    """
    Calculates the corner coordinates of each footprint for a given orbit track. The
    corner coordinates are given as upper right, lower right, upper left, lower left, were
    the "upper coordinates" are situated in the direction of the orbit (might south of lower)

    along-track footprint : spacing of footprint center coordinates
    across-track footprint : 1650m (pulse limited flat-surface footprint at CS-2 bandwidth)
    """
    def __init__(self):
        self.g = pyproj.Geod(ellps='WGS84')
        self.lat_limit = [-90.1, 90.1]
        self.lon_limit = [-180.1, 180.1]
        self._pointdist_tresh = 303.0
        self._sarin_dist_treshold = 1.0
        self._dist_max_filter_width = 11
        self._sin_detection = True

    def from_file(self, filename, corner_coords=False):
        """
        filename : str
            Link to file
        source : bool (optional, default=True)
             true: assume that only center coordinates are in file
             false: Precalculated center & corner coordinates
                    (output of pyoval.Reforbit.save_fp_corner_coords)
         """
        self.fullpath = filename
        self.filename = os.path.basename(filename)

        if corner_coords:
            self.orbit = int(self.filename.split('_')[2])
            self.__readCCFile()
        else:
            self.orbit = int(self.filename.split('_')[2])
            self.__readFile()
            self.__calcTime()
            self.__calc_pointdist()
            self.__correct_sarin_center_coords()
            self.__dist_filter()
            self.__calc_fp_corner_coords()


    def limit_region(self, lon_limit=None, lat_limit=None):
        """
        User inferface to Region of Interest setting
        Keywords:
        lon_limit : numeric list [lon min, lon max] in degree
            longitude range between -180 deg and 180 deg
        lat_limit : numeric list [lon min, lon max] in degree
        """
        if isinstance(lon_limit, list):
            self.lon_limit = lon_limit
        if isinstance(lat_limit, list):
            self.lat_limit = lat_limit

    def set_sin_detection(self, value):
        """
        Switch for automatic detection of SARIn mode (default: True)
        value : Bool
        """
        if not isinstance(value, bool): raise AssertionError('value not of type bool')
        self._sin_detection = value


    def __roi_index(self, lon, lat):
        """ Returns elements and number of given lon, lat points that satisfy ROI criteria"""
        inLonROI = np.logical_and(lon >= self.lon_limit[0],
                                  lon <= self.lon_limit[1])
        inLatROI = np.logical_and(lat >= self.lat_limit[0],
                                  lat <= self.lat_limit[1])
        inROI = np.where(np.logical_and(inLonROI, inLatROI))[0]
        return inROI, len(inROI)


    def __readFile(self):
        """ Reads reference orbit file (only keeps data in ROI if defined) """
        # Read File Content
        content = np.loadtxt(self.fullpath, dtype={'names': ('day', 'secs', 'lat', 'lon'),
                                                   'formats': ('i4', 'f8', 'f8', 'f8')})
        # Limit file content to region of interset (if applicable, default: entire globe)
        inROI, n_records = self.__roi_index(content['lon'], content['lat'])
        # Retrieve data in ROI
        self.lon, self.lat, self.day = content['lon'][inROI], content['lat'][inROI], content['day'][inROI]
        self.__secs = content['secs'][inROI]
        self.n_records = n_records


    def __readCCFile(self):
        # Definition of corner coordinate ascii file
        names = ('_yy', '_mm', '_dd', '_hh', '_mi', '_sec',
                 'lon', 'lat', 'lon_ur', 'lat_ur', 'lon_ul', 'lat_ul', 'lon_lr', 'lat_lr', 'lon_ll', 'lat_ll')
        formats = ('i4', 'i4', 'i4', 'i4', 'i4', 'f4',
                   'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')
        # Read data
        content = np.loadtxt(self.fullpath, dtype={'names': names, 'formats': formats})
        # Set data fields (not time)
        for key in content.dtype.names:
            setattr(self, key, content[key])
        # set time
        self.n_records = len(self.lon)
        self.time = np.ndarray(shape=(self.n_records), dtype=object)
        for i in np.arange(self.n_records):
            self.time[i] = datetime.datetime(self._yy[i], self._mm[i], self._dd[i],
                                             self._hh[i], self._mi[i], int(self._sec[i]),
                                             int( 1e6*(self._sec[i] - np.floor(self._sec[i]))),
                                             tzinfo=pytz.utc)

    def __calcTime(self, start_time_index = 3):
        """ Creates datetime object for each CS-2 point """
        # Calculate time stamp:
        #  1.1.2000 (ESA definition) + number of days (file) + seconds (file)
        self.time = np.ndarray(shape=(self.n_records), dtype=object)
        self.__msecs = self.__secs - np.floor(self.__secs)
        self.__msecs *= 1e6
        self.__secs = np.floor(self.__secs)
        for i in np.arange(self.n_records):
            self.time[i] = datetime.datetime(2000, 1, 1, tzinfo=pytz.utc) + \
                           datetime.timedelta(days=int(self.day[i]),
                                              seconds=int(self.__secs[i]),
                                              microseconds=int(self.__msecs[i]))

    def __correct_sarin_center_coords(self):
        """
        Applies a correction to center points, that are identified
        as SIN (SARin) measurment positions
        XXX: This is a hot fix and in future orbit positions should be used
        """
        if self._sin_detection:
            self.__classify_sarin()
            #self.sarin_debug_map()
            sarin = np.where(self.cs2_is_sarin)[0]
            # Check if SIN section in data
            if True not in self.cs2_is_sarin:
                return
            # Get mean great circle track for SIN secion
            gc0, gc1 = self.__sarin_mean_gc()
            meanorbit = GreatCircleTrack(gc0, gc1)
            meanorbit.add_track(self.lon[sarin], self.lat[sarin])
            meanorbit.project()
            #meanorbit.debug_map()
            self.lon[sarin] = meanorbit.gcp[:, 0]
            self.lat[sarin] = meanorbit.gcp[:, 1]
        else:
            # For compability reason set classification
            self.cs2_is_sarin = np.zeros(shape=(self.n_records), dtype=bool)

    def debug_map(self):
        """ Creates a quick & dirty map of SAR/SIN classification and corner coordinates """
        lat_0 = np.median(self.lat)
        lon_0 = np.median(self.lon)
        plt.figure("SAR/SIN Debug Map")
        m = Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0, lat_ts=lat_0,
                    resolution = 'i', width=5e6, height=5e6)
        m.drawcoastlines(color="black", linewidth=0.1)
        m.fillcontinents(color='#AAAAAA')
        m.drawmapboundary()
        xc, yc = m(self.lon, self.lat)
        m.scatter(xc[np.where(~self.cs2_is_sarin)], yc[np.where(~self.cs2_is_sarin)],
                     s=20, color='blue', edgecolors='none', zorder=100, label='SAR / Unclassified')
        m.scatter(xc[np.where(self.cs2_is_sarin)], yc[np.where(self.cs2_is_sarin)],
                     s=20, color='red', edgecolors='none', zorder=100, label='SIN')
        plt.legend()
        plt.show()

    def __classify_sarin(self):
        """ Classify SARin points based on the variation from the median distance """
        orbit_typical_dist = np.median(self.pointdist)
        dist_max_window = maximum_filter1d(self.pointdist, self._dist_max_filter_width)
        self.cs2_is_sarin = dist_max_window > orbit_typical_dist + self._sarin_dist_treshold
        sarin_list = np.where(self.cs2_is_sarin)[0]
        # Correct for window size of maximum filter
        filter_buffer = (self._dist_max_filter_width-1)/2
        if sarin_list[0] != 0:
            i0 = sarin_list[0]
            i1 = sarin_list[0]+filter_buffer
            self.cs2_is_sarin[i0:i1] = False
        if sarin_list[-1] != len(self.cs2_is_sarin)-1:
            i0 = sarin_list[-1]-filter_buffer+1
            i1 = sarin_list[-1]+1
            self.cs2_is_sarin[i0:i1] = False

    def __sarin_mean_gc(self):
        """
        Approximates SARIn track with great circle
        returns start [lon,lat] and end [lon, lat]
        2) SIN at end of data subset
        3) SIN in the middle of data subset
        4) Subset only SIN (e.g. Wingham Box)
        """
        sarin_list = np.where(self.cs2_is_sarin)[0]
        n = self.n_records-1
        m = len(sarin_list)-1
        # Case 1: SIN at beginning of data subset
        # 0                        n
        # === SIN == -------SAR----
        # 2          10
        if self.cs2_is_sarin[0] and ~self.cs2_is_sarin[n]:
            i0 = sarin_list[m]+2
            i1 = sarin_list[m]+1
            i2 = 0
            gc0, gc1 = self.__project_gc_point(i0, i1, i2)
            return gc0, gc1
        # Case 2: SIN at end of data subset
        # 0                        n
        # ------SAR--- - === SIN ===
        #            0 1           2
        if ~self.cs2_is_sarin[0] and self.cs2_is_sarin[n]:
            i0 = sarin_list[0]-2
            i1 = sarin_list[0]-1
            i2 = n
            gc0, gc1 = self.__project_gc_point(i0, i1, i2)
            return gc0, gc1
        # Case 3: SIN in between
        # 0                        n
        # --SAR-- == SIN == --SAR---
        #       0           1
        if ~self.cs2_is_sarin[0] and ~self.cs2_is_sarin[n]:
            gc0 = [self.lon[sarin_list[0]-1], self.lat[sarin_list[0]-1]]
            gc1 = [self.lon[sarin_list[m]+1], self.lat[sarin_list[m]+1]]
            return gc0, gc1
        # Case 4: Only SIN
        # 0                        n
        # =========== SIN ==========
        # 0                        1
        if self.cs2_is_sarin[0] and self.cs2_is_sarin[n]:
            gc0 = [self.lon[0], self.lat[n]]
            gc1 = [self.lon[0], self.lat[n]]
            return gc0, gc1

    def __project_gc_point(self, i0, i1, i2):
        """
        Get mean circle end points based on two points that define
        initial points and heading
        i0 = index of heading start point
        i2 = index of heading end point
        i3 = point that defines the distance from
        """
        faz, __, __ = self.g.inv(self.lon[i0], self.lat[i0], self.lon[i1], self.lat[i1])
        gc0 = [self.lon[i0], self.lat[i0]]
        __, __, dist = self.g.inv(self.lon[i0], self.lat[i0], self.lon[i2], self.lat[i2])
        lon1, lat1, __ = self.g.fwd(self.lon[i0], self.lat[i0], faz, dist)
        gc1 = [lon1, lat1]
        return gc0, gc1

    def __calc_fp_corner_coords(self, actr_fp='default'):
        """
        Calculate the corner coordinates from the center coordinates
        Keywords:
        actr_fp : float (optional)
            'default'| Custom across track footprint value in meter
        """
        # Set custom across track footprint (if required)
        if np.isreal(actr_fp):
            self.__actr_fp = actr_fp # Custom across track footprint
        else:
            self.__actr_fp = 1650.0  # Standard across track footprint
        self.corner_pos = np.ndarray(shape=(2, 4, self.n_records), dtype=np.float64)
        for i in np.arange(self.n_records):
            self.corner_pos[:,:,i] = geobox(self.lon[i], self.lat[i], # center position
                                            self.faz[i],              # forward azimuth
                                            self.__actr_fp,           # across track fp
                                            self.pointdist[i])        # along track fp
        # Unwrap in defined corner positions
        l1 = ['lon', 'lat']
        l2 = ['ur','ul', 'lr', 'll']
        for i in [0, 1]:
            for j  in [0, 1, 2, 3]:
                setattr(self, l1[i]+'_'+l2[j], self.corner_pos[i,j,:])


    def __calc_pointdist(self):
        """
        Distance in m from ith point to ith+1 point
        (will be used as along-track footprint)
        """
        self.pointdist = np.ndarray(shape=self.n_records, dtype=np.float64)
        self.faz = np.ndarray(shape=self.n_records, dtype=np.float64)
        self.baz = np.ndarray(shape=self.n_records, dtype=np.float64)
        for i in np.arange(self.n_records-1):
            self.faz[i], self.baz[i], self.pointdist[i] = self.__geo_inv(self.lon[i], self.lat[i],
                                                                         self.lon[i+1], self.lat[i+1])
        # Assume contant value for last point
        self.pointdist[-1] = self.pointdist[-2]
        self.baz[-1] = self.baz[-2]
        self.faz[-1] = self.faz[-2]

    def __geo_inv(self, lon1, lat1, lon2, lat2):
        faz, baz, dist = self.g.inv(lon1, lat1, lon2, lat2)
        return faz, baz, dist

    def __dist_filter(self):
        # Filter bad values
        invalid_list = np.where(self.pointdist > self._pointdist_tresh)[0]
        self.pointdist[invalid_list] = self._pointdist_tresh

    def to_CCfile(self, folder):
        """
        Write calculated corner coordinates to file in ASCII format for later re-use

        Arguments:
        folder : str
            Output folder (filename is automatically generated)
        """
        str_fmt_dt_filename = '%Y%m%dT%H%M%S'
        str_fmt_filename = "cs2_reforbit_%06g_%s_%s.dat"
        str_fmt_dt_rec = '%Y %m %d %H %M %S'
        str_fmt_rec = '%s.%03d'+' %14.8f'*10+'\n'
        # Create filename
        start_time = self.time[0].strftime(str_fmt_dt_filename)
        stop_time = self.time[-1].strftime(str_fmt_dt_filename)
        filename = str_fmt_filename% (self.orbit, start_time, stop_time)
        fileout = os.path.join(folder, filename)
        # Write data to file
        fhandle = open(fileout,'w')
        for i in np.arange(len(self.time)):
            fhandle.write(str_fmt_rec % (self.time[i].strftime(str_fmt_dt_rec),
                                         self.time[i].microsecond/1000,
                                         self.lon[i], self.lat[i],
                                         self.lon_ur[i], self.lat_ur[i],
                                         self.lon_ul[i], self.lat_ul[i],
                                         self.lon_lr[i], self.lat_lr[i],
                                         self.lon_ll[i], self.lat_ll[i] ))
        fhandle.close()


class CS2OrbitResData(object):
    """
    pyoval.CS2OrbitResData : Class for the Resampling of
    arbitrary validation data on CryoSat-2 ground tracks
    """

    def __init__(self):
        self.__flag_orbit_pos = False
        self.__flag_calval_data = False
        self.__flag_pdf_bins = False
        self.g = pyproj.Geod(ellps='WGS84')
        self.data_label = 'datalabel'
        self.data_extlabel = 'data_extlabel'
        self.data_unit = 'data_unit'


    def add_calval_data(self, time, lon, lat, val):
        """
        Add the Cal/Val data to be resampled

        Arguments:
        time : datetime object (numpy) array, mandatory
            timestamps of each cal/val data point
        lon :  float or double (numpy) array, mandatory
            longitude of each cal/val data point
        lat : float or double (numpy) array, mandatory
            latitude of each cal/val data point
        val : float or double (numpy) array, mandatory
            data sample at lon, lat
        """
        self.cv_time = time
        self.cv_lon = lon
        self.cv_lat = lat
        self.cv_data = val
        self.__flag_calval_data = True


    def add_cs2orbit(self, reforbit):
        """Add information on the corner coordinates of the CryoSat-2 footprint"""
        if not  isinstance(reforbit, CS2RefOrbit):
            raise ValueError('Function argument needs to be pyoval.RefOrbit class')
        self.cs2fp = reforbit
        self.__flag_orbit_pos = True


    def set_data_type(self, datatype):
        """
        Use predefined defaults for pdf setting of different data types

        Freeboard PDF    0cm  -   50cm, bin size = 2cm
                        50cm  -  400cm, bin size = 10cm
                       400cm  - 1000cm, bin size = 50cm

        Thickness PDF    0cm  -  500cm, bin size = 10cm
                       500cm  - 2000cm, bin size = 50cm
        Arguments:
        datatype : str
            Valid choices: 'freeboard' | 'thickness'
        """
        if datatype.lower() == 'freeboard':
            self.set_pdf_bins(np.concatenate((np.linspace(0, 50, 50/2+1),
                                              np.linspace(60, 400, 350/10),
                                              np.linspace(450, 1000, 600/50)))*0.01)
        elif datatype.lower() == 'thickness':
            self.set_pdf_bins(np.concatenate((np.linspace(0, 500, 500/10+1),
                                              np.linspace(550, 2000, 1500/50)))*0.01)
        else:
            raise ValueError('Wrong datatype string:\nValid choices: \'freeboard\', \'thickness\'')


    def set_pdf_bins(self, bins):
        """
        Manually set bin locations of PDF

        Arguments:
        bins : value (numpy) array
            bin locations (see numpy histogram method for details)
        """
        self.bins = bins
        self._nbins = len(self.bins)-1
        self.bin_center = np.ndarray(shape=(self._nbins))
        for i in np.arange(self._nbins):
            self.bin_center[i] = (self.bins[i+1]-self.bins[i])/2.0 + self.bins[i]
        self.__flag_pdf_bins = True


    def set_data_label(self, label):
        """ Sets the short label (for filename) of the cal/val data type"""
        if not isinstance(label, str):
            raise ValueError('Argument needs to be of type string')
        self.data_label = label


    def set_data_extlabel(self, label):
        """ Sets the extended label (for preview figure) of the cal/val data type"""
        if not isinstance(label, str):
            raise ValueError('Argument needs to be of type string')
        self.data_extlabel = label


    def set_parameter_extlabel(self, label):
        """ Sets the extended label (for preview figure) of the cal/val data type"""
        if not isinstance(label, str):
            raise ValueError('Argument needs to be of type string')
        self.parameter_extlabel = label


    def set_data_unit(self, unit):
        """ Sets the unit (for preview figure) of the cal/val data type"""
        if not isinstance(unit, str):
            raise ValueError('Argument needs to be of type string')
        self.data_unit = unit


    def resample(self):
        """
        Perform the data resampling
        """
        # Check preconditions:
        if not self.__flag_orbit_pos:
            raise ValueError('Use method add_cs2orbit to add footprint loctions before resampling')
        if not self.__flag_calval_data:
            raise ValueError('Use method add_calval_data to add high resolution data before resampling')
        if not self.__flag_pdf_bins:
            raise ValueError('Use method set_data_type of set_pdf_bins to add binning information before resampling')

        # Prepare data fiels
        self.n_samples = np.zeros(shape=(self.cs2fp.n_records), dtype=np.int32)
        self.time_offset_seconds = np.ones(shape=(self.cs2fp.n_records), dtype=np.float32)*np.nan
        self.res_mean = np.zeros(shape=(self.cs2fp.n_records), dtype=np.float32)*np.nan
        self.res_median = np.zeros(shape=(self.cs2fp.n_records), dtype=np.float32)*np.nan
        self.res_sdev = np.zeros(shape=(self.cs2fp.n_records), dtype=np.float32)*np.nan
        self.res_mode = np.zeros(shape=(self.cs2fp.n_records), dtype=np.float32)*np.nan
        self.pdf_scaled = np.zeros(shape=(self.cs2fp.n_records, self._nbins), dtype=np.float32)
        self.pdf_scale_fact = np.zeros(shape=(self.cs2fp.n_records), dtype=np.float32)
        self.cs2fp_points_list = np.ndarray(shape=(self.cs2fp.n_records), dtype=object)

        # All clear, start resampling
        for i in np.arange(self.cs2fp.n_records):
            # Define CS-2 footprint as polygon
            cs2fp_x = [self.cs2fp.lon_ur[i], self.cs2fp.lon_ul[i], self.cs2fp.lon_ll[i], self.cs2fp.lon_lr[i]]
            cs2fp_y = [self.cs2fp.lat_ur[i], self.cs2fp.lat_ul[i], self.cs2fp.lat_ll[i], self.cs2fp.lat_lr[i]]
            # CryoSat-2 tracks are mostly norths-south in cal/val regions
            # -> Use a latitude filter to speed things up again
            # -> look only at values +/- 0.02 degree latitude around footprint center coordinate
            data_subset = np.where(np.abs(self.cs2fp.lat[i] - self.cv_lat) < 0.02)[0]

            # Check if any Cal/Val data point near footprint center coordinate
            if len(data_subset) == 0:
                continue
            # Simple Points in Polygon calculation
            # (assumuning linearity for lon/lat at footprint scale)
            cs2fp_subset_points = np.where(points_in_polygon(self.cv_lon[data_subset],
                                                             self.cv_lat[data_subset],
                                                             cs2fp_x, cs2fp_y))[0]
            # Sanity check: Number of points in footprint > 0
            if len(cs2fp_subset_points) == 0:
                continue
            cs2fp_points = data_subset[cs2fp_subset_points]

            # Sanity check: Not only NaN in the footprint points
            if len(np.where(np.isfinite(self.cv_data[cs2fp_points]))[0]) == 0:
                continue
            self.cs2fp_points_list[i] = cs2fp_points
            # Analyze statistics of footprint cal/val data
            self.n_samples[i] = len(cs2fp_points)
            self.res_mean[i] = np.nanmean(self.cv_data[cs2fp_points])
            self.res_median[i] = np.median(self.cv_data[cs2fp_points])
            self.res_sdev[i] = np.nanstd(self.cv_data[cs2fp_points])
            # Calculate maximum time difference
            time_delta = self.cs2fp.time[i] - self.cv_time[cs2fp_points]
            max_index = np.argmax(np.abs(time_delta))
            self.time_offset_seconds[i] = time_delta[max_index].total_seconds()
            # Calculate pdf
            # Scale the pdf result that maximum will be on
            # Number of points per bin = pdf_scaled * pdf_scale_fact
            freq, bins = np.histogram(self.cv_data[cs2fp_points], bins=self.bins, density=True)
            self.pdf_scale_fact[i] = np.amax(freq)
            self.pdf_scaled[i,:] = freq.astype(np.float32)/self.pdf_scale_fact[i]
            self.res_mode[i] = self.bin_center[np.argmax(freq)]

        # Get overlap between orbit and cal/cal data
        # (remove empty footprints at end & beginning of orbit)
        self._get_overlap()

    def to_file(self, folder):
        """
        Write resampled data to file in ASCII format

        Output format
        Time, Lon, Lat, Time offset, Number of samples, Mean, Median, Std.Dev., PDF scale factor, PDF Bins

        Arguments:
        folder : str
            Output folder (filename is automatically generated)
        """
        str_fmt_dt_filename = '%Y%m%dT%H%M%S'
        str_fmt_filename = "cs2_resdata_%s_%06g_%s_%s.dat"
        str_fmt_dt_rec = '%Y %m %d %H %M %S'
        str_fmt_rec = '%s.%03d'+' %14.8f'*2+' %9g %6g'+' %14.5f'*5

        # Create filename
        i0 = self.overlap[0]
        i1 = self.overlap[-1]
        start_time = self.cs2fp.time[i0].strftime(str_fmt_dt_filename)
        stop_time = self.cs2fp.time[i1].strftime(str_fmt_dt_filename)
        filename = str_fmt_filename% (self.data_label, self.cs2fp.orbit, start_time, stop_time)
        fileout = os.path.join(folder, filename)
        # Write data to file
        fhandle = open(fileout,'w')
        for i in np.where(self.overlap)[0]:
            line = str_fmt_rec % (self.cs2fp.time[i].strftime(str_fmt_dt_rec),
                                  self.cs2fp.time[i].microsecond/1000,
                                  self.cs2fp.lon[i], self.cs2fp.lat[i],
                                  self.time_offset_seconds[i],
                                  self.n_samples[i],
                                  self.res_mean[i],
                                  self.res_median[i],
                                  self.res_sdev[i],
                                  self.res_mode[i],
                                  self.pdf_scale_fact[i])

            pdf_str = " ".join('%12.6f'% bin_val for bin_val in self.pdf_scaled[i,:])
            line = " ".join([line, pdf_str])
            line = " ".join([line,'\n'])
            fhandle.write(line)
        fhandle.close()

    def summary_plot(self, folder, **kwargs):
        """ Create Overview Plot in specified folder"""
        plot = CS2OrbitResPlot(self)
        plot.save(folder, **kwargs)

    def _get_overlap(self):
        self.overlap = np.zeros(shape=(self.cs2fp.n_records), dtype=bool)
        overlap_index = np.where(self.n_samples > 0)[0]
        self.overlap[np.arange(overlap_index[0], overlap_index[-1]+1)] = True


class CS2OrbitResPlot(object):
    """ Plot class for summary plots of resampled data"""

    def __init__(self, data):

        # Style thingies
        rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        self.fig = plt.figure(figsize=(8.26771654, 11.6929134))
        plt.subplots_adjust(hspace=0.5, wspace=0.6)
        self.ax_map = plt.subplot2grid((7, 1), (0, 0), rowspan=4)
        self.ax_data = plt.subplot2grid((7, 1), (4, 0), rowspan=2)
        self.ax_tdiff = plt.subplot2grid((7, 1), (6, 0))


        self.g = pyproj.Geod(ellps='WGS84')
        self.data = data
        self.overlap = np.where(self.data.overlap)[0]
        self.dateformatter = DateFormatter('%H:%M:%S')
        self.dates_float = date2num(self.data.cs2fp.time[self.overlap])
        self.__get_filename()
        self.__get_proj_info()

        self._fontsize_small = 6
        self._fontsize_medium = 8
        self._fontsize_large = 12

        self._fontsize_map_dt_labels = 8
        self._fontsize_map_tick_labels = 6

        self._color_map_wet = '#CAE1FF'  # Ocean Color
        self._color_map_dry = '#FFF8DC'  # Land Color
        self._color_map_cl = '#AAAAAA'   # Coast Line
        self._color_bg = '#e5e5e5'

        # Color map for number of samples
        self._cmap_nsamples = make_cmap([(255,99,71), (255,215,0), (127,255,0)],
                                         lut=3, bit=True)
        self._cmap_nsamples.set_under('#AAAAAA')
        self._cmap_nsamples.set_over('#00FF00')

        # Start plot
        self.__projdata()
        self.__plot_time_delta()
        self.__plot_data()
        self.__plot_map()
        self.__add_labels()


    def save(self, folder, dpi=300, bbox_inches=0.0, pad_inches=0.00, auto_crop=False):
        fullpath = os.path.join(folder, self.filename)
        if auto_crop:
            bbox_inches='tight'
            pad_inches = 0.05
        try:
            plt.savefig(fullpath, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
        except IOError: # pdf open and cannot be overwritten?
            fileName, fileExtension = os.path.splitext(self.filename)
            for i in np.arange(10):
                filename_conflict = fileName+'_conflict'+"%02g"%i+fileExtension
                fullpath_conflict = os.path.join(folder, filename_conflict)
                if not os.path.exists(fullpath_conflict):
                    break
            plt.savefig(fullpath_conflict, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
        # Free Memory
        # Debug Purposes only: plt.show(block=True)
        plt.close(self.fig)


    def __get_filename(self):
        str_fmt_dt_filename = '%Y%m%dT%H%M%S'
        str_fmt_filename = "cs2_resdata_%s_%06g_%s_%s.pdf"
        i0 = self.data.overlap[0]
        i1 = self.data.overlap[-1]
        start_time = self.data.cs2fp.time[i0].strftime(str_fmt_dt_filename)
        stop_time = self.data.cs2fp.time[i1].strftime(str_fmt_dt_filename)
        self.filename = str_fmt_filename % (self.data.data_label, self.data.cs2fp.orbit, start_time, stop_time)


    def __get_proj_info(self):
        self.lon_0 = np.median(self.data.cs2fp.lon[self.overlap])
        self.lat_0 = np.median(self.data.cs2fp.lat[self.overlap])
        self.lat_ts = self.lat_0
        faz, baz, length = self.g.inv(self.data.cs2fp.lon[self.overlap[0]], self.data.cs2fp.lat[self.overlap[0]],
                                      self.data.cs2fp.lon[self.overlap[-1]], self.data.cs2fp.lat[self.overlap[-1]])
        self.width = length*1.02
        self.height = length*1.02
        self.proj = 'stere'
        self.res = 'h'

        self.basemap = Basemap(width=self.width, height=self.height,
                         resolution=self.res, projection=self.proj,
                         lat_ts=self.lat_ts, lat_0=self.lat_0 ,lon_0=self.lon_0,
                         ax=self.ax_map)

        # Angle of orbit in projection coordinates -> Orientation of time labels
        x0, y0 = self.basemap(self.data.cs2fp.lon[self.overlap[0]], self.data.cs2fp.lat[self.overlap[0]])
        x1, y1 = self.basemap(self.data.cs2fp.lon[self.overlap[-1]], self.data.cs2fp.lat[self.overlap[-1]])
        self._orbit_pj_angle = np.rad2deg(np.arctan((y1-y0)/(x1-x0)))


    def __projdata(self):

        # Calculate projection coordinates for map plotting
        self.ul_pj_x, self.ul_pj_y = self.basemap(self.data.cs2fp.lon_ul, self.data.cs2fp.lat_ul)
        self.ur_pj_x, self.ur_pj_y = self.basemap(self.data.cs2fp.lon_ur, self.data.cs2fp.lat_ur)
        self.ll_pj_x, self.ll_pj_y = self.basemap(self.data.cs2fp.lon_ll, self.data.cs2fp.lat_ll)
        self.lr_pj_x, self.lr_pj_y = self.basemap(self.data.cs2fp.lon_lr, self.data.cs2fp.lat_lr)
        self.cr_pj_x, self.cr_pj_y = self.basemap(self.data.cs2fp.lon, self.data.cs2fp.lat)


    def __plot_map(self):

        plt.sca(self.ax_map)
        self.ax_map.set_axis_bgcolor(self._color_map_wet)
        self.basemap.drawmapboundary(linewidth=0.2, zorder=1000)
        self.basemap.drawcoastlines(color=self._color_map_cl, linewidth=0.1, zorder=1000)
        self.basemap.fillcontinents(color=self._color_map_dry)

        # Plot footprints color coded for number of points
        fp_patches = []
        for i in self.overlap:
            xy_pj = [[self.ul_pj_x[i], self.ul_pj_y[i]],
                     [self.ur_pj_x[i], self.ur_pj_y[i]],
                     [self.lr_pj_x[i], self.lr_pj_y[i]],
                     [self.ll_pj_x[i], self.ll_pj_y[i]]]
            polygon = Polygon(xy_pj, True)
            fp_patches.append(polygon)
        cs2fp_faces = PatchCollection(fp_patches, cmap=self._cmap_nsamples,
                                      edgecolors='none', zorder=900)#, rasterized=True)


        color_id = np.floor(np.log10(self.data.n_samples[self.overlap]+0.1))+0.5
        cs2fp_faces.set_array(color_id)
        cs2fp_faces.set_clim(vmin=0, vmax=3)
        self.ax_map.add_collection(cs2fp_faces)

        # Label time steps on track (according to xticks in data and time_delte plots)
        dt_ticks = self.ax_data.xaxis.get_majorticklocs()
        for dt_tick in dt_ticks:
            index = np.argmin(np.abs(self.dates_float - dt_tick))
            index = self.overlap[index]
            dt_str = num2date(dt_tick).strftime('%H:%M:%S')
            plt.gca().annotate(dt_str, (self.cr_pj_x[index], self.cr_pj_y[index]),
                    xytext=(5, 5), textcoords='offset points', ha='left', va='center',
                    rotation=self._orbit_pj_angle+90, fontsize=self._fontsize_map_dt_labels)

        self.__add_grid(latstep=0.025, lonstep=0.25, linewidth=0.05,  color='#ECFAFF')
        self.__add_grid(latlabel=True, lonlabel=True)
        # Shift map plot to the left
        a = self.ax_map.get_position().bounds
        a = np.array(a)
        a[0] = 0.05
        self.ax_map.set_position(a)
        # Add the colorbar
        axins = inset_axes(self.ax_map,
                           width="3%", # width = 10% of parent_bbox width
                           height="50%", # height : 50%
                           loc=3,
                           bbox_to_anchor=(1.02, 0., 1, 1),
                           bbox_transform=self.ax_map.transAxes,
                           borderpad=0,
                           )
        cb = plt.colorbar(cs2fp_faces, cax=axins, extend='min',
                          ticks=[0,1,2,3], extendrect=True)
        cb.ax.tick_params(labelsize=self._fontsize_map_tick_labels)
        cb.solids.set_edgecolor("1.0")
        cb.outline.set_linewidth(0.5)
        cb.outline.set_alpha(0.0)
        for t in cb.ax.get_yticklines(): t.set_color("1.0")
        cb.ax.set_yticklabels(['1', '10', '100', '1000+ Samples'])
        cb.ax.tick_params('both', length=0.1, which='major', pad=-5)

    def __plot_data(self):

        plt.sca(self.ax_data)
        # Calculate quadriliteratl coordinates for pcolor represenations of the bins
        step = (self.dates_float[-1] - self.dates_float[0])/(len(self.dates_float))
        n = len(self.overlap)
        px = np.ndarray(shape=n+1)
        steps = np.ndarray(shape=n+1)
        ref_step = self.dates_float[1]-self.dates_float[0]
        for i in np.arange(n):
            if i < n-1:
                step = self.dates_float[i+1]-self.dates_float[i]
            if step > 2.0 * ref_step:  # Catch some data gaps in CS-2 data (SAR | SARin)
                step = ref_step
            px[i] = self.dates_float[i]-step/2.
            steps[i]= step
        px[-1] = px[-2]+step
        py = self.data.bins
        xx = np.array([px,]*len(py)).transpose()
        yy = np.array([py,]*len(px))
        # Plot bins in back ground
        pdf_scaled = self.data.pdf_scaled[self.overlap,:]
        cmap_pdf = plt.get_cmap('hot_r')
        cmap_pdf.set_under(self._color_bg)
        self.ax_data.pcolor(xx, yy, pdf_scaled, cmap=cmap_pdf, vmin=1e-8, vmax=1.2,
                            edgecolor='none', rasterized=True, zorder=50 )
        # Plot statistical Parameters as symbol overlay
        self.ax_data.plot_date(self.data.cs2fp.time[self.overlap],
                               self.data.res_mean[self.overlap], 'o',
                               markeredgewidth=0.2, markersize=2,
                               markerfacecolor='none', markeredgecolor='0.0',
                               label='Mean',
                               zorder=100, rasterized=True)
        self.ax_data.plot_date(self.data.cs2fp.time[self.overlap],
                               self.data.res_mode[self.overlap], '^',
                               markeredgewidth=0.2, markersize=2,
                               markerfacecolor='none', markeredgecolor='0.0',
                               label='Mode', zorder=100, rasterized=True)

        self.ax_data.plot_date(self.data.cs2fp.time[self.overlap],
                               self.data.res_median[self.overlap], 'x',
                               markeredgewidth=0.2, markersize=2,
                               markerfacecolor='none', markeredgecolor='0.0',
                               label='Median', zorder=100, rasterized=True)

        self.ax_data.vlines(self.dates_float,
                            self.data.res_mean[self.overlap]-0.5*self.data.res_sdev[self.overlap],
                            self.data.res_mean[self.overlap]+0.5*self.data.res_sdev[self.overlap],
                            linewidth=0.2, zorder=100, label='Standard Deviation')


        self.ax_data.hlines(self.data.bins,
                            self.dates_float[0], self.dates_float[-1],
                            linewidth=0.2, zorder=90, color='1.0', alpha=0.5)

        plt.legend(bbox_to_anchor=(0.2, 1.01, 0.6, .102), loc=3,
                    ncol=4, mode="expand", borderaxespad=0., fontsize=self._fontsize_small,
                    frameon=False, markerscale=2, scatterpoints=1)

        self.ax_data.set_ylabel(self.data.data_unit)

        self.__axis_style(self.ax_data)
        self.__axis_label(self.ax_data, self.data.parameter_extlabel)

    def __plot_time_delta(self):

        plt.sca(self.ax_tdiff)

        self.ax_tdiff.plot_date(self.data.cs2fp.time[self.overlap],
                                self.data.time_offset_seconds[self.overlap], '.',
                                markersize=2, markeredgecolor=None, color='#d02090', zorder=100)

        # Replace yticks with custom strings
        yticks = self.ax_tdiff.yaxis.get_majorticklocs()
        ytick_str = []
        for ytick in yticks: ytick_str.append(timedelta_str(datetime.timedelta(seconds=ytick)))
        self.ax_tdiff.set_yticklabels(ytick_str)


        self.__axis_style(self.ax_tdiff)
        self.__axis_label(self.ax_tdiff, u'Temporal Offset')



    def __add_grid(self, latstep=0.1, lonstep=2, latlabel=False, lonlabel=False,
                         linewidth=0.2, color='#FFFFFF'):

        latlabels = [0,0,0,0]
        lonlabels = [0,0,0,0]
        if latlabel: latlabels = [1,0,0,0]
        if lonlabel: lonlabels = [0,0,0,1]

        self.basemap.drawparallels(np.arange(70, 88, latstep), labels=latlabels, color=color,
                                   dashes=[], linewidth=linewidth,
                                   fontsize=self._fontsize_map_tick_labels)

        self.basemap.drawmeridians(np.arange(0, 360, lonstep), labels=lonlabels,
                                   latmax=88, color=color,
                                   dashes=[], linewidth=linewidth,
                                   fontsize=self._fontsize_map_tick_labels)


    def __axis_style(self, ax):

        ax.patch.set_facecolor(self._color_bg)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_major_formatter(self.dateformatter)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=self._fontsize_medium, width=0.75, color='0.5', length=2)
        ax.xaxis.label.set_color('0.5')
        ylocs = ax.yaxis.get_majorticklocs()
        yminorticks = MultipleLocator( (ylocs[1]-ylocs[0]) / 2.0 )
        ax.yaxis.set_minor_locator(yminorticks)
        ax.grid(True, 'major', color='1.0', linestyle='-', linewidth=0.75, zorder=900)
        ax.grid(True, 'minor', color='1.0', linestyle='-', linewidth=0.25, zorder=900, alpha=0.5)
        for t in ax.xaxis.get_ticklabels():
            t.set_color('0.5')
        for t in ax.yaxis.get_ticklabels():
            t.set_color('0.5')
        #remove axis border
        for child in ax.get_children():
            if isinstance(child, Spine):
                child.set_alpha(0)
        # restyle the tick lines
        for line in ax.get_xticklines() + ax.get_yticklines():
            line.set_markersize(5)
            line.set_color("0.5")
            line.set_markeredgewidth(0.8)
        #remove the minor tick lines
        for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
            line.set_markersize(0)
        #only show bottom left ticks, pointing out of axis
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    def __axis_label(self, ax, label):

        # Get location in figure coordinates
        fig_coord = ax.get_position()

        rect = Rectangle((fig_coord.x0+fig_coord.width+0.001, fig_coord.y0),
                          0.02, fig_coord.height,
                          facecolor='#aaaaaa', edgecolor='none',
                          transform=self.fig.transFigure, alpha=0.5, zorder=100)
        self.fig.patches.append(rect)

        plt.annotate(unicode(label),
                     (fig_coord.x0+fig_coord.width+0.01, fig_coord.y0+0.5*fig_coord.height),
                     xycoords='figure fraction',
                     ha='center', va='center', rotation=-90,
                     fontsize=8, fontstretch='extra-expanded', fontweight='bold', zorder=900 )


    def __add_labels(self):
        # Document Title
        plt.annotate('CryoSat-2 Sea-Ice Validation Data Set',
                     (0.5, 0.95), xycoords='figure fraction',
                     ha='center', va='center',
                     fontsize=18)
        # Data Extended Label
        plt.annotate(self.data.data_extlabel,
                     (0.5, 0.92), xycoords='figure fraction',
                     ha='center', va='center',
                     fontsize=14)
        i0 = self.overlap[0]
        i1 = self.overlap[-1]
        void, void, dist_val = self.g.inv(self.data.cs2fp.lon[i0], self.data.cs2fp.lat[i0],
                                          self.data.cs2fp.lon[i1], self.data.cs2fp.lat[i1])
        labels = ['Orbit', 'Day', 'Start Time', 'End Time', 'Length']
        values = [unicode(self.data.cs2fp.orbit),
                  self.data.cs2fp.time[self.overlap[0]].strftime("%Y-%m-%d"),
                  self.data.cs2fp.time[self.overlap[0]].strftime('%H:%M:%S %Z'),
                  self.data.cs2fp.time[self.overlap[-1]].strftime('%H:%M:%S %Z'),
                  unicode("%5.1f km" % (1e-3 * dist_val)) ]
        lx = 0.77
        ly = 0.89
        dly = -0.04
        i = 0
        for label, value in zip(labels, values):
            plt.annotate(label,
                         (lx, ly+i*dly), xycoords='figure fraction',
                         ha='left', va='center',
                         fontsize=7, color="0.6", fontweight='bold' )
            plt.annotate(value,
                         (lx, ly+i*dly-0.013), xycoords='figure fraction',
                         ha='left', va='center',
                         fontsize=11)
            i += 1
        # Project Acknowledgement
        plt.annotate(u'Created for ESA CryoVal-SI (AO/1-­‐7340/12/NL/CT)',
                     (0.9, 0.05), xycoords='figure fraction',
                     ha='right', va='center',
                     fontsize=self._fontsize_small, fontstretch='extra-expanded', fontstyle='italic' )

class GreatCircleTrack():
    """
    Class for calculting across track errors (XTE) and
    along track distances (ATD) in spherical geometry
    Based on Ed Williams Aviation Formulary V1.43
    """
    def __init__(self, gcpos0, gcpos1):
        self.gcpos0 = gcpos0
        self.gcpos1 = gcpos1
        self.g = pyproj.Geod(ellps='WGS84')
        self.n_records = 0

    def add_track(self, lon, lat):
        """ Add track data (Arrays of longitude, latitude) """
        self.lon = lon
        self.lat = lat
        self.n_records = len(lon)
        self.r = self.__geo_getr(self.lat)

    def project(self):
        """ Calculates XTE and ATD for track data with respect to great circle """
        self.__geo_commons()
        self.__geo_xte()
        self.__geo_atd()
        self.__geo_gcp()

    def __geo_commons(self):
        """ Creates parameter arrays used for calculation of xte, atd & gcp """
        self.crs_ab = self.__geo_crs(self.gcpos0[0], self.gcpos0[1],
                                     self.gcpos1[0], self.gcpos1[1])
        self.dist_ad = np.ndarray(shape=(self.n_records), dtype=np.float64)
        self.crs_ad = np.ndarray(shape=(self.n_records), dtype=np.float64)
        self.uxte = np.ndarray(shape=(self.n_records), dtype=np.float64)
        for i in np.arange(self.n_records):
              self.dist_ad[i] = self.__geo_udst(self.gcpos0[0], self.gcpos0[1],
                                                self.lon[i], self.lat[i])
              self.crs_ad[i] = self.__geo_crs(self.gcpos0[0], self.gcpos0[1],
                                              self.lon[i], self.lat[i])
              self.uxte[i] = np.arcsin(np.sin(self.dist_ad[i]) * \
                                       np.sin(self.crs_ad[i]-self.crs_ab) )

    def __geo_xte(self):
          """ Calculate Across Track Error (xte) in meter"""
          self.xte = np.ndarray(shape=(self.n_records), dtype=np.float64)
          self.xte = self.uxte * self.__geo_getr(self.lat)

    def __geo_atd(self):
        """ Calculate Along Track Distance (atd) in meter """
        self.atd = np.ndarray(shape=(self.n_records), dtype=np.float64)
        self.eang = np.ndarray(shape=(self.n_records), dtype=np.float64)
        for i in np.arange(self.n_records):
            self.atd[i] = np.arccos(np.cos(self.dist_ad[i])/np.cos(self.uxte[i]))
            self.atd[i] *= self.r
            self.eang[i] = self.__geo_eang(self.crs_ab, self.crs_ad[i])
            if self.eang[i] > 2.0*np.pi: self.atd[i] *= -1.0

    def __geo_gcp(self):
        """
        Calculates the coordinates of the projected track
        points on the great circle
        """
        faz, __, dist = self.g.inv(self.gcpos0[0], self.gcpos0[1],
                                   self.gcpos1[0], self.gcpos1[1])
        tc = -1.0*np.deg2rad(faz)
        lon1 = np.deg2rad(self.gcpos0[0])
        lat1 = np.deg2rad(self.gcpos0[1])
        self.gcp = np.ndarray(shape=(self.n_records, 2), dtype=np.float64)
        for i in np.arange(self.n_records):
            d = self.atd[i]/self.r
            lat = np.arcsin(np.sin(lat1)*np.cos(d)+np.cos(lat1)*np.sin(d)*np.cos(tc))
            dlon = np.arctan2(np.sin(tc)*np.sin(d)*np.cos(lat1),np.cos(d)-np.sin(lat1)*np.sin(lat))
            lon = np.mod( lon1-dlon + np.pi, 2.0*np.pi ) - np.pi
            self.gcp[i,:] = [np.rad2deg(lon), np.rad2deg(lat)]

    def __geo_udst(self, deg_lon1, deg_lat1, deg_lon2, deg_lat2):
        """ Distance on sphere with unity radius """
        lon1 = np.deg2rad(deg_lon1)
        lon2 = np.deg2rad(deg_lon2)
        lat1 = np.deg2rad(deg_lat1)
        lat2 = np.deg2rad(deg_lat2)
        d = 2.0 * np.arcsin(np.sqrt( (np.sin((lat1-lat2)/2.0))**2.0 +\
            np.cos(lat1)*np.cos(lat2)*(np.sin((lon1-lon2)/2.0))**2.0) )
        return d

    def __geo_crs(self, deg_lon1, deg_lat1, deg_lon2, deg_lat2):
        """ Heading between two points in radians """
        lon1 = np.deg2rad(deg_lon1)
        lon2 = np.deg2rad(deg_lon2)
        lat1 = np.deg2rad(deg_lat1)
        lat2 = np.deg2rad(deg_lat2)
        crs = np.arctan2(np.sin(lon1-lon2)*np.cos(lat2), \
                         np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon1-lon2))
        crs = crs % ( 2.0 * np.pi )
        return crs

    def __geo_getr(self, lat):
        """ Mean earth ellipsoidial radius at latitude range in meter """
        sma  = 6378137.0
        smn  = 6356752.3142755
        mlat = np.nanmean(lat)
        r = np.sqrt( ((sma**2.0 / np.cos(mlat))**2.0 + (smn**2.0 / np.sin(mlat))**2.0) / \
                     ((sma / np.cos(mlat))**2.0 + (smn / np.sin(mlat))**2.0) )
        return r

    def __geo_eang(self, crs1, crs2):
        """ Enclosed angle between two headings in radian """
        eang = np.abs(crs2-crs1)
        if eang > np.pi: eang = 2.0*np.pi - eang
        return eang

    def debug_map(self):
        """ Create a quick & dirty verification map """
        lat_0 = np.median(self.lat)
        lon_0 = np.median(self.lon)
        plt.figure("class GreatCircleTrack - Debug Map")
        m = Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0, lat_ts=lat_0,
                    resolution = 'i', width=5e5, height=5e5)
        m.drawcoastlines(color="black", linewidth=0.1)
        m.fillcontinents(color='#AAAAAA')
        m.drawmapboundary()
        xc, yc = m(self.lon, self.lat)
        x0, y0 = m(self.gcpos0[0], self.gcpos0[1])
        x1, y1 = m(self.gcpos1[0], self.gcpos1[1])
        xp, yp = m(self.gcp[:,0], self.gcp[:,1])
        m.scatter(xc, yc, s=20, color='blue', edgecolors='none', zorder=100)
        m.scatter(xp, yp, s=20, color='red', edgecolors='none', zorder=100)
        m.scatter(x0, y0, s=60, color='black', edgecolors='none', zorder=100)
        m.scatter(x1, y1, s=60, color='orange', edgecolors='none', zorder=100)
        for i in np.arange(self.n_records):
            m.plot([xc[i], xp[i]], [yc[i], yp[i]], color='black', zorder=99)
        m.drawgreatcircle(self.gcpos0[0], self.gcpos0[1],
                          self.gcpos1[0], self.gcpos1[1],
                          del_s=10, color='violet', lw=2,
                          alpha=0.5, zorder=98)
        plt.show()

def geobox(clon, clat, azimuth, width, height):
    """
    Calculates geopgraphical corner coordinates of a box
    from the center coordinates, the azimuth and width & height

                  azimuth
                | /
    ul          |/          ur
    *-----------------------*   -
    |                       |
    |           C           |   height
    |                       |
    *-----------------------*   -
    ll                      lr

    |        width          |

    Arguments:

    clon : float
        Center point longitude
    clat : float
        Center point latitude
    azimuth : float
        orientation of box: north in 'height'-direction
    width : float
        Width of box in meter
    height : float
        height of box in meter

    Returns:
      pos : array(2, 4)
          longitude and latiude of corner coordinates
          lon, lat for ur, ul, lr, ll
    """
    g = pyproj.Geod(ellps='WGS84')
    radius = np.sqrt((width/2.0)**2.0 + (height/2.0)**2.0)
    phi = np.degrees(np.arctan(height/width))
    alpha = np.zeros(4, dtype=np.float32)
    alpha[0] =  90.0 - phi # upper right
    alpha[1] = 270.0 + phi # upper left
    alpha[2] =  90.0 + phi # lower right
    alpha[3] = 270.0 - phi # lower left
    pos = np.zeros((2, 4), dtype=np.float32)
    for i in np.arange(4):
        pos[0, i], pos[1, i], baz = g.fwd(clon, clat, alpha[i]+azimuth, radius)
    return pos


def points_in_polygon(x, y, px, py, return_type='bool'):
    """
    Calculates which points x, y are inside a polygon, defined by px, py
    Arguments:
        x, y : float/double array
            positions of points to be tested
        px, py : float/double array with minimum of 3 points
    Keywords:
        return_type : str
            Indicates the type of return value
            'bool' (default) - Array of type bool of same length as x/y
                               True: in polygon, False: outside polygon
            'masked' - Masked numpy array (x, y)
    Return:
        Depending on keyword return either bool array of masked (x,y) arrays
    """
    from shapely.geometry import Polygon
    from shapely.geometry import Point

    # Sanity checks
    if np.shape(x) != np.shape(y):
        raise ValueError('Position of input points are of different shape')
    if np.shape(px) != np.shape(py):
        raise ValueError('Position of polygon points are of different shape')
    if len(px) < 3:
        raise ValueError('Not enough points to span polygon (px, py)')
    n_points = len(x)
    n_polygon_points = len(px)
    # Set up the polygon
    px_closed = np.concatenate((px, [px[0]]))
    py_closed = np.concatenate((py, [py[0]]))
    polygon_points = np.ndarray(shape=(n_polygon_points+1, 2), dtype=np.float32)
    for j in np.arange(n_polygon_points+1):
        polygon_points[j, 0] = px_closed[j]
        polygon_points[j, 1] = py_closed[j]
    polygon = Polygon(polygon_points)
    # Check if each point is in Polygon
    in_polygon = np.ndarray(shape=(n_points), dtype=np.bool)
    for i in np.arange(n_points):
        in_polygon[i] = polygon.contains(Point(x[i], y[i]))
    if return_type == 'bool':
        return in_polygon
    elif return_type == 'masked':
        return (np.ma.masked_where(~in_polygon, x),
                np.ma.masked_where(~in_polygon, y))
    else:
        raise ValueError('return_type '+str(return_type)+' unkown. \'bool\' (default) or \'masked\'')


def make_cmap(colors, position=None, bit=False, lut=256):
    '''
	NAME
		Custom Colormaps for Matplotlib
	PURPOSE
		This program shows how to implement make_cmap which is a function that
		generates a colorbar.  If you want to look at different color schemes,
		check out https://kuler.adobe.com/create.
	PROGRAMMER(S)
		Chris Slocum
	REVISION HISTORY
		20130411 -- Initial version created
		20140313 -- Small changes made and code posted online
		20140320 -- Added the ability to set the position of each color

    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, lut)
    return cmap
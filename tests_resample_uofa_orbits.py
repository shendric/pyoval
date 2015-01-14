# -*- coding: utf-8 -*-
import os
import glob


# Import Data Resampling Class from user library
from lib.pyoval import CS2RefOrbit
from lib.pyoval import CS2OrbitResData
from lib.data import read_embird_dat
from lib.data import read_OIB
from lib.helpers import EmptyObject


class PathInfo():
    """ Definition of pathes and filenames for tests """
    def __init__(self):
        # Path to reference orbit files
        self.dir = EmptyObject()
        self.dir.workspace = '/Volumes/Data/OneDrive/CryoVal-SI/GoldenDays/'
        self.dir.SampleData = os.path.join(self.dir.workspace,'20130424')
        self.dir.RefOrbit = self.dir.SampleData
        self.dir.RefCCOrbit = self.dir.SampleData
        self.dir.SampleFigures = self.dir.SampleData
        # Files
        self.file = EmptyObject()
        #self.file.example_orbit = r'Location_file_21011_20140326T090003_20140326T102353.txt'
        #self.file.example_des_orbit = r'Location_file_10520_20120402T110853_20120402T124403_L1B_vBC.txt'
        #self.file.example_aem = r'20140321_allfinal.dat'
        self.file.example_oib = r'OIB_20130424_IDCSI2.txt'

def example_calc_corner_coordinates(orbit=10520,lon_limit=[-180.1, 180.1], lat_limit=[60., 88.]):
    """Calculate Corner Coordinates from reference orbit files"""
    # Read examplary reference orbit data (center locations of footprint)
    # Initialize RefOrbit Class
    # Required for setting projection, initial ROI etc.
    info = PathInfo()
    reforbit = CS2RefOrbit()

    # Limit reference orbit data to region of interest
    # (Speed things up for corner coordinate caluclation and resampling)
    # ! Needs to be called before from_file
    # This example: Lincoln Sea (CryoVex 2011)
    # Note: longitude limits in range between -180 and 180
    reforbit.set_sin_detection(False)
    reforbit.limit_region(lat_limit=lat_limit, lon_limit=lon_limit)
    reforbit_file= str(os.path.join(info.dir.RefOrbit, info.file.example_orbit))
    print reforbit_file
    # Read center coordinates (files provided from UCL for CryoVal-SI)
    # Input parameter: Full file name
    reforbit.from_file(reforbit_file)
    reforbit.debug_map()
    # Save corner coordinates to file (in specified folder)
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)
    
def example_resample_oib(orbit=10520):
    """
    Example of how to resample arbitrary validation data onto
    CryoSat-2 footprints (no drift correction applied)
    """

    #===================================================================================
    # 1) Get Corner coordinates
    #===================================================================================
    # Read Reference coorner coordinates from ascii file
    # Retrieve corner coordinates from orbit 5399
    # (Example, can be customized)
    info = PathInfo()
    cc_file = glob.glob(os.path.join(info.dir.RefCCOrbit, '*'+str(orbit)+'*'))[0]

    # Initialize reference orbit data object and read corner coordinate file
    reforbit = CS2RefOrbit()
    reforbit.from_file(cc_file, corner_coords=True)

    # Example: Read EM-Bird data
    #   Data structure of Cal/Val data can be arbitrary as long as
    #       * data is organized as vectors/array with identical shape
    #       * arrays exist for time, longitude, latitude, value (e.g. thickness, freeboard, ....)
    #       * time information is available as python datetime object (time zone UTC())
    oib_file = os.path.join(info.dir.SampleData, info.file.example_oib)
    print oib_file
    oib = read_OIB(oib_file)

    #===================================================================================
    # 2) Resample Cal/Val Data to CryoSat-2 orbits
    #===================================================================================
    # Initialize object
    # Required for projection and sanity checks
    resdata = CS2OrbitResData()

    # Add the Orbit Corner Coordinates
    resdata.add_cs2orbit(reforbit)

    # Add the Cal/Val data
    resdata.add_calval_data(oib.dt_time,   # time info as datetime object vector
                            oib.longitude, # longitude vector
                            oib.latitude,  # latitude vector
                            oib.freeboard) # cal/val data vector

    # make use of pre-defined bins from progress meeting 1
    # user bins can be added by the method pyoval.CS2OrbitResData.set_pdf_bins(bins)
    resdata.set_data_type('freeboard')

    # Resample Cal/Val Data
    resdata.resample()

    #===================================================================================
    # 3) Generate output (data files and summary plots)
    #===================================================================================
    # Define Labels etc for output
    resdata.set_data_label('oibql')
    resdata.set_data_extlabel("OIB freeboard")
    resdata.set_parameter_extlabel('Freeboard')
    resdata.set_data_unit('Meter')

    # Write ASCII data output
    resdata.to_file(info.dir.SampleFigures)

    # Create summary plot (in specified folder)
    resdata.summary_plot(info.dir.SampleFigures)

def example_resample_embird(orbit=10520):
    """
    Example of how to resample arbitrary validation data onto
    CryoSat-2 footprints (no drift correction applied)
    """

    #===================================================================================
    # 1) Get Corner coordinates
    #===================================================================================
    # Read Reference coorner coordinates from ascii file
    # Retrieve corner coordinates from orbit 5399
    # (Example, can be customized)
    info = PathInfo()
    cc_file = glob.glob(os.path.join(info.dir.RefCCOrbit, '*'+str(orbit)+'*'))[0]

    # Initialize reference orbit data object and read corner coordinate file
    reforbit = CS2RefOrbit()
    reforbit.from_file(cc_file, corner_coords=True)

    # Example: Read EM-Bird data
    #   Data structure of Cal/Val data can be arbitrary as long as
    #       * data is organized as vectors/array with identical shape
    #       * arrays exist for time, longitude, latitude, value (e.g. thickness, freeboard, ....)
    #       * time information is available as python datetime object (time zone UTC())
    aem_file = os.path.join(info.dir.SampleData, info.file.example_aem)
    aem = read_embird_dat(aem_file)

    #===================================================================================
    # 2) Resample Cal/Val Data to CryoSat-2 orbits
    #===================================================================================
    # Initialize object
    # Required for projection and sanity checks
    resdata = CS2OrbitResData()

    # Add the Orbit Corner Coordinates
    resdata.add_cs2orbit(reforbit)

    # Add the Cal/Val data
    resdata.add_calval_data(aem.dt_time,   # time info as datetime object vector
                            aem.longitude, # longitude vector
                            aem.latitude,  # latitude vector
                            aem.thickness) # cal/val data vector

    # make use of pre-defined bins from progress meeting 1
    # user bins can be added by the method pyoval.CS2OrbitResData.set_pdf_bins(bins)
    resdata.set_data_type('thickness')

    # Resample Cal/Val Data
    resdata.resample()

    #===================================================================================
    # 3) Generate output (data files and summary plots)
    #===================================================================================
    # Define Labels etc for output
    resdata.set_data_label('aem')
    resdata.set_data_extlabel("EM-Bird sea-ice thickness")
    resdata.set_parameter_extlabel('Sea-Ice Thickness')
    resdata.set_data_unit('Meter')

    # Write ASCII data output
    resdata.to_file(info.dir.SampleFigures)

    # Create summary plot (in specified folder)
    resdata.summary_plot(info.dir.SampleFigures)


def test_asscending_orbit():
    """ Test of SARin position corrections for descending orbits """
    info = PathInfo()
    reforbit = CS2RefOrbit()
    reforbit.set_sin_detection(False)
    reforbit.limit_region(lat_limit=[82., 88.], lon_limit=[-120., -40.])
    reforbit.from_file(os.path.join(info.dir.RefOrbit, info.file.example_orbit))
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)


def test_descending_orbit():
    """ Test of SARin position corrections for descending orbits """
    info = PathInfo()
    reforbit = CS2RefOrbit()
    reforbit.set_sin_detection(False)
    reforbit.limit_region(lat_limit=[82., 88.], lon_limit=[-120., -40.])
    reforbit.from_file(os.path.join(info.dir.RefOrbit, info.file.example_des_orbit))
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)


if __name__ == '__main__':
    #example_calc_corner_coordinates(orbit=21011,
    #                                lat_limit=[60., 88.], 
    #                                lon_limit=[-180.1, 180.1])
    #example_resample_embird(orbit=16139)
    #example_resample_embird(orbit=21077)
    example_resample_oib(orbit=16139)
    #example_resample_oib(orbit=21092)
    #test_descending_orbit()
    #test_asscending_orbit()
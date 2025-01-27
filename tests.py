# -*- coding: utf-8 -*-
import os
import glob


# Import Data Resampling Class from user library
from lib.pyoval import CS2RefOrbit
from lib.pyoval import CS2OrbitResData
from lib.data import read_embird
from lib.helpers import EmptyObject


class PathInfo():
    """ Definition of pathes and filenames for tests """
    def __init__(self):
        # Path to reference orbit files
        self.dir = EmptyObject()
        self.dir.workspace = os.path.dirname(__file__)
        self.dir.SampleData = os.path.join(self.dir.workspace, 'testdata')
        self.dir.RefOrbit = os.path.join(self.dir.SampleData, 'cs2')
        self.dir.RefCCOrbit = os.path.join(self.dir.SampleData, 'cs2cc')
        self.dir.SampleFigures = os.path.join(self.dir.workspace, 'testoutput')
        # Files
        self.file = EmptyObject()
        self.file.example_orbit = r'Location_file_5399_20110415T140524_20110415T154438.txt'
        self.file.example_des_orbit = r'Location_file_10520_20120402T110305_20120402T124219.txt'
        self.file.example_aem = os.path.join('aem',
                                             'HEM_PAM11_20110415T152401_20110415T164026.nc')


def example_calc_corner_coordinates():
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
    reforbit.limit_region(lat_limit=[82., 88.], lon_limit=[-120., -40.])

    # Read center coordinates (files provided from UCL for CryoVal-SI)
    # Input parameter: Full file name
    reforbit.from_file(os.path.join(info.dir.RefOrbit, info.file.example_orbit))

    # Save corner coordinates to file (in specified folder)
    reforbit.to_CCfile(folder=info.dir.refCCOrbit)


def example_resample_embird():
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
    cc_file = glob.glob(os.path.join(info.dir.RefCCOrbit, '*5399*'))[0]

    # Initialize reference orbit data object and read corner coordinate file
    reforbit = CS2RefOrbit()
    reforbit.from_file(cc_file, corner_coords=True)

    # Example: Read EM-Bird data
    #   Data structure of Cal/Val data can be arbitrary as long as
    #       * data is organized as vectors/array with identical shape
    #       * arrays exist for time, longitude, latitude, value (e.g. thickness, freeboard, ....)
    #       * time information is available as python datetime object (time zone UTC())
    aem_file = os.path.join(info.dir.SampleData, info.file.example_aem)
    aem = read_embird(aem_file)

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
    reforbit.limit_region(lat_limit=[82., 88.], lon_limit=[-120., -40.])
    reforbit.from_file(os.path.join(info.dir.RefOrbit, info.file.example_orbit))
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)


def test_descending_orbit():
    """ Test of SARin position corrections for descending orbits """
    info = PathInfo()
    reforbit = CS2RefOrbit()
    reforbit.limit_region(lat_limit=[82., 88.], lon_limit=[-120., -40.])
    reforbit.from_file(os.path.join(info.dir.RefOrbit, info.file.example_des_orbit))
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)


if __name__ == '__main__':
    example_resample_embird()
    test_descending_orbit()
    test_asscending_orbit()
# -*- coding: utf-8 -*-
import os
import glob


# Import Data Resampling Class from user library
from lib.pyoval import CS2RefOrbit
from lib.helpers import EmptyObject

class PathInfo():
    """ 
    Definition of pathes (UofA orbits!) and filenames for tests 
    """
    def __init__(self):
        # Path to reference orbit files
        self.dir = EmptyObject()
        self.dir.workspace = os.path.dirname(__file__)
        self.dir.SampleData = os.path.join(self.dir.workspace, 'testdata')
        self.dir.RefOrbit = os.path.join(self.dir.SampleData, 'cs2_orbit')
        self.dir.RefCCOrbit = os.path.join(self.dir.SampleData, 'cs2cc')
        self.dir.SampleFigures = os.path.join(self.dir.workspace, 'testoutput')
                                             

def test_uofa_orbits(orbit=5399, lon_limit=[82., 88.], lat_limit=[-120., -40.]):
    """ Test of UofA orbits with CS2RefOrbit debug map """
    # Get directory information
    # Note: RefOrbit set to cs2_orbit (UofA generated orbit positions)
    info = PathInfo()
    # Initialize reference orbit with given regions of interest
    reforbit = CS2RefOrbit()
    reforbit.limit_region(lat_limit=lat_limit, lon_limit=lon_limit)
    # This is new
    reforbit.set_sin_detection(False)
    # Locate orbit file for given orbit
    reforbit_file = glob.glob(os.path.join(info.dir.RefOrbit, '*'+str(orbit)+'_*'))
    reforbit.from_file(reforbit_file[0])
    # Display orbit position debug map
    reforbit.debug_map()
                                             
if __name__ == '__main__':
    # test with lincoln sea standard orbit
    test_uofa_orbits(orbit=5399, 
                     lat_limit=[82., 88.], 
                     lon_limit=[-120., -40.])                                             
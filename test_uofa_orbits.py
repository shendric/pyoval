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
        self.dir.workspace = '/Volumes/Data/OneDrive/CryoVal-SI/GoldenDays/'
        self.dir.SampleData = os.path.join(self.dir.workspace,'20120402')
        self.dir.RefOrbit = os.path.join(self.dir.workspace, '20120402')
        self.dir.RefCCOrbit = os.path.join(self.dir.workspace, '20120402')
        self.dir.SampleFigures = self.dir.SampleData
                                             

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
    print reforbit_file
    reforbit.from_file(reforbit_file[0])
    # Display orbit position debug map
    reforbit.debug_map()
    reforbit.to_CCfile(folder=info.dir.RefCCOrbit)
                                             
if __name__ == '__main__':
    # test with lincoln sea standard orbit

    # test with problematic (NPI) orbit                 
    test_uofa_orbits(orbit = 10520, 
                     lat_limit=[60., 88.], 
                     lon_limit=[-180.1, 180.1])
                                                
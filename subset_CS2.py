__version__ = 1.0
__author__ = "Justin Beckers, YorkU/UAlberta"
def subset_CS2(cs2file,ccfile,cc=0,version='V001'):
    '''
    subset_cs2
        Takes the CryoSat-2 L2 output from IDL, and the corner coordinate file
        produced by pyoval. The CC file is cut to the size of the CryoSat-2 
        file. Furthermore, if you give it both an AWI and an ESA CS2 file, it 
        makes sure that both have the same points and same extent (in case any
        differences do exist, we use the ESA CS2 file as the standard (since
        the CryoVal project is designed to examine the ESA CS2 data.)
        module: subset_cs2
            parameters:
                cs2file:    the .csv file with the CS2 L2 Baseline C data.
                            this can be an ESA file or AWI CS2 file, or both.
                ccfile:     the .dat file of corner coordinates
                cc:         flag should be set to 1 if you wish to output the 
                            CS2ORB corner coordinate file that is now sized to
                            the extent of the CS2 file.
                version:    processing version code string.
    Author: Justin Beckers
    Author Contact: beckers@ualberta.ca
    
    Verion: 1.0
    Version Notes:
        1.0: Initial release.
        
    Release Date: 2016/03/17
    
    Usage Notes: 
        Case 1: Using Main 
            From commandline/terminal, can simply run:
                >>>python subset_CS2.py
            This executes the program with the parameters set on 
            lines 327 - 349. Modify this section in the script file text to 
            meet your requirements
                if __name__=='__main__':
    '''


    #Import python modules.
    import os
    import numpy as np
    from shapely.geometry import Polygon,Point
    from rtree import index
    from shapely import speedups
    speedups.enable()
    
    #Add a quick converter  to convert longitude to 0/360 from -180/180.
    cv = lambda x: float(x)+360.0 if float(x) < 0.0 else float(x)
    
    #Load in the Cryosat-2 Corner Coordinates file. Convert the Longitude to 
    #0-360
    ccdata = np.genfromtxt(ccfile,converters={6:cv,8:cv,10:cv,12:cv,14:cv})
    
    #Load the Cryosat-2 L2 Data from AWI or ESA. Convert the Coordinates to
    #0-360
    if os.path.basename(cs2file).startswith('CS2AWI'): #If it is an AWI file
        cs2data=np.genfromtxt(cs2file,delimiter=',',skip_header=1,
        missing_values='nan',filling_values=np.nan,converters={1:cv})
    
    else: # It is assumed to be an ESA CS L2 file
        cs2data=np.genfromtxt(cs2file,delimiter=',',skip_header=1,
        missing_values='nan',filling_values=np.nan,converters={10:cv})

    #The AWI CS L2 data contains NaNs in the Lon,Lat, and Time fields that 
    #result from the L1B missing waveforms. These NaN values cause problems for
    #the polygon generation later on, and we don't really need them as they are
    #not in the ESA CS2 output from IDL (already cut out) so they are cut.
    foo=[] #A temporary variable
    #If of the longitude,latitude, or time are NaN, keep a list of the row.
    for j in np.arange(len(cs2data)):
        if (np.isnan(cs2data[j,1]) or np.isnan(cs2data[j,2]) or 
            np.isnan(cs2data[j,0])): 
            foo.append(j)
    cs2data=np.delete(cs2data,foo,axis=0) #Cut the bad rows.
   
    #Calculating a polygon over the Prime Meridian or the International Date 
    #Line can be troublesome (polygons sometimes wrap around the earth in the
    #wrong direction, depending on your coordinate system). Here we check
    #if we cross either the Prime Meridian or the I.D.L., and either convert to
    #0/360 coordinate system (already done), or back to -180/180
    converted=1 #Flag to indicate if we needed to convert the longitudes.
    if os.path.basename(cs2file).startswith('CS2AWI'): #If AWI CS2 file:
        if (np.any(cs2data[:,1] >=359.)==1 and np.any(cs2data[:,1] <=1.)==1 and
            np.any(np.logical_and(179.<cs2data[:,10],cs2data[:,10]<181.))==0):
            print "Crosses Prime Meridian but not IDL. Coordinates will be " 
            print "processed in -180/180 system."
            for i in np.arange(len(cs2data)):
                cs2data[i,1]=deconvert_lon(cs2data[i,1])
            for i in np.arange(len(ccdata)): 
                ccdata[i][6]=deconvert_lon(ccdata[i][6])
                ccdata[i][8]=deconvert_lon(ccdata[i][8])
                ccdata[i][10]=deconvert_lon(ccdata[i][10])
                ccdata[i][12]=deconvert_lon(ccdata[i][12])
                ccdata[i][14]=deconvert_lon(ccdata[i][14])
            converted=0
        elif (np.any(cs2data[:,1] >=359.)==0 and np.any(cs2data[:,1] <=1.)==0 
                and np.any(np.logical_and(
                179.<cs2data[:,10],cs2data[:,10]<181.))==1):
            print "Does not cross prime meridian but crosses IDL. Coordinates " 
            print "will be processed in 0/360 system."
            converted=1
        elif (np.any(cs2data[:,1] >=359.)==1 and np.any(cs2data[:,1] <=1.)==1 
                and np.any(np.logical_and(179.<cs2data[:,10],
                cs2data[:,10]<181.))==1):
            print "Crosses both the Prime Meridian and the IDL. Coordinates " 
            print "will be processed in 0/360 system."
            converted=1
        elif (np.any(cs2data[:,1] >=359.)==0 and np.any(cs2data[:,1] <=1.)==0 
                and np.any(np.logical_and(179.<cs2data[:,10],
                    cs2data[:,10]<181.))==0):
            print "Does not cross the IDL or the Prime Meridian. Coordinates "
            print "will be processed in 0/360 system."
            converted=1
    else: #If ESA CS2 file
        if (np.any(cs2data[:,10] >=359.) and np.any(cs2data[:,10] <=1.) and
                np.any(np.logical_and(179.<cs2data[:,10],
                cs2data[:,10]<181.))==0):
            print "Crosses Prime Meridian but not IDL. Coordinates will be " 
            print "processed in -180/180 system."
            for i in np.arange(len(cs2data)):
                cs2data[i,10]=deconvert_lon(cs2data[i,10]) 
            for i in np.arange(len(ccdata)):
                ccdata[i][6]=deconvert_lon(ccdata[i][6])
                ccdata[i][8]=deconvert_lon(ccdata[i][8])
                ccdata[i][10]=deconvert_lon(ccdata[i][10])
                ccdata[i][12]=deconvert_lon(ccdata[i][12])
                ccdata[i][14]=deconvert_lon(ccdata[i][14])
            converted=0
            pass
        elif (np.any(cs2data[:,10] >=359.)==0 and np.any(cs2data[:,10] <=1.)==0
                and np.any(np.logical_and(179.<cs2data[:,10],
                cs2data[:,10]<181.))==1):
            print "Does not cross prime meridian but crosses IDL. Coordinates "
            print "will be processed in 0/360 system."
            converted=1
        elif (np.any(cs2data[:,10] >=359.)==1 and np.any(cs2data[:,10] <=1.)==1
                and np.any(np.logical_and(179.<cs2data[:,10],
                cs2data[:,10]<181.))==1):
            print "Crosses both the Prime Meridian and the IDL. Coordinates "
            print "will be processed in 0/360 system."
            converted=1
        elif (np.any(cs2data[:,10] >=359.)==0 and np.any(cs2data[:,10] <=1.)==0 
                and np.any(np.logical_and(179.<cs2data[:,10],
                cs2data[:,10]<181.))==0):
            print "Does not cross the IDL or the Prime Meridian. Coordinates "
            print "will be processed in 0/360 system."
            converted=1
    
    #Setup some variables for later       
    n_records=len(ccdata)#Number of polygons
    idx=index.Index() #Setup a spatial index
    p=np.ndarray((n_records),dtype=object) #array to hold the polygons
    newpolys=[] #holds the index positions of the cc polygons
    newcs2=[] #holds the index positions of the cs2 data
    outpoly=[] #output polygon holder
    t=np.ndarray(0) # a temporary variable 
    nanline = np.empty(len(cs2data[0]))*np.nan
    
    #Fill the spatial index with the CS2 data
    if os.path.basename(cs2file).startswith('CS2AWI'): #If AWI CS2
        for i in np.arange(len(cs2data)):
            idx.insert(i,Point(cs2data[i,1],cs2data[i,2]).bounds)
    else: #If ESA CS2
        for i in np.arange(len(cs2data)):
            idx.insert(i,Point(cs2data[i,10],cs2data[i,9]).bounds)
    
    #Calculate the polygons
    n_polygon_points=5
    polygon_points = np.ndarray(shape=(n_polygon_points, 2), dtype=np.float64)        
    for i in np.arange(n_records):
        #self.cs2fp.lon_ur[i], self.cs2fp.lon_ul[i], self.cs2fp.lon_ll[i], self.cs2fp.lon_lr[i]]
        fp_x = [ccdata[i,8],ccdata[i,10],ccdata[i,14],ccdata[i,12],ccdata[i,8]]
        fp_y = [ccdata[i,9],ccdata[i,11],ccdata[i,15],ccdata[i,13],ccdata[i,9]]
        polygon_points[:,0]=fp_x[j] #Polygon X coordinates
        polygon_points[:,1]=fp_y[j] #Polygon Y coordinates
        p[i]=Polygon(polygon_points) #Builds the polygons

    #The work: go through each polygon, find the CS2 points that belong 
    #There should only be 1. #Intersection does not imply containment so need 
    #to test for intersection (grabs only the possible data), then test for 
    #containment.
    for i,poly in enumerate(p):
        for j in idx.intersection(poly.bounds): #Test for intersection
            if poly.area>0.01: #because polygon is crossing prime meridian but 
                #in the wrong direction. Let's transform the polygon 
                #coordinates back to -180/180 system, transform the points to 
                #-180/180 and test for containment of points.
                points=list(poly.exterior.coords)
                points_x,points_y=zip(*points)
                newx=[] #holds the retransformed points
                for i in np.arange(len(points_x)):
                    if points_x[i]>180.0:
                        newx.append(points_x[i]-360.0)
                    else:
                        newx.append(points_x[i])
                points_x=newx
                polygon_points[:,0]=points_x[:]
                polygon_points[:,1]=points_y[:]
                newpoly=Polygon(polygon_points)
                #Do the actual test for containment.
                if os.path.basename(cs2file).startswith('CS2AWI'):
                    if (newpoly.contains(Point(deconvert_lon(cs2data[j,1]),
                            cs2data[j,2]))):
                        newpolys.append(ccdata[i])
                        newcs2.append(cs2data[j])
                        outpoly.append(poly)
                        t=np.append(t,cs2data[j,0])
                    else:
                        newcs2.append(nanline)
                        newpolys.append(ccdata[i])
                        outpoly.append(poly)
                        pass
                else:
                    if (newpoly.contains(Point(deconvert_lon(cs2data[j,10]),
                            cs2data[j,9]))==1):
                        newpolys.append(ccdata[i])
                        newcs2.append(cs2data[j])
                        outpoly.append(poly)    
                    else:
                        pass
            else:
                if os.path.basename(cs2file).startswith('CS2AWI'):
                    if poly.contains(Point(cs2data[j,1],cs2data[j,2]))==1:
                        newpolys.append(ccdata[i])
                        newcs2.append(cs2data[j])
                        outpoly.append(poly)
                        t=np.append(t,cs2data[j,0])  
                    else:
                        pass #So if no points are found in the polygon, then we
                        #do not write out the polygon. This limits the output
                        #to the CS2 file extent. So if you wanted to keep each
                        #polygon, despite no CS2 data, you can do so here.
                else: #ESA CS2 file
                    if poly.contains(Point(cs2data[j,10],cs2data[j,9]))==1:
                        newpolys.append(ccdata[i])
                        newcs2.append(cs2data[j])
                        outpoly.append(poly)
                    else:
                        pass

    #Now we do some back conversion of lat/lon to -180/180 system
    if os.path.basename(cs2file).startswith('CS2AWI'):
        if converted==1:
            for i in np.arange(len(newcs2)):
                newcs2[i][1]=deconvert_lon(newcs2[i][1])
                newpolys[i][6]=deconvert_lon(newpolys[i][6])
                newpolys[i][8]=deconvert_lon(newpolys[i][8])
                newpolys[i][10]=deconvert_lon(newpolys[i][10])
                newpolys[i][12]=deconvert_lon(newpolys[i][12])
                newpolys[i][14]=deconvert_lon(newpolys[i][14])  
        else: #Was not converted out of -180/180 so we can leave it
            pass
    else:
        if converted ==1:
            for i in np.arange(len(newcs2)) :
                newcs2[i][10]=deconvert_lon(newcs2[i][10]) 
                newpolys[i][6]=deconvert_lon(newpolys[i][6])
                newpolys[i][8]=deconvert_lon(newpolys[i][8])
                newpolys[i][10]=deconvert_lon(newpolys[i][10])
                newpolys[i][12]=deconvert_lon(newpolys[i][12])
                newpolys[i][14]=deconvert_lon(newpolys[i][14]) 
        else:#Was not converted out of -180/180 so we can leave it
            pass
   
    #Print some information to the user. This is a useful check that the file
    #was correctly sized.
    print "Length of CCFile to start: ", len(ccdata)
    print "Length of CS2data: ",len(newcs2)," Length of CCs: ",len(newpolys)    
    print "Length of CCs == Length of CS2data: ",len(newcs2)==len(newpolys)
    print ""
    
    #Let's write the new corner coordinate file to CS2ORB_*.dat?
    if cc==1:#Yes let's write it
        if os.path.basename(cs2file).startswith('CS2AWI'):
            ccoutfile=os.path.join(os.path.dirname(ccfile),
                'CS2ORB_'+
                os.path.basename(ccfile).split('_')[2].zfill(6)+'_'+
                os.path.basename(cs2file).split('_')[2]+'_'+
                os.path.basename(cs2file).split('_')[3]+'_'+
                os.path.basename(cs2file).split('_')[4][:-4]+
                '_'+version+os.path.basename(ccfile)[-4:])    
        else:
            

            ccoutfile=os.path.join(os.path.dirname(ccfile),
                'CS2ORB_'+
                os.path.basename(ccfile).split('_')[2].zfill(6)+'_'+
                os.path.basename(cs2file).split('_')[2]+'_'+
                os.path.basename(cs2file).split('_')[3]+'_'+
                os.path.basename(cs2file).split('_')[4][:-4]+
                '_'+version+os.path.basename(ccfile)[-4:])
    with open(ccoutfile,'w') as ccout:
            np.savetxt(ccout,newpolys,delimiter=' ',fmt='%14.8f', newline='\n')        
    
    #Setup the output CS2 file names.
    if os.path.basename(cs2file).startswith('CS2AWI'):
        csoutfile=os.path.join(os.path.dirname(ccfile),
            os.path.basename(cs2file)[:-14]+'_'+version+
            os.path.basename(cs2file)[-4:])
    else:          
        csoutfile=os.path.join(os.path.dirname(ccfile),
            os.path.basename(cs2file).split('_')[0]+'_'+
            os.path.basename(cs2file).split('_')[1].zfill(6)+'_'+
            os.path.basename(cs2file).split('_')[2]+'_'+
            os.path.basename(cs2file).split('_')[3]+'_'+
            os.path.basename(cs2file).split('_')[4][:-4]+
            '_'+version+os.path.basename(cs2file)[-4:])
    #Write the cs2 output data
    #Open the file, read in the header line   
    with open(cs2file)as r:
        header=r.readline()
    #Now write the output file, first writing in the header line.
    with open(csoutfile,'w') as csout:
        csout.write(header)
        np.savetxt(csout,newcs2,delimiter=',',fmt='%14.8f', newline='\n')

def deconvert_lon(x):
    '''
    Converts longitude in 0-360 system to -180/180 system
    '''
    import numpy as np
    x=np.where(x>180.0,x-360.0,x)
    return x  
                                     
if __name__ == '__main__':
    '''
    If called from the command-line execute a recursive search of the current 
    working directory for all matches of corner coordinate and airborne data
    files and execute the code above for both shp and dat files for each file.
    '''
    import os
    import traceback
    
    ccfileList=[u'/Volumes/Data/Projects/CryoVal-SI_CCN/20140326_021016/REFERENCE/cs2_reforbit_021016_20140326T171524_20140326T173439_C001.dat']
    cs2fileList=["/Volumes/Data/Projects/CryoVal-SI_CCN/20140326_021016/CS2L2I_21016_20140326T172805_20140326T172902_C001.csv"]

    for i,cfile in enumerate(ccfileList):
        print "File: ",i,": ","GoldenDay: ",os.path.dirname(cfile).split("/")[3]
        print "Orbit File: ",os.path.basename(cfile)
        print "CS2 file: ",os.path.basename(cs2fileList[i])
        orb = os.path.basename(ccfileList[i]).split('_')[2]	
        print "Orbit: ",orb
        try:
            subset_CS2(cs2fileList[i],ccfileList[i],cc=1,version='A_V001')
        except:
            print traceback.format_exc()
            pass
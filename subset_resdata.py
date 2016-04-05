__version__= 1.1
__author__= "Justin Beckers, YorkU/UAlberta"

def subset_resdata(ccfile,subdatafile,version="A_V001"):
    '''Takes the resampled airborne data (cs2_resdata_ABCDEFn*) and cuts and
    extends it to the extent of the relevent CryoSat-2 file.The input ccfile
    is therefore the CS2ORB files produced by subset_cs2 code that ensures that
    the ESA and AWI files are of the size extent.
    
    Author: Justin Beckers
    Author Contact: beckers@ualberta.ca
    
    Verion: 1.1
    Version Notes:
        1.0: Initial release.
        1.1: Added in check for more than one resampled point per polygon.
             This has been done post CryoVal-SI analysis work.
    Release Date: 2016/03/14
    
    Inputs:
        ccfile: CryoSat-2 corner coordinate file (CS2ORB*, cs2_reforbit*)
        subdatafile: File to be subset to ccfile extent.
        Version: file version output code string.
    Outputs:
        AEMSIT*
        ALSLFB*
        ASRRFB*
        ASRSSH*
        ASRSSA*
        
    Usage Notes: 
        Case 1: Using Main 
            From commandline/terminal, can simply run:
                >>>python susbset_resdata.py
            This executes the program with the parameters set on 
            lines 205 - 213. Modify this section in the script file text to 
            meet your requirements
                if __name__=='__main__':
        Case 2: Import function:
        import subset_resdata
        subset_resdata(ccfile,subdatafile,version="versionstr")     
    '''

    
    #Load in the required python libraries
    import os
    import numpy as np
    import pandas as pd
    from shapely.geometry import Point,Polygon
    from shapely import speedups
    from rtree import index
    import re
    speedups.enable() #enable faster analysis of polygons and points with shapely
    
    #Convert from -180/180 to 0-360 longitude systems
    cv = lambda x: float(x)+360.0 if float(x) <0.0 else float(x)
    
    #Read the corner coordinate file.
    ccdata = np.genfromtxt(ccfile,converters={6:cv,8:cv,10:cv,12:cv,14:cv})
    
    #Quick message just to indicate that
    print "Assumes that airborne data has already been resampled."\
    " If not, you need to run pyoval resampling scripts first."
    
    #Extract the type of airborne data from the subdatafile filename.
    t = '_'+os.path.basename(subdatafile).split('_')[2]+'_'
    v = re.sub("[^A-Za-z]","",t)
    z='_'+v+'_'
    
    #Use the type of airborne data extracted on lines 47-49 above to setup the
    #output filenames and field names.
    if "asr" in v:
        if z=='_asr_':
            z='_fb_'
            x="ASRRFB_"
        elif z=='_asrssha_':
            z='_ssha_'
            x='ASRSSA_'
        elif z=='_asrssh_':
            z='_ssh_'
            x='ASRSSH_'       
        nam=['ASR'+z+'Year','ASR'+z+'Month','ASR'+z+'Day','ASR'+z+'Hour',
            'ASR'+z+'Minute','ASR'+z+'Seconds','ASR'+z+'lon','ASR'+z+'lat',
            'ASR'+z+'dT','ASR'+z+'nazg','ASR'+z+'Mean','ASR'+z+'Median',
            'ASR'+z+'stdez','ASR'+z+'Mode']
    elif "als" in v:
        if z=='_als_':
            z='_fb_'
            x="ALSLFB_"
        elif z=='_alsssha_':
            z='_ssha_'
            x='ALSSSA_'
        elif z=='_alsssh_':
            z='_ssh_'
            x='ALSSSH_' 
        nam=['ALS'+z+'Year','ALS'+z+'Month','ALS'+z+'Day','ALS'+z+'Hour',
            'ALS'+z+'Minute','ALS'+z+'Seconds','ALS'+z+'lon','ALS'+z+'lat',
            'ALS'+z+'dT','ALS'+z+'nazg','ALS'+z+'Mean','ALS'+z+'Median',
            'ALS'+z+'stdez','ALS'+z+'Mode']
    elif "aem" in v:
         nam=['AEM_Year','AEM_Month','AEM_Day','AEM_Hour',
            'AEM_Minute','AEM_Seconds','AEM_lon','AEM_lat',
            'AEM_dT','AEM_navg','AEM_Mean','AEM_Median',
            'AEM_stdev','AEM_Mode']
    else:
        print "Error input file has not been resampled, or is unsupported."
    
    
    #Read in the subdatafile with pandas. Convert the longitude field to 0-360.
    #Here is also where we could expand future output to include the
    #probability distribution by setting usecols=None.
    subdata=pd.read_table(subdatafile,sep='\s+',header=None,
                usecols=np.arange(14),converters={6:cv},index_col=False,
                na_values=['nan'],names=nam)
    
    #Setup the spatial index and relationship tree.    
    idx=index.Index()
    p=[] #empty list to hold the polygons
    n_records=len(ccdata) #Number of records in the corner coordinate file.
    #Empty array to hold results
    newsub=np.empty((n_records,len(subdata.iloc[0])))*np.nan 
    #A line of nans the size of the input airborne resampled data
    nanline=np.empty(len(subdata.iloc[0]))*np.nan
    n_polygon_points=5 #number of points in each polygon
    #Empty array to hold the polygon points
    polygon_points = np.ndarray(shape=(n_polygon_points, 2), dtype=np.float64)
    #For every coorner coordinate, generate the polygon.
    for i in np.arange(n_records):
        fp_x = [ccdata[i,8],ccdata[i,10],ccdata[i,14],ccdata[i,12],
            ccdata[i,8]]
        fp_y = [ccdata[i,9],ccdata[i,11],ccdata[i,15],ccdata[i,13],
            ccdata[i,9]]
        for j in np.arange(n_polygon_points):
            polygon_points[j,0]=fp_x[j]
            polygon_points[j,1]=fp_y[j]
        p.append(Polygon(polygon_points))
    #Generate a spatial index
    for j in np.arange(len(subdata[nam[6]])):
        idx.insert(j,Point(subdata.iloc[j][nam[6]],
            subdata.iloc[j][nam[7]]).bounds)
    '''
    There should only be one resampled datapoint for each footprint. 
    V1.1: Let's test that this is true, and let the user know. During CryoVal-
    SI several concatenated files (include ALS files split into 5ths during 
    resampling) were tested and none showed overlapping footprints. ASR files
    were also tested.
    
    Since this code is used to subset the resampled data to the extent of 
    other data, we want to keep a record for each polygon. Polygons with no
    resampled data should get a NaN value.
    '''
    newpts=[] #holds the rows where there are airborne points within a polygon
    out=[] #holds the polygon where airborne data contains points.
    nanlines=[] #holds the rows with no airborne data points in a polygon
    for i,poly in enumerate(p):
        n_pts=0
        for j in idx.intersection(poly.bounds):
            if poly.contains(Point(subdata.iloc[j][nam[6]],
                subdata.iloc[j][nam[7]]))==1:
                    n_pts+=1
                    newpts.append(j)
                    out.append(i)
                        
        #If no points in polygon, we will fill it with nans, so keep track
        if n_pts == 0:
            nanlines.append(i)
        elif n_pts >1:
            print "Hmm, found %i points in this polygon: %i " %n_pts,i
            
    newsub[out]=subdata.iloc[newpts] 
    newsub[nanlines]=nanline #not entirely sure this is needed.
    outdf=pd.DataFrame(data=newsub,columns=nam)
    outdf[nam[6]]=deconvert_lon(outdf[nam[6]]) #convert the longitude back.
    
    #Print some statistics so we know that everything is okay.
    nfinite = np.count_nonzero(np.isfinite(outdf[nam[10]]))
    print v+" file: ",os.path.basename(subdatafile)
    print "CC file: ",os.path.basename(ccfile)
    print "Number of input "+v+" data points: ",len(subdata)
    print "Number of footprints in Orbit file: ",len(ccdata)
    print "Number of output "+v+" Data Points: ", outdf.shape[0]
    print "Number of valid "+v+" footprints: ",nfinite
   
   #Create the output filename
    outfile=os.path.join(os.path.dirname(subdatafile),x+
        os.path.basename(ccfile).split('_')[1]+'_'+
        os.path.basename(ccfile).split('_')[2]+'_'+
        os.path.basename(ccfile).split('_')[3]+'_'+
        os.path.basename(ccfile).split('_')[4]+'_'+version+
        '.csv')
    print "Output file name: ",os.path.basename(outfile)
    
    #Write out the data.
    outdf.to_csv(path_or_buf=outfile,sep=',',na_rep='nan',header=True,
        index=False,mode='w',float_format='%14.8f',index_label='ID')  
            
def deconvert_lon(x):
    '''Converts from 0-360 to -180/180 longitude systems'''
    import numpy as np
    x=np.where(x>180.0,x-360.0,x)
    return x     

if __name__=='__main__':
    import time #Used to roughly time the run time of the code.
    
    ccfile=['/Volumes/Data/Projects/CryoVal-SI/GoldenDays_V001/20110415_005399/CS2ORB_005399_20110415T142828_20110415T143357_C001_A_V001.dat']
    alsfile=['/Volumes/Data/Projects/CryoVal-SI/GoldenDays_V001/20110415_005399/ASIRAS/cs2_resdata_asr_005399_20110415T142828_20110415T142828_concat.csv']
    a=time.clock()
    subset_resdata(ccfile[0],alsfile[0],version="A_V001")
    b=time.clock()
    print "Time to run the script (s): ",b-a   
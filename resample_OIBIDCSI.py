import traceback
import numpy as np
from rtree import index
import pandas as pd
__version__=1.0
__author__="Justin Beckers, YorkU/UAlberta"

'''
    
    Resamples the NASA Operation IceBridge IDCSI data product.
    Takes the IDCSI file and resamples it to the footprint corner coordinates
    file, providing statistics for each footprint.
    
    This code requires that the CS2ORB file (orbit file sized to Cryosat2 file)
    has already been generated. If not, the file will produce a file for the
    extent of any cc file (but then contains many polygons with no data and no
    CryoSat-2 data). However this functionality may be useful.
    The code performs similar operations to the pyoval CS2OrbitResData.resample
    function. However, as of version 1.0 this script does not produce a map or 
    other pdf summary. Furthermore, because the NASA OIB IDCSI product is a 40m
    average product, and this code is intended to resample the data to CryoSat
    footprints which are 300m along track (thus ~7-9 points per footprint), a
    statistical probability distribution is not calculated, thus we have no
    calculation of the modal value in each footprint. 
    Mode/PDF calculations may be implemented in a future release.
    
    
    Author: Justin Beckers
    Author Contact: beckers@ualberta.ca
    
    Verion: 1.0
    Version Notes:
        1.0: Initial release.
        
    Release Date: 2016/03/17
    
    Inputs:
        oibfile:  Operation IceBridge IDCSI data file.
        ccfile:   CryoSat-2 corner coordinate file (CS2ORB*, cs2_reforbit*)
        resample: Flag set to 1 if you wish to resample the OIB IDCSI file
        exp_dat:  Flag set to 1 if writing the resampled data to a file
        version:  File version output code string.
    Outputs:
        OIBSIX*.csv: where X is either 2,4 or Q for the NASA OIB IDCSI 
                     quicklook product
        
    Usage Notes: 
        Case 1: Using Main 
            From commandline/terminal, can simply run:
                >>>python process_OIB.py
            This executes the program with the parameters set on 
            lines 497 - 511. Modify this section in the script file text to 
            meet your requirements
                if __name__=='__main__':
        Case 2: Import functions:
        import process_OIB
        import deconvert_lon
        import resample_OIB
        import read_OIB
        import export_dat
             
    '''

def process_OIB(oibfile,ccfile, resample=1,exp_dat=1,version=''):
    '''
    The main processing script. Calls the other functions to read, resample,
    export the data.
    '''
    try:
        print "Processing OIB file: ", oibfile
        print "Processing Orbit file: ",ccfile
        oibdata, ccdata = read_OIB(oibfile,ccfile)
        print "Length of Orbit file: ", len(ccdata)
        print "Original Length of IceBridge File: ", len(oibdata)
        res=resample
        if res == 1:
            if "ql" in oibfile:
                print "OIB file is a quicklook file"
                outdata = resample_OIB(oibdata,ccdata,ql=1)
            else:
                print "OIB file is a final version file"
                outdata = resample_OIB(oibdata,ccdata,ql=0)

        if outdata.shape[0] == len(ccdata):
            print "OIB output file and orbit file are same length: SUCCESS!!"
        else:
            print "Strange!? The OIB output and orbit file are not the same length. Error!","\n"
        if exp_dat == 1:
            if "ql" in oibfile:
                    outfile = os.path.join(os.path.dirname(ccfile),'OIBSIQ_'+
                        os.path.basename(ccfile).split('_')[1]+'_'+
                        os.path.basename(ccfile).split('_')[2]+'_'+
                        os.path.basename(ccfile).split('_')[3]+'_'+
                        os.path.basename(ccfile).split('_')[4]+version+
                        '_V001.csv')
            else:
                outfile = os.path.join(os.path.dirname(ccfile),'OIBSI'+
                    os.path.basename(oibfile).split('_')[0][-1]+'_'+
                    os.path.basename(ccfile).split('_')[1]+'_'+
                    os.path.basename(ccfile).split('_')[2]+'_'+
                    os.path.basename(ccfile).split('_')[3]+'_'+
                    os.path.basename(ccfile).split('_')[4]+version+
                    '_V001.csv')
            export_dat(outdata,outfile)
            print "Output path: ",outfile,"\n" 
    except:
        print traceback.format_exc()
        pass    
def export_dat(outdata,filename):
    '''Exports the resampled data to a csv file'''
    
    a=list(outdata.columns.values)
   
    del a[-1] #needed because we were hanging onto the polygon objects.
    b=[x for x in a if not "empty" in x] #needed to remove the polygon objects
    #Try Exporting the data:
    try:
        outdata.to_csv(path_or_buf=filename,sep=',',columns = b,
            na_rep='nan',header=True,index=False,mode='w',index_label='ID')
        Error = "Success"
    #Catch an errors and print them
    except:
        Error = "Error"
        print traceback.format_exc()
        pass
    
    return Error
#
def deconvert_lon(x):
    '''Converts 0-360 longitudes to -180/180 longitudes'''
    x=np.where(x>180.0,x-360.0,x)
    return x  

def read_OIB(oibfile,ccfile):
    '''Read in an oib IDCSI file'''
    import numpy as np
    #These correspond to the NASA OIB IDCSI columns of interest.
    cols =np.array([0,1,2,3,4,5,6,7,8,15,16,24,26,29,39])

    #Had great difficulty reading in the NASA OIB IDCSI NAN = -99999, -999 
    #as NaNs.
    oibdata = pd.read_table(oibfile,sep=',',header=0,
    usecols=cols,na_values=['-99999','-99999.','-99999.0','-99999.00','-999',
        '-999.','-99999.000','-99999.0000','-99999.00000','-99999.000000',
        '-99999.0000000','-99999.00000000','-99999.000000000',
        str(np.int64(-99999.000000)),'-99999.0000000000','nan',
        str(np.int64(-99999)),' -99999',np.int64(-99999),-99999.,-99999,
        -99999.0,-99999.00,-99999.000,-99999.0000,-99999.00000,-99999.000000,
        -99999.0000000,-99999.00000000,-99999.000000000,int(-99999),
        np.int64(-99999),np.float64(-99999),np.int64(-99999.000000)]
        )
    #Had an issue with a few of the ICDSI2 files or QL files where the latitude
    #column was not being read correctly, so we had to read in the data another
    #way.
    if oibdata['lat'].max() >90.0:
        print "Error reading data using Pandas. Lat column not being read"
        print "(affects 20110317 and 20110415 data files)."
        print "Reading data from numpy.genfromtxt instead."
        oibdata=np.empty(1)
        #read using numpy
        bar=np.genfromtxt(oibfile,delimiter=',',names=True,usecols=cols,
        missing_values=['-99999','-99999.','-99999.0','-99999.00','-99999.000',
            '-999.','-999','-99999.0000','-99999.00000','-99999.000000',
            '-99999.0000000','-99999.00000000','-99999.000000000',
            '-99999.0000000000','nan',' -99999',str(np.int64(-99999.000000)),
            str(np.int64(-99999)),np.int64(-99999),-99999.,-99999,-99999.0,
            -99999.00,-99999.000,-99999.0000,-99999.00000,-99999.000000,
            -99999.0000000,-99999.00000000,-99999.000000000,int(-99999),
            np.float64(-99999),np.int64(-99999.000000)],
            filling_values=np.nan,usemask=False)
        #Convert loaded numpy data array to Pandas dataframe
        oibdata=pd.DataFrame(data=bar,columns = list(bar.dtype.names))
        oibdata.replace(['-99999','-99999.','-99999.0','-99999.00','-99999.000',
            '-99999.0000','-99999.00000','-99999.000000','-99999.0000000',
            '-999.','-999','-99999.00000000','-99999.000000000',
            '-99999.0000000000','nan',str(np.int64(-99999.000000)),
            str(np.int64(-99999)),np.int64(-99999),' -99999',' -99999',
            -99999.,-99999,-99999.0,-99999.00,-99999.000,-99999.0000,
            -99999.00000,-99999.000000,-99999.0000000,-99999.00000000,
            -99999.000000000,int(-99999),np.float64(-99999),
            np.int64(-99999.000000)],value=np.nan,inplace=True)
 
    #Read in the cc file.
    cv = lambda x: float(x)+360.0 if float(x) < 0.0 else float(x)
    ccdata = np.genfromtxt(ccfile,converters={6:cv,8:cv,10:cv,12:cv,14:cv})

    return oibdata,ccdata

def resample_OIB(oibdata,ccdata,ql=1):
    '''
    Does the resampling of the OIB IDCSI files.
    Uses a spatial index to search points fast. Does not use prepared polygons
    '''
    from shapely.geometry import Point, Polygon
    from shapely import speedups 
    import datetime
    import numpy as np
    import pytz
    #Setup which columns we want to calculate which statistics for.
    meancols = [0,1,2,3,4,5,6,7,8,11,12,13]
    stdevcols = [2,3,4,5,6,7,8,11,12]
    mediancols = [2,3,4,5,6,7,8,11,12]
    maxmincols = [2,3,4,5,6,7,8,11,12,13] 
    speedups.enable()
    #Check for crossings of IDL/prime meridian
    if ((oibdata['lon'] >=359.).any(skipna=True)==1 and 
        (oibdata['lon']<1.).any(skipna=True)==1 and 
        np.logical_and((179.<oibdata['lon']).any(skipna=True),
        (oibdata['lon']<181.).any(skipna=True)==0)):
        print "Crosses Prime Meridian but not IDL. Coordinates will be" 
        print "processed in -180/180 system."
        oibdata['lon']=deconvert_lon(oibdata['lon'])
        for i in np.arange(len(ccdata)): 
            ccdata[i][6]=deconvert_lon(ccdata[i][6])
            ccdata[i][8]=deconvert_lon(ccdata[i][8])
            ccdata[i][10]=deconvert_lon(ccdata[i][10])
            ccdata[i][12]=deconvert_lon(ccdata[i][12])
            ccdata[i][14]=deconvert_lon(ccdata[i][14])
            converted=0
    elif ((oibdata['lon'] >=359.).any(skipna=True)==0 and
        (oibdata['lon']<1.).any(skipna=True)==0 and 
        np.logical_and((179.<oibdata['lon']).any(skipna=True),
        (oibdata['lon']<181.).any(skipna=True)==1)):
        print "Does not cross prime meridian but crosses IDL. Coordinates will" 
        print "be processed in 0/360 system."
        converted=1
    elif ((oibdata['lon'] >=359.).any(skipna=True)==1 and
        (oibdata['lon']<1.).any(skipna=True)==1 and 
        np.logical_and((179.<oibdata['lon']).any(skipna=True),
        (oibdata['lon']<181.).any(skipna=True)==1)):
        print "Crosses both the Prime Meridian and the IDL. Coordinates will"
        print "be processed in -180/180 system."
        oibdata['lon']=deconvert_lon(oibdata['lon'])
        for i in np.arange(len(ccdata)): 
            ccdata[i][6]=deconvert_lon(ccdata[i][6])
            ccdata[i][8]=deconvert_lon(ccdata[i][8])
            ccdata[i][10]=deconvert_lon(ccdata[i][10])
            ccdata[i][12]=deconvert_lon(ccdata[i][12])
            ccdata[i][14]=deconvert_lon(ccdata[i][14])
            converted=0
    elif ((oibdata['lon'] >=359.).any(skipna=True)==0 and 
        (oibdata['lon']<1.).any(skipna=True)==0 and 
        np.logical_and((179.<oibdata['lon']).any(skipna=True),
        (oibdata['lon']<181.).any(skipna=True)==0)):
        print "Does not cross the IDL or the Prime Meridian. Coordinates will" 
        print "be processed in 0/360 system."
        converted=1

    #Try calculating polygons, filling the index.
    try:
        p=[]
        idx=index.Index() #Setup a spatial indec
        n_records=len(ccdata)
        n_polygon_points=5
        polygon_points =np.ndarray(shape=(n_polygon_points,2),dtype=np.float64)
        p=np.ndarray(shape=(n_records),dtype=object)
        #Calculate the polygons for the CryoSat-2 (cc) footprints
        for i in np.arange(n_records):
            fp_x = [ccdata[i,8],ccdata[i,10],ccdata[i,14],ccdata[i,12],
                ccdata[i,8]]
            fp_y = [ccdata[i,9],ccdata[i,11],ccdata[i,15],ccdata[i,13],
                ccdata[i,9]]
            n_polygon_points=len(fp_x)
            polygon_points[:,0]=fp_x[:]
            polygon_points[:,1]=fp_y[:]
            p[i](Polygon(polygon_points))
        #Fill the spatial index with the icebridge data points
        for j in np.arange(len(oibdata['lon'])):
            idx.insert(j,Point(oibdata.iloc[j]['lon'],
                oibdata.iloc[j]['lat']).bounds)
        #Setup our variables to hold the resampled data. Preallocation is fast.
        newpolys = []
        polys = []
        navg=np.empty((n_records,1))*np.nan
        foobar = np.empty((n_records))*np.nan
        myiflag = np.empty((n_records,1))*np.nan #mode of my ice flag (0 or 1)
        u=np.empty((n_records,len(meancols)))*np.nan #mean
        s=np.empty((n_records,len(stdevcols)))*np.nan #stdev
        m=np.empty((n_records,len(mediancols)))*np.nan #median
        h=np.empty((n_records,len(maxmincols)))*np.nan #max
        l=np.empty((n_records,len(maxmincols)))*np.nan #min
        out = np.empty((n_records,(len(stdevcols)+len(mediancols)+(2*len(maxmincols))+len(meancols)+1)))
        #Let's try the resampling.
        try:
            for i,poly in enumerate(p):   
                newpts=[]
                n_pts = 0
                #Use spatial index to find points that intersect the polygon.
                #Just because a point intersects the polygon, does not mean 
                #that it is inside the polygon, so we test that further on.
                #j is an iterative point in the list of points (in the spatial)
                #index that does intersect with the polygon bounds.
                for j in idx.intersection(poly.bounds): 
                    
                    #We discovered an issue with the polygons near the IDL/PM
                    #where they wrapped the wrong way around the Earth. We test
                    #for them here and fix them. 
                    if poly.area>0.1: 
                    #because polygon is crossing prime meridian area is large
                    #Let's transform the polygon coordinates back to -180/180 
                    #Transform the points to -180/180 and test
                        print poly 
                        #left in so that we notice when these cases happen.
                        points=list(poly.exterior.coords)
                        points_x,points_y=zip(*points)
                        newx=[]
                        for i in np.arange(len(points_x)):
                            if points_x[i]>180.0:
                                newx.append(points_x[i]-360.0)
                            else:
                                newx.append(points_x[i])
                        points_x=newx
                        n_polygon_points=len(points_x)
                        polygon_points = np.ndarray(shape=(n_polygon_points, 2), dtype=np.float64)
                        for k in np.arange(n_polygon_points):
                            polygon_points[k,0]=points_x[k]
                            polygon_points[k,1]=points_y[k]
                        poly=Polygon(polygon_points)
                        #Test that the point is in the polygon
                        if poly.contains(
                                Point(deconvert_lon(oibdata.iloc[j]['lon']),
                                oibdata.iloc[j]['lat'])):
                            newpts.append(j)
                            n_pts+=1
                    #If the area is okay, let's move on and check if the point
                    #is in the polygon
                    else:    
                        if poly.contains(Point(oibdata.iloc[j]['lon'],
                            oibdata.iloc[j]['lat'])):
                            newpts.append(j)
                            n_pts+=1
                #let's append the polygon number (is keeping all polygons).
                #Relies on having already calculated orbit extent for CS2 file.
                newpolys.append(i) 
                #Based on the number of points, calculate the statistics
                #If no points, values are nan or 0 (for navg)
                if n_pts == 0:
                    u[i] = np.empty(len(meancols))*np.nan
                    s[i] = np.empty(len(stdevcols))*np.nan
                    m[i] = np.empty(len(mediancols))*np.nan
                    h[i] = np.empty(len(maxmincols))*np.nan
                    l[i] =np.empty(len(maxmincols))*np.nan
                    myiflag[i] = np.empty(1)*np.nan
                    #foobar will hold the median point of the points in a poly
                    foobar[i]= 0.0 #
                    navg[i]=0
                #If 1 point, calculate the statistics, but std will be nan,
                #others will just be the value of the point
                elif n_pts == 1: 
                    u[i] = oibdata.iloc[newpts,meancols].mean()
                    s[i] = oibdata.iloc[newpts,stdevcols].std()
                    m[i] = oibdata.iloc[newpts,mediancols].median()
                    h[i] = oibdata.iloc[newpts,maxmincols].max()
                    l[i] = oibdata.iloc[newpts,maxmincols].min()
                    myiflag[i] = oibdata.iloc[newpts,14].mean()
                    navg[i] = n_pts
                    foobar[i]=np.floor(np.median(newpts))
                #If more than one point (what we expect), calculate the 
                #statistics. Note that the mode (most common value) is 
                #calculated for all fields, but only the value for the 
                #MY_ICE_FLAG is kept.
                elif n_pts > 1:
                    foo = oibdata.iloc[newpts].mode()
                    u[i] = oibdata.iloc[newpts,meancols].mean()
                    s[i] = oibdata.iloc[newpts,stdevcols].std()
                    m[i] = oibdata.iloc[newpts,mediancols].median()
                    h[i] = oibdata.iloc[newpts,maxmincols].max()
                    l[i] = oibdata.iloc[newpts,maxmincols].min()
                    myiflag[i]=foo['my_ice_flag'].iloc[0]
                    navg[i] = n_pts
                    foobar[i]=np.floor(np.median(newpts))
        except:
            print traceback.format_exc()
        #Newpolys sould be unique anyways, but doublecheck and return only the
        #unique polygon records (no doubles)  
        polys=ccdata[np.unique(newpolys),:]
        #Concatenate the variables holding the resampled data.
        out=np.concatenate((u,s,m,h,l,myiflag),axis=1)
        #Out should have lenght of newpolys so again check for uniqueness
        out=out[np.unique(newpolys),:]
        #Same with navg.
        navg=navg[np.unique(newpolys),:]
        #Lets rewrite the header here to separate out IDCSI and IDCSI quickooks
        if ql==1:
            outhead=['OIBql_lat_mean','OIBql_lon_mean','OIBql_thickness_mean','OIBql_thickness_unc_mean','OIBql_mean_fb_mean',
            'OIBql_ATM_fb_mean','OIBql_fb_unc_mean','OIBql_snow_depth_mean','OIBql_snow_depth_unc_mean',
            'OIBql_ssh_mean','OIBql_ssh_sd_mean','OIBql_ssh_tp_dist_mean',
            'OIBql_thickness_stdev','OIBql_thickness_unc_stdev','OIBql_mean_fb_stdev',
            'OIBql_ATM_fb_stdev','OIBql_fb_unc_stdev','OIBql_snow_depth_stdev','OIBql_snow_depth_unc_stdev',
            'OIBql_ssh_stdev','OIBql_ssh_sd_stdev', 
            'OIBql_thickness_median','OIBql_thickness_unc_median','OIBql_mean_fb_median','OIBql_ATM_fb_median',
            'OIBql_fb_unc_median','OIBql_snow_depth_median','OIBql_snow_depth_unc_median',
            'OIBql_ssh_median','OIBql_ssh_sd_median',
            'OIBql_thickness_max','OIBql_thickness_unc_max','OIBql_mean_fb_max','OIBql_ATM_fb_max','OIBql_fb_unc_max',
            'OIBql_snow_depth_max','OIBql_snow_depth_unc_max','OIBql_ssh_max','OIBql_ssh_sd_max','OIBql_ssh_tp_dist_max',
            'OIBql_thickness_min','OIBql_thickness_unc_min',
            'OIBql_mean_fb_min','OIBql_ATM_fb_min','OIBql_fb_unc_min','OIBql_snow_depth_min','OIBql_snow_depth_unc_min',
            'OIBql_ssh_min','OIBql_ssh_sd_min','OIBql_ssh_tp_dist_min','OIBql_MYIflag_mode']
        else:
             outhead=['OIB_lat_mean','OIB_lon_mean','OIB_thickness_mean','OIB_thickness_unc_mean','OIB_mean_fb_mean',
            'OIB_ATM_fb_mean','OIB_fb_unc_mean','OIB_snow_depth_mean','OIB_snow_depth_unc_mean',
            'OIB_ssh_mean','OIB_ssh_sd_mean','OIB_ssh_tp_dist_mean',
            'OIB_thickness_stdev','OIB_thickness_unc_stdev','OIB_mean_fb_stdev',
            'OIB_ATM_fb_stdev','OIB_fb_unc_stdev','OIB_snow_depth_stdev','OIB_snow_depth_unc_stdev',
            'OIB_ssh_stdev','OIB_ssh_sd_stdev',
            'OIB_thickness_median','OIB_thickness_unc_median','OIB_mean_fb_median','OIB_ATM_fb_median',
            'OIB_fb_unc_median','OIB_snow_depth_median','OIB_snow_depth_unc_median',
            'OIB_ssh_median','OIB_ssh_sd_median',
            'OIB_thickness_max','OIB_thickness_unc_max','OIB_mean_fb_max','OIB_ATM_fb_max','OIB_fb_unc_max',
            'OIB_snow_depth_max','OIB_snow_depth_unc_max','OIB_ssh_max','OIB_ssh_sd_max','OIB_ssh_tp_dist_max',
            'OIB_thickness_min','OIB_thickness_unc_min',
            'OIB_mean_fb_min','OIB_ATM_fb_min','OIB_fb_unc_min','OIB_snow_depth_min','OIB_snow_depth_unc_min',
            'OIB_ssh_min','OIB_ssh_sd_min','OIB_ssh_tp_dist_min','OIB_MYIflag_mode']

        #Let's create some of the other variables, like the timestamp
        oib_time = np.ndarray(shape=(len(oibdata)),dtype=object)
        secs = np.empty(len(oibdata))*np.nan
        msecs = np.empty(len(oibdata))*np.nan
        secs = np.floor(oibdata['elapsed'].values.tolist())
        msecs = np.floor(1e6 * (oibdata['elapsed'].values.tolist() - np.floor(secs)))
        date = oibdata['date'].values.tolist()
        date=map(int,date)
        date=map(str,date)
        #Lets calculate the OIB timestamp.
        for i in np.arange(len(date)):
            oib_time[i] = datetime.datetime.strptime(date[i],"%Y%m%d").replace(
                tzinfo=pytz.utc) +\
                datetime.timedelta(seconds=int(secs[i]), microseconds=int(
                msecs[i]))
        foobar=foobar.astype(int)       
        #Get the for the middle point of all OIB points in each footprint
        oib_time=np.where(foobar==0,0,oib_time[foobar])

        #Let's calculate the CS2 timestamp from the CC data.
        cc_time=np.ndarray(shape=(len(polys)),dtype=object)
        for i in np.arange(len(cc_time)):
            cc_time[i] = datetime.datetime(int(ccdata[i,0]), int(ccdata[i,1]),int(ccdata[i,2]),
            int(ccdata[i,3]), int(ccdata[i,4]), int(ccdata[i,5]),int(1e6*(ccdata[i,5] - np.floor(ccdata[i,5]))),
            tzinfo=pytz.utc)
        #Let's calculate the difference between the CS2 time and OIB.
        dt_time=np.ndarray(shape=(len(cc_time)),dtype=object)
        for i in np.arange(len(cc_time)):
            if oib_time[i]==0:
                dt_time[i]=np.nan
            else:
                dt_time[i]=(cc_time[i]-oib_time[i]).total_seconds()
        
        #Check for uniqueness in the shapely polygon objects.
        g = np.unique(newpolys)
        c = [p[x] for x in g]
        #Setup the output dataframe
        outdf=pd.DataFrame(data=out,columns=outhead)

        #Add in the delta time, footprint latitude,longitude, and npts 
        if ql == 1:
            outdf['OIB_dt_time']=dt_time
            outdf['OIBql_fp_lat']=polys[:,7]
            outdf['OIBql_fp_lon']=polys[:,6]
            outdf['OIBql_n_pts']=navg[:,0]
            if converted == 1:
                outdf['OIBql_fp_lon']=deconvert_lon(outdf['OIBql_fp_lon'])
                outdf['OIBql_lon_mean']=deconvert_lon(outdf['OIBql_lon_mean'])
            nfinite = np.count_nonzero(np.isfinite(outdf['OIBql_mean_fb_mean']))
        else:
            outdf['OIB_dt_time']=dt_time
            outdf['OIB_fp_lat']=polys[:,7]
            outdf['OIB_fp_lon']=polys[:,6]
            outdf['OIB_n_pts']=navg[:,0]
            if converted == 1:
                outdf['OIB_fp_lon']=deconvert_lon(outdf['OIB_fp_lon'])
                outdf['OIB_lon_mean']=deconvert_lon(outdf['OIB_lon_mean'])
            nfinite = np.count_nonzero(np.isfinite(outdf['OIB_mean_fb_mean']))
        #Let's add in the polygon geometry objects
        
        outdf['OIB_geometry']=c
        #Print out a bit of info about the resampling result and return data.
        print "Number of Resampled IceBridge Data Points: ", outdf.shape[0]
        print "Number of finite freeboard values: ",nfinite
        return outdf
        
    except ValueError:
        pass
        
if __name__=='__main__':
    import os
    import time
    try:
        ccfileList=[u'/Volumes/Data/OneDrive/CryoVal-SI/GoldenDays_V001/20120315_010262/CS2ORB_010262_20120315T165011_20120315T165246_C001_V001.dat']
        oibfileList=[u'/Volumes/Data/OneDrive/CryoVal-SI/GoldenDays_V001/20120315_010262/IDCSI4_20120315.txt']
        a=time.clock()
        for i,ofile in enumerate(oibfileList):
            try:
                process_OIB(oibfileList[i],ccfileList[i], resample=1,exp_dat=1,version='')
            except:
                print traceback.format_exc()
                pass
        b=time.clock()
        print "Time to execute the script (s): ",b-a
    except:
        print traceback.format_exc()


__version__ = 1.0
__author__ = "Justin Beckers, YorkU/UAlberta"
def concatenate_output(indir,completeflag=1,version="_A_V001"):
    '''
    Takes a directory containing the individual output, resampled, subset files
    and joins them together in a single file.
    
    Author: Justin F. Beckers,
    Release Date: 2016/03/18
    
    Verions: 1.0
    Version Notes:
        1.0: Initial release
    
    Inputs:
        indir:        directory to search for input files. non-recursive.
                      AEM*,ASR*,ALS*,OIB*,CS2*. Explicitly does not include
                      CVLDAT* (the output of this script) and CS2ORB* (corner
                      coordinate file sized to CS2 files)
        completeflag: Set to 0 if producing concatenated output at some inter-
                      mediate stage where not all the data is available for a 
                      particular golden day yet (just incase). Flag is default
                      set to 1, which indicates you have all necessary files.
                      Only affects the filename (if 0, _incomplete is appended)
        version:      Processing version output code string. Used to only get
                      the same version of the data.
        
    Outputs: concatenated file contain all found files in a directory.
             Format of the files is csv.
    
    Usage:
        1).Can import function and run with string of indir as 
            import concatenat_output
            indir = 'somedirectory'
            concatenate_output(indir)
        2). Can run as main.
                Modify code on lines: 308-319 to your usage needs.
                Currently, input directory is set to current working directory.
                Furthermore we have implemented a directory walk to recursively
                walk through all subdirectories and look for the ESA CS2 file.
                If found, then we use that directory and run the concatente 
                script before moving to the next directory. This lets one 
                process an entire directory containing multiple subdirectories
                for each golden day in one go.
                
                python concatenate_output.py
    '''
    
    #Import Python Modules.
    import os
    import glob
    import pandas as pd
    import numpy as np  
     
    #Find relevent files in the input directory    
    cs2file = glob.glob(indir+'*CS2L2I*'+version+'*')
    AWIfile = glob.glob(indir+'*CS2AWI*'+version+'*')
    resdata = set(glob.glob(indir+'*'+version+'*')) - \
                set(glob.glob(indir+'*CS2*')) - \
                set(glob.glob(indir+'*CS2ORB*'))-\
                set(glob.glob(indir+'*CVLDAT*'))
    #Convert the sets back to a list            
    resdata=list(resdata)
    resdata.sort() #Sort the list into alphabetical order (AEM,ALS,ASR,OIB)
    
    #If a golden day does not include all data types, then we need to add in
    #columns of nans with the correct column names.
    addcols=0
    addhead=[]
    
    #If file is not file in the list of files, count the columns we should have
    #and create a header string for those columns.
    if any("CS2AWI" in s for s in AWIfile)==False: #No AWI CS2 file
        addcols+=25
        addhead+=['AWICS2_time','AWICS2_lon','AWICS2_lat','AWICS2_fb_stat',
        'AWICS2_fb_syst','AWICS2_fb_unc','AWICS2_fb','AWICS2_sit_stat',
        'AWICS2_sit_syst','AWICS2_sit_unc','AWICS2_sit_fb_syst',
        'AWICS2_sit_snow_syst','AWICS2_sit_snow_syst_iav',
        'AWICS2_sit_rho_i_syst','AWICS2_sit_rho_s_syst','AWICS2_sit',
        'AWICS2_ssha_unc','AWICS2_ssh','AWICS2_ssha','AWICS2_sd',
        'AWICS2_rho_i','AWICS2_rho_s','AWICS2_ice_type',
        'AWICS2_surface_type','AWICS2_ice_conc']
    if any("CS2L2I" in s for s in cs2file)==False: #No ESA CS2 file
        addcols+=49
        addhead+=['ID','CryoSatDateTime','Year','Month','Day','Hour','Minute',
        'Second','Orbit','Latitude','Longitude','Freeboard','FreeboardFlag',
        'Height_1','Peakiness','Sigma0_1','SHA','Interp_SHA',
        'Interp_OceanHt_Error','Interp_OceanHt_CNT_FWD',
        'Interp_OceanHt_CNT_BKWD','N_Avg','IceConc','SnowDepth','SnowDensity',
        'SAR_DISCRIM','MSS_Mod','Geoid_Mod','ODLE_Mod','DryT_C','WetT_C',
        'IB_C','DAC_C','IONO_GIM','IONO_Mod','H_OT_C','H_LPEOT_C','H_OLT_C',
        'H_SET_C','H_GPT_C','Beam_Center','Beam_Amp','Beam_Kurt','Beam_SD',
        'Beam_SDAng','Beam_Skew','Beam_Center_ang','BadDataFlag',
        'Corr_Error_F']
    if any("AEMSIT" in s for s in resdata)==False: #No AEM file
        addcols+=14
        addhead+=['AEM_Year','AEM_Month','AEM_Day','AEM_Hour','AEM_Minute',
        'AEM_Seconds','AEM_lon','AEM_lat','AEM_dT','AEM_navg','AEM_Mean',
        'AEM_Median','AEM_stdev','AEM_Mode']
    if any("ALS_nodes" in s for s in resdata)==False: #No ALS nodes file
        addcols+=1
        addhead+=['ALS_leads']
    if any("ASR_nodes" in s for s in resdata)==False: #No ASR nodes file
        addcols+=1
        addhead+=['ASR_leads']
    if any("ASRRFB" in s for s in resdata)==False: #No ASR FB file
        addcols+=14
        addhead+=['ASR_fb_Year','ASR_fb_Month','ASR_fb_Day','ASR_fb_Hour',
        'ASR_fb_Minute','ASR_fb_Seconds','ASR_fb_lon','ASR_fb_lat','ASR_fb_dT'
        ,'ASR_fb_navg','ASR_fb_Mean','ASR_fb_Median','ASR_fb_stdev',
        'ASR_fb_Mode']
    if any("ASRSSH" in s for s in resdata)==False: #No ASR SSH file
        addcols+=14
        addhead+=['ASR_ssh_Year','ASR_ssh_Month','ASR_ssh_Day','ASR_ssh_Hour',
        'ASR_ssh_Minute','ASR_ssh_Seconds','ASR_ssh_lon','ASR_ssh_lat',
        'ASR_ssh_dT','ASR_ssh_navg','ASR_ssh_Mean','ASR_ssh_Median',
        'ASR_ssh_stdev','ASR_ssh_Mode']
    if any("ASRSSA" in s for s in resdata)==False: #No ASR SSA file
        addcols+=14
        addhead+=['ASR_ssha_Year','ASR_ssha_Month','ASR_ssha_Day',
        'ASR_ssha_Hour','ASR_ssha_Minute','ASR_ssha_Seconds','ASR_ssha_lon',
        'ASR_ssha_lat','ASR_ssha_dT','ASR_ssha_navg','ASR_ssha_Mean',
        'ASR_ssha_Median','ASR_ssha_stdev','ASR_ssha_Mode']
    if any("ALSLFB" in s for s in resdata)==False: #No ALS FB file
        addcols+=14
        addhead+=['ALS_fb_Year','ALS_fb_Month','ALS_fb_Day','ALS_fb_Hour',
        'ALS_fb_Minute','ALS_fb_Seconds','ALS_fb_lon','ALS_fb_lat','ALS_fb_dT',
        'ALS_fb_navg','ALS_fb_Mean','ALS_fb_Median','ALS_fb_stdev',
        'ALS_fb_Mode']  
    if any("OIBSI" in s for s in resdata)==False: #No OIB file.
        #Note that here we also rename the oib headers for the OIB quicklook
        #to match the OIB IDCSI product.
        addcols+=55
        addhead+=['OIB_lat_mean','OIB_lon_mean','OIB_thickness_mean',
        'OIB_thickness_unc_mean','OIB_mean_fb_mean','OIB_ATM_fb_mean',
        'OIB_fb_unc_mean','OIB_snow_depth_mean','OIB_snow_depth_unc_mean',
        'OIB_ssh_mean','OIB_ssh_sd_mean','OIB_ssh_tp_dist_mean',
        'OIB_thickness_stdev','OIB_thickness_unc_stdev','OIB_mean_fb_stdev',
        'OIB_ATM_fb_stdev','OIB_fb_unc_stdev','OIB_snow_depth_stdev',
        'OIB_snow_depth_unc_stdev','OIB_ssh_stdev','OIB_ssh_sd_stdev',
        'OIB_thickness_median','OIB_thickness_unc_median','OIB_mean_fb_median',
        'OIB_ATM_fb_median','OIB_fb_unc_median','OIB_snow_depth_median',
        'OIB_snow_depth_unc_median','OIB_ssh_median','OIB_ssh_sd_median',
        'OIB_thickness_max','OIB_thickness_unc_max','OIB_mean_fb_max',
        'OIB_ATM_fb_max','OIB_fb_unc_max','OIB_snow_depth_max',
        'OIB_snow_depth_unc_max','OIB_ssh_max','OIB_ssh_sd_max',
        'OIB_ssh_tp_dist_max','OIB_thickness_min','OIB_thickness_unc_min',
        'OIB_mean_fb_min','OIB_ATM_fb_min','OIB_fb_unc_min',
        'OIB_snow_depth_min','OIB_snow_depth_unc_min','OIB_ssh_min',
        'OIB_ssh_sd_min','OIB_ssh_tp_dist_min','OIB_MYIflag_mode',
        'OIB_dt_time','OIB_fp_lat','OIB_fp_lon','OIB_n_pts']
    
    #Print the files we found
    print "CS2file: ",cs2file
    print "AWIfile: ",AWIfile
    print "Airborne Files: ",resdata
    
    #Set the output file name. We need to check for existence of AWI/ESA CS2 
    #files. We should always have at least one of these (usually both), so we 
    #use their filename structure to generate the output filename structure.
    if cs2file != []: #If ESA CS2 file exists we use it.
        if completeflag == 1:
            foo = os.path.join(indir,'CVLDAT_'+
                os.path.basename(cs2file[0]).split('_')[1]+"_"+
                os.path.basename(cs2file[0]).split('_')[2]+"_"+
                os.path.basename(cs2file[0]).split('_')[3]+'_'+
                os.path.basename(cs2file[0]).split('_')[4]+'_'+
                os.path.basename(cs2file[0]).split('_')[5]+'_'+
                version+'.csv')
        elif completeflag == 0:
            foo = os.path.join(indir,'CVLDAT_'+
                os.path.basename(cs2file[0]).split('_')[1]+"_"+
                os.path.basename(cs2file[0]).split('_')[2]+"_"+
                os.path.basename(cs2file[0]).split('_')[3]+'_'+
                os.path.basename(cs2file[0]).split('_')[4]+'_'+
                os.path.basename(cs2file[0]).split('_')[5]+'_'+
                version+'_incomplete.csv')
    #If ESA file does not exist but AWI does            
    elif cs2file == [] and AWIfile[0] != []: 
        if completeflag ==1:
            foo = os.path.join(indir,'CVLDAT_'+
                os.path.basename(AWIfile[0]).split('_')[1]+"_"+
                os.path.basename(AWIfile[0]).split('_')[2]+'_'+
                os.path.basename(AWIfile[0]).split('_')[3]+'_'+
                os.path.basename(AWIfile[0]).split('_')[4]+'_'+
                os.path.basename(AWIfile[0]).split('_')[5]+'_'
                +version+'.csv')
        elif completeflag == 0:
            foo = os.path.join(indir,'CVLDAT_'+
                os.path.basename(AWIfile[0]).split('_')[1]+"_"+
                os.path.basename(AWIfile[0]).split('_')[2]+'_'+
                os.path.basename(AWIfile[0]).split('_')[3]+'_'+
                os.path.basename(AWIfile[0]).split('_')[4]+'_'+
                os.path.basename(AWIfile[0]).split('_')[5]+'_'+
                version+'_incomplete.csv')
    
    #Let's start reading in the data files.
    if cs2file !=[]:   #ESA file exists
        cs2=pd.read_table(cs2file[0],sep=',',header=0,index_col=False,
        na_values=['nan'])
    if AWIfile !=[]: #AWI file exists
        AWI = pd.read_table(AWIfile[0],sep=',',header=0,index_col=False,
        na_values=['nan'])
    #Join the two CS2 files
    if AWIfile != [] and cs2file !=[]: #Both AWI and ESA CS2 files exist
        out=pd.concat([cs2,AWI],ignore_index=False,axis=1)
    elif AWIfile ==[] and cs2file !=[]: #ESA only 
        out=pd.concat([cs2],ignore_index=False,axis=1)
    elif AWIfile !=[] and cs2file ==[]: #AWI only
        out=pd.concat([AWI],ignore_index=False,axis=1)	
    	
    #Read in the resampled, subset airborne data.
    for rdata in resdata: #iterate through list of filenames
        data=pd.read_table(rdata,sep=',',header=0,index_col=False,
            na_values=['nan'])
        #If the OIB file is a quicklook, let's rename the headers in this
        #product, the headers in the OIBSIQ file stay the same.  
        if "OIBSIQ" in rdata:
            data.columns=['OIB_lat_mean','OIB_lon_mean','OIB_thickness_mean',
        'OIB_thickness_unc_mean','OIB_mean_fb_mean','OIB_ATM_fb_mean',
        'OIB_fb_unc_mean','OIB_snow_depth_mean','OIB_snow_depth_unc_mean',
        'OIB_ssh_mean','OIB_ssh_sd_mean','OIB_ssh_tp_dist_mean',
        'OIB_thickness_stdev','OIB_thickness_unc_stdev','OIB_mean_fb_stdev',
        'OIB_ATM_fb_stdev','OIB_fb_unc_stdev','OIB_snow_depth_stdev',
        'OIB_snow_depth_unc_stdev','OIB_ssh_stdev','OIB_ssh_sd_stdev',
        'OIB_thickness_median','OIB_thickness_unc_median','OIB_mean_fb_median',
        'OIB_ATM_fb_median','OIB_fb_unc_median','OIB_snow_depth_median',
        'OIB_snow_depth_unc_median','OIB_ssh_median','OIB_ssh_sd_median',
        'OIB_thickness_max','OIB_thickness_unc_max','OIB_mean_fb_max',
        'OIB_ATM_fb_max','OIB_fb_unc_max','OIB_snow_depth_max',
        'OIB_snow_depth_unc_max','OIB_ssh_max','OIB_ssh_sd_max',
        'OIB_ssh_tp_dist_max','OIB_thickness_min','OIB_thickness_unc_min',
        'OIB_mean_fb_min','OIB_ATM_fb_min','OIB_fb_unc_min',
        'OIB_snow_depth_min','OIB_snow_depth_unc_min','OIB_ssh_min',
        'OIB_ssh_sd_min','OIB_ssh_tp_dist_min','OIB_MYIflag_mode',
        'OIB_dt_time','OIB_fp_lat','OIB_fp_lon','OIB_n_pts']
        out=pd.concat([out,data],ignore_index=False,axis=1)
    
    #Now let's add in any missing data records
    if "Latitude" in out.columns: #If ESA exists
        addarr=np.empty((len(out["Latitude"]),addcols))
    elif "AWICS2_lat" in out.columns: #Else if the AWI exists
        addarr=np.empty((len(out["Latitude"]),addcols))
    #Change those fields to nans
    addarr[:]=np.nan
    #Change to dataframe
    addpd = pd.DataFrame(data=addarr,columns=addhead)
    #Add them to the output dataframe
    out=pd.concat([out,addpd],ignore_index=False,axis=1)
    
    #We want the data sorted into ESA,AWI, then airborne sensors alphabetically
    sortorder=['ID','CryoSatDateTime','Year','Month','Day','Hour','Minute',
        'Second','Orbit','Latitude','Longitude','Freeboard','FreeboardFlag',
        'Height_1','Peakiness','Sigma0_1','SHA','Interp_SHA',
        'Interp_OceanHt_Error','Interp_OceanHt_CNT_FWD',
        'Interp_OceanHt_CNT_BKWD','N_Avg','IceConc','SnowDepth','SnowDensity',
        'SAR_DISCRIM','MSS_Mod','Geoid_Mod','ODLE_Mod','DryT_C','WetT_C',
        'IB_C','DAC_C','IONO_GIM','IONO_Mod','H_OT_C','H_LPEOT_C','H_OLT_C',
        'H_SET_C','H_GPT_C','Beam_Center','Beam_Amp','Beam_Kurt','Beam_SD',
        'Beam_SDAng','Beam_Skew','Beam_Center_ang','BadDataFlag',
        'Corr_Error_F','AWICS2_time','AWICS2_lon','AWICS2_lat','AWICS2_fb_stat',
        'AWICS2_fb_syst','AWICS2_fb_unc','AWICS2_fb','AWICS2_sit_stat',
        'AWICS2_sit_syst','AWICS2_sit_unc','AWICS2_sit_fb_syst',
        'AWICS2_sit_snow_syst','AWICS2_sit_snow_syst_iav',
        'AWICS2_sit_rho_i_syst','AWICS2_sit_rho_s_syst','AWICS2_sit',
        'AWICS2_ssha_unc','AWICS2_ssh','AWICS2_ssha','AWICS2_sd',
        'AWICS2_rho_i','AWICS2_rho_s','AWICS2_ice_type',
        'AWICS2_surface_type','AWICS2_ice_conc','AEM_Year','AEM_Month',
        'AEM_Day','AEM_Hour','AEM_Minute','AEM_Seconds','AEM_lon','AEM_lat',
        'AEM_dT','AEM_navg','AEM_Mean','AEM_Median','AEM_stdev','AEM_Mode',
        'ALS_fb_Year','ALS_fb_Month','ALS_fb_Day','ALS_fb_Hour',
        'ALS_fb_Minute','ALS_fb_Seconds','ALS_fb_lon','ALS_fb_lat','ALS_fb_dT',
        'ALS_fb_navg','ALS_fb_Mean','ALS_fb_Median','ALS_fb_stdev',
        'ALS_fb_Mode','ALS_leads','ASR_fb_Year','ASR_fb_Month','ASR_fb_Day',
        'ASR_fb_Hour','ASR_fb_Minute','ASR_fb_Seconds','ASR_fb_lon',
        'ASR_fb_lat','ASR_fb_dT','ASR_fb_navg','ASR_fb_Mean','ASR_fb_Median',
        'ASR_fb_stdev','ASR_fb_Mode','ASR_ssh_Year','ASR_ssh_Month',
        'ASR_ssh_Day','ASR_ssh_Hour','ASR_ssh_Minute','ASR_ssh_Seconds',
        'ASR_ssh_lon','ASR_ssh_lat','ASR_ssh_dT','ASR_ssh_navg','ASR_ssh_Mean',
        'ASR_ssh_Median','ASR_ssh_stdev','ASR_ssh_Mode','ASR_ssha_Year',
        'ASR_ssha_Month','ASR_ssha_Day','ASR_ssha_Hour','ASR_ssha_Minute',
        'ASR_ssha_Seconds','ASR_ssha_lon','ASR_ssha_lat','ASR_ssha_dT',
        'ASR_ssha_navg','ASR_ssha_Mean','ASR_ssha_Median','ASR_ssha_stdev',
        'ASR_ssha_Mode','ASR_leads','OIB_lat_mean','OIB_lon_mean','OIB_thickness_mean',
        'OIB_thickness_unc_mean','OIB_mean_fb_mean','OIB_ATM_fb_mean',
        'OIB_fb_unc_mean','OIB_snow_depth_mean','OIB_snow_depth_unc_mean',
        'OIB_ssh_mean','OIB_ssh_sd_mean','OIB_ssh_tp_dist_mean',
        'OIB_thickness_stdev','OIB_thickness_unc_stdev','OIB_mean_fb_stdev',
        'OIB_ATM_fb_stdev','OIB_fb_unc_stdev','OIB_snow_depth_stdev',
        'OIB_snow_depth_unc_stdev','OIB_ssh_stdev','OIB_ssh_sd_stdev',
        'OIB_thickness_median','OIB_thickness_unc_median','OIB_mean_fb_median',
        'OIB_ATM_fb_median','OIB_fb_unc_median','OIB_snow_depth_median',
        'OIB_snow_depth_unc_median','OIB_ssh_median','OIB_ssh_sd_median',
        'OIB_thickness_max','OIB_thickness_unc_max','OIB_mean_fb_max',
        'OIB_ATM_fb_max','OIB_fb_unc_max','OIB_snow_depth_max',
        'OIB_snow_depth_unc_max','OIB_ssh_max','OIB_ssh_sd_max',
        'OIB_ssh_tp_dist_max','OIB_thickness_min','OIB_thickness_unc_min',
        'OIB_mean_fb_min','OIB_ATM_fb_min','OIB_fb_unc_min',
        'OIB_snow_depth_min','OIB_snow_depth_unc_min','OIB_ssh_min',
        'OIB_ssh_sd_min','OIB_ssh_tp_dist_min','OIB_MYIflag_mode',
        'OIB_dt_time','OIB_fp_lat','OIB_fp_lon','OIB_n_pts']   
    out=out[sortorder]
    
    #Write the output data
    out.to_csv(path_or_buf=foo,sep=',',na_rep='nan',header=True,index=False,
        mode='w',float_format='%f')
    
if __name__=='__main__':
    import os
    import fnmatch
    indir=os.getcwd()
    dirList=[]
    version="A_V001"
    for root,dirnames,filenames in os.walk(indir):
            for filename in fnmatch.filter(filenames,'CS2L2I*'+version+'.csv'):
                dirList.append(os.path.join(root,os.path.dirname(filename)))
    
    for direct in dirList:
        concatenate_output(direct,flag=1,version=version,append=version)
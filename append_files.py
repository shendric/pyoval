# -*- coding: utf-8 -*-
__version__= 1.0
__author__= "Justin Beckers, YorkU/UAlberta"

def append_files(outdir,custom="_concat",fileList,header=0,skiprows=0,sep=',',
        dropna=0):
    '''
    Generic file appending. For a directory, find all files matching certain
    file name characteristics and append the files. Files must be delimited.
    Files may or may not have a header. If the header flag is set, then the
    header from the first file is used for the output.
    
    Author: Justin F. Beckers
    Author Contact: beckers@ualberta.ca
    
    Version: 1.0
    Version Notes: Initial release
    
    Release Date: March 15,2016
    
    Inputs:
        outdir: output directory
        custom: custom string to add into filename. Default is "_concat"
        fileList: list of files to append
        header: Sets whether or not file has a header row or if it should be 
                read. Set to None if you don't want a header to be read, or if
                the file has no header. Set to header row if the header is 
                present and you want it to be read. Default is None.
       skiprows: rows to skip before the header
        sep: file delimiter. Default is comma (','). use '\s+' for space
             delimited files (*.dat).
        dropna: if dropna is set to 1, drop rows that have any nans.
    Outputs:
        Outputs a delimited file with delimiter equal to input file delimiter,
        name equal to the start of the name of the first file with the custom
        string added. Output file has a header if the input header is specified
        or used. Output file is put into outdir directory.
    Usage Notes: 
        Case 1: Using Main 
            From commandline/terminal, can simply run:
                >>>python append_files.py
            This executes the program with the parameters set on 
            lines 110 - 124. Modify this section in the script file text to 
            meet your requirements
                if __name__=='__main__':
        Case 2: Import function:
        import append_files
        append_files(outdir,custom,files,header=0,skiprows=0,dropna=0,sep=',')    
        
    '''   
    #Load the required libraries
    import os
    import pandas as pd
    import numpy as np
    
    #Set the output filename
    outname = os.path.join(outdir,os.path.basename(
            fileList[0]).split('_')[0]+custom+'.csv')
              
    #If files have no header (user specified) or is being ignored:
    if header == None:
        #Read the first file to establish the data record.
        out = pd.read_table(fileList[0],sep=sep,header=None,skiprows=skiprows,
                index_col=None,na_values=['nan','NaN',str(np.nan),np.nan])
        if len(fileList)>1: #If the fileList has more than one file in it.
            #For the rest of the filenames in the fileList
            for x in np.arange(1,len(fileList)):
                data=pd.read_table(fileList[x],sep=sep,header=None,
                    index_col=None,skiprows=skiprows,na_values=['nan','NaN',
                    str(np.nan),np.nan])
                #Concatenate the data files
                out=pd.concat([out,data],ignore_index=True,axis=0)
        #If you want to remove rows where any value is a nan value.
        if dropna==1:
            out=out.dropna()
        '''
        If you need to do anything else to the data, i.e. calculate a new
        column, sort the data, etc., this is the place to insert that code.
        '''
        #Write the output to a new file with the same delimiter as the input.
        #Output file also has no header (like the input files)
        out.to_csv(path_or_buf=outname,sep=sep,header=None,index=False,
            mode='w')
    else: #If the input files have a header (user-specified flag)
        #Read the first file to establish the output data record
        out = pd.read_table(fileList[0],sep=sep,header=header,index_col=None,
                skiprows=skiprows,na_values=['nan','NaN',np.nan])
        #Extract the header in case we need it again further down.
        head=out.columns
        if len(fileList)>1:#If fileList has more than one file.
            #For the rest of the filenames in the fileList
            for x in np.arange(1,len(fileList)):
                data=pd.read_table(fileList[x],sep=sep,index_col=None,
                    na_values=['nan','NaN',np.nan],header=header,
                    skiprows=skiprows)
                #Concatenate the data records.
                out=pd.concat([out,data],ignore_index=True,axis=0)
        #If you want to remove rows where any value is a nan value.
        if dropna==1:
            out=out.dropna()
        '''
        If you need to do anything else to the data, i.e. calculate a new
        column, sort the data, etc., this is the place to insert that code.
        '''
        #Write the output to a new file with the same delimiter as the input.
        #Output file also has a header (like the input files)
        out.to_csv(path_or_buf=outname,sep=sep,na_rep='nan',header=True,
                    index=False,mode='w')
        
if __name__=='__main__':
    import glob
    import time #Used to roughly time the run time of the code.
    
    #The setup.
    indir='/Volumes/JustinBeckers_2TBPortable/CryoSat_Lakes/PARLSAR_CS2Data/'
    files=glob.glob(indir+"*_positions.csv")
    outdir='/Volumes/JustinBeckers_2TBPortable/CryoSat_Lakes/PARLSAR_CS2Data'
    custom='concat'
  
    #Run.
    a=time.clock()
    append_files(outdir,custom,files,header=0,skiprows=0,dropna=0,sep=',')
    b=time.clock()
    print "Time to run the script (s): ",b-a
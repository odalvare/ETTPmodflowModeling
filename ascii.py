# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:19:43 2019

@author: The Bodrio

Object for reading ascii files, used for creating the ibound array for 
MODFLOW execution.
"""

import numpy as np
import sys
import logging

class ascii:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,workspace,name_file):
        """
        Construction of an instance for reading and writing ascii file
        
        Parameters:
        ----------
            workspace : str
                The folder location of the ascii file (private)
            name : str
                The name of the ascii file, without extension      
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=workspace
        self.__name=name_file
        self.__ncol=int(0)
        self.__nfil=int(0)
        self.__xia=float(0.0)
        self.__yia=float(0.0)
        self.__res=float(0.0)
        self.__ndata=float(-9999.0)        
    #Ending of constructor ----------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def write(self,mtrx,res,pres,xia,yia,ndata):
        """
        Method for writing raster data on an ascii file

        Parameters:
        ----------
            mtrx : float
                2D matrix of data
            res : float
                Spatial extent of each cell on the raster (resolution)
            pres : integer
                Number of decimal places used when printing the number    
            xia : float
                Western coordinate of the raster boundary
            yia : float
                Southern coordinate of the raster boundary
            ndata : float
                No data value
        
        Returns:
        -------
            null
            
        """        
        #code begins here
        n=np.size(mtrx,0)
        m=np.size(mtrx,1)
        
        fl=open(self.__workspace+'\\'+self.__name+'.asc','w')
        
        fl.write('NCOLS         '+str(m))
        fl.write('\n')
        fl.write('NROWS         '+str(n))
        fl.write('\n')
        fl.write('XLLCORNER     '+str(xia))
        fl.write('\n')
        fl.write('YLLCORNER     '+str(yia))
        fl.write('\n')
        fl.write('CELLSIZE      '+str(res))
        fl.write('\n')
        fl.write('NODATA_VALUE  '+str(ndata))
        fl.write('\n')    
        
        for i in range(0,n,1):
            for j in range(0,m,1):
                fl.write(str(round(mtrx[i,j],pres))+' ')
            fl.write('\n')
        
        fl.close     
        print('*---El archivo '+self.__name+'.asc se escribio correctamente---*')
        
    #Ending of method ---------------------------------------------------------
    
    #Begining of method -------------------------------------------------------    
    def header(self): 
        """
        Method for reading the header of an ascii file

        Parameters:
        ----------        
            null
        
        Returns:
        -------        
            (col,nfil,xia,yia,res,ndata) : tuple
                Header of *.asc file (ind=1).           
        """ 
        
        #code begins here
        M=np.loadtxt(fname=self.__workspace+'\\'+self.__name+'.asc',dtype='str',
                     comments='#',delimiter=None,converters=None,skiprows=0,
                     usecols=(0,1),unpack=False,ndmin=0)                     
            
        self.__ncol=int(M[0,1])
        self.__nfil=int(M[1,1])
        self.__xia=float(M[2,1])
        self.__yia=float(M[3,1])
        self.__res=float(M[4,1])
        self.__ndata=float(M[5,1])
        
        print('*---Succesfull reading of the header in file '+self.__name+'.asc---*')
        
        return self.__ncol,self.__nfil,self.__xia,self.__yia,self.__res,self.__ndata            
    #Ending of method ---------------------------------------------------------
        
    #Begining of method -------------------------------------------------------    
    def readmtrx(self,tip,filini=6): 
        """
            Method for reading raster data on an ascii file

        Parameters:
        ----------        
            tip : integer
                Type of data used for the number representation, as follows:
                    0: Short integer.
                    1: Long integer.
                    2: Single float precision
                    3: Double float presicion
            filini : integer
                Row of the *.asc file on which matrix data reading starts
            ind : integer
                Indicarot of which data are requested from the reading:
                    0: Matrix data only
                    1: *.asc files header only
                    2: Both, matrix and header
        
        Returns:
        -------
            M : array
                Matrix of data (ind=0)        
        """ 
        
        #code begins here
        if tip==0:
            M=np.loadtxt(fname=self.__workspace+'\\'+self.__name+'.asc',dtype='int32',
                         comments='#',delimiter=None,converters=None,
                         skiprows=filini,usecols=None,unpack=False,ndmin=0)                        
            print('*---Succesfull reading of the matrix in file '+self.__name+'.asc---*')
        elif tip==1:
            M=np.loadtxt(fname=self.__workspace+'\\'+self.__name+'.asc',dtype='int64',
                         comments='#',delimiter=None,converters=None,
                         skiprows=filini,usecols=None,unpack=False,ndmin=0)              
            print('*---Succesfull reading of the matrix in file '+self.__name+'.asc---*')
        elif tip==2:
            M=np.loadtxt(fname=self.__workspace+'\\'+self.__name+'.asc',dtype='float32',
                         comments='#',delimiter=None,converters=None,
                         skiprows=filini,usecols=None,unpack=False,ndmin=0)              
            print('*---Succesfull reading of the matrix in file '+self.__name+'.asc---*')
        elif tip==3:
            M=np.loadtxt(fname=self.__workspace+'\\'+self.__name+'.asc',dtype='float64',
                         comments='#',delimiter=None,converters=None,
                         skiprows=filini,usecols=None,unpack=False,ndmin=0)              
            print('*---Succesfull reading of the matrix in file '+self.__name+'.asc---*')
        else:
            print('*---Invalid type of number ('+str(tip)+')---*')  
        
        mtrx=M
        return mtrx            
    #Ending of method ---------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def checkprmtrs(self,nc,nf,xi,yi,r,miss):
        """
        Method for checking if the headers of the ascii equals a set of parameters
        
        Parameters:
        ----------        
            nc : integer
                Number columns for checking.
            nf : integer
                Number of rows for checking.
            xi: float 32 bits
                x coordinate of the western border of the ascii.
            yi : float 32 bits
                y coordinate of the southern border of the ascii.
            r : float 32 bits
                Spatial resolution for checking.
            miss : float 32 bits
                Missing value for checking.
           
        Returns:
        ------
            equal : boolean
                False if the parameters are different of the actual ascci.
            
        Exception:
        ------
            RuntimeError: inconsistent parameters.
        """
        
        #code begins here
        logging.basicConfig(level=logging.WARNING)
        log = logging.getLogger('In ascii checkprmtrs()')
        try:
            equal = False
            if (nc!=self.__ncol):
                equal = True
                print("*---Error in number of columns---*")                
            if(nf!=self.__nfil):
                equal = True                
                print("*---Error in number of rows---*")                
            if(xi!=self.__xia):
                equal=True
                print("*---Error in y coordinate---*")                
            if(yi!=self.__yia):
                equal=True                
                print("*---Error in y coordinate---*")                
            if(r!=self.__res):
                equal=True                
                print("*---Error in spatial resolution---*")                
            if(miss!=self.__ndata):                
                equal=True                
                print("*---Error in missing value---*")                
            if(equal==True):
                print('*---Check metadata of file '+self.__workspace+'\\'+self.__name+'.asc---*')
                self.throwsasc()
        except Exception:
            log.exception('Error from throws():')
            sys.exit()
        
        print('*---Succefull checking of ascii file---*')
        return equal
    #Ending of method ---------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def throwsasc(self):
        """Method that throws a RuntimeException for checking the parameters of the ascii map
        
        Parameters
        ----------
        
        Raises
        ----------
            RuntimeError: Unsuccefull checking of ascii file
        
        """
        
        #code begins here
        raise RuntimeError('Unsuccefull checking of ascii file')
    #Ending of method ---------------------------------------------------------
        
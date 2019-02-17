# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:52:42 2019

@author: alvareo
"""

import flopy
import numpy as np
from Sgemsimp import Sgemsimp
from rundata import rundata
from configBAS import bldBAS
from configWELL import bldWELL

class runMF:
    
    #Beginning of constructor -------------------------------------------------
    def __init__(self,exefile,workspace,dataspace,rea,thrs,csim):
        """
        Construction of an instance for running MODFLOW

        Parameters
        ----------
            workspace : string
                The folder location of the modeling space
             dataspace : string
                The folder location of the data space, where running info is saved      
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__exefile=exefile
        self.__workspace=workspace
        self.__dataspace=dataspace
        
        #Acquisition of data for running MODFLOW
        self.__rd=rundata(dataspace)
        self.__modelname,Lx,Ly,Lz,zbot,nlay,delcx,delry,nper,steps,Lt,grid,wlay=self.__rd.readdata()
        self.__grid=grid        
        
        #This is the Modflow model
        self.__mf = flopy.modflow.Modflow(self.__modelname, exe_name=self.__exefile, model_ws=self.__workspace)
        
        #Discretization module
        nrow=self.__rd.getNrows()
        ncol=self.__rd.getNcols()
        self.__dis,nstp,perlen=self.__rd.createdis(self.__mf)                  #nstp is an array wit the time steps
        self.__nper=nper
        self.__steps=steps                                                     #Constant time step
        self.__perlen=perlen
        self.__nstp=nstp
        
        #BAS package
        basb=bldBAS(self.__dataspace,self.__grid)
        self.__bas=basb.createbas(self.__mf)
        
        #Forcing packages
        r = 1.0*np.ones((nrow, ncol), dtype=np.float32)/1000000.0              #Recharge
        self.__rcg = flopy.modflow.mfrch.ModflowRch(self.__mf,rech=r)
        self.__wellbld=bldWELL(dataspace,self.__grid,self.__rd,csim)           #Wells
        self.__wel=self.__wellbld.createwel(self.__mf)
        
        # Add LPF package to the MODFLOW model
        layt=[0]        #Type of layer
        sg=Sgemsimp(self.__dataspace,self.__grid)
        hkd = 200.0*np.ones((nlay, nrow, ncol), dtype=np.float32)
        hkd=sg.readsgems1rea(rea)
        hkd[hkd<thrs]=thrs
        vkd = hkd
        self.__lpf = flopy.modflow.ModflowLpf(self.__mf,laytyp=layt,hk=hkd,vka=vkd,ipakcb=53)
        
        # PCG package for matrix computation
        self.__pcg = flopy.modflow.ModflowPcg(self.__mf, mxiter=200, iter1=80)
        
        # Add OC package to the MODFLOW model
        spdata = {}
        for i in range(0,nper,1):
            for j in range(0,nstp[i],1):
                spdata[i,j]=['save head','save budget']
        self.__oc = flopy.modflow.ModflowOc(self.__mf, stress_period_data=spdata, compact=True)
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def runMODFLOW(self,nmflt):
        """
        Runs a simulation of MODFLOW

        Parameters:
        ----------        
            nmflt: String
                Location to save the connection with MT#D
        
        Returns:
        -------
            null
        """
        # linkage to mt3dms LMT package (ATTENTION!!!)
        lmt = flopy.modflow.ModflowLmt(self.__mf, output_file_name=nmflt)
        # Write the MODFLOW model input files
        self.__mf.write_input()
        # Run the MODFLOW model
        success, buff = self.__mf.run_model()  #Flow processes are finished
        
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getGrid(self): 
        """
        Get the grid of the model
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            grid : dictionary
                Contains the information of the domain's geometry (extracted form read data)
        """
        #code begins here
        return self.__grid
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getMfName(self): 
        """
        Get the name of the model
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            modelname : string
                Returns the name of the model
        """
        #code begins here
        return self.__modelname
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getPerlen(self): 
        """
        Get the stress period duration
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            perlen : array float
                Duration of each forcing period
        """
        #code begins here
        return self.__perlen
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNper(self): 
        """
        Get the number of stress periods
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nper : int
                Number of strss periods
        """
        #code begins here
        return self.__nper
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getNstp(self): 
        """
        Get the number of steps inside each stress period
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            nstp : int
                Number of steps on each stress period
        """
        #code begins here
        return self.__nstp
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getWelBld(self): 
        """
        Get the objectes for reading well information
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            wellbld : configWELL object
                Object for reading the well information
        """
        #code begins here
        return self.__wellbld
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------  
    def getModel(self): 
        """
        Get the mf model object.
        
        Parameters
        ----------        
            null
        
        Returns
        -------        
            mf : flopy.modflow.Modflow mf
                mf object of the simuation
        """
        #code begins here
        return self.__mf
    #End of method ------------------------------------------------------------
    
    
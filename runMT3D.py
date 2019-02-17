# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 18:11:30 2019

@author: alvareo
"""

import flopy
from runMF import runMF

class runMT3D:
    
    #Beginning of constructor -------------------------------------------------
    def __init__(self,mf,rnmf,exemt3dms,modelname,workspace,dataspace,fltfile):
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
        self.__exefile=exemt3dms
        self.__workspace=workspace
        self.__dataspace=dataspace
        self.__modelname=modelname
        self.__ftlnm=fltfile
        
        # create mt3dms model object
        self.__mt = flopy.mt3d.Mt3dms(modflowmodel=mf, modelname=self.__modelname, 
                               exe_name=self.__exefile , ftlfilename=self.__ftlnm,
                               model_ws=self.__workspace)
        
        # basic transport package
        perlen=rnmf.getPerlen()
        nper=rnmf.getNper()
        nstp=rnmf.getNstp()
        self.__btn = flopy.mt3d.Mt3dBtn(self.__mt, prsity=0.3, icbund = 1, sconc=0.0, ncomp=1, 
                             perlen = perlen, nper=nper, nstp = nstp, tsmult = 1.0, 
                             nprs = -1, nprobs = 10, cinact = -1, chkmas=True)
        
        # advaction package
        self.__adv = flopy.mt3d.Mt3dAdv(self.__mt, mixelm=-1, percel=0.75)
        # dispersion package
        self.__dsp = flopy.mt3d.Mt3dDsp(self.__mt, al=1000, trpt=0.1, trpv=0.1, dmcoef=1e-09)
        #Forcing objects
        itype = flopy.mt3d.Mt3dSsm.itype_dict()
        print(itype)
        print(flopy.mt3d.Mt3dSsm.get_default_dtype())
        wellbld=rnmf.getWelBld()
        self.__ssm=wellbld.createssm(self.__mt,itype)
        # matrix solver package
        self.__gcg = flopy.mt3d.Mt3dGcg(self.__mt, cclose=1e-6)
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def runMT3D(self):
        """
        Runs a simulation of MT3D

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            null
        """
        
        # write mt3dms input
        self.__mt.write_input()
        # run mt3dms
        self.__mt.run_model()
        
    #End of method ------------------------------------------------------------
    
    
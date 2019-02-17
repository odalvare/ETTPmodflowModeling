# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 17:39:29 2018

@author: The Bodrio
"""

import numpy as np

"""
*******************************************************************************
This routine writes a 3D array of data in a *.dat geoeas file.
"""
def wr_3Dar(hdr,mtrx,var,rx,ry,rz,pres,xia,yia,zia,ndata,rt,name):
    fl=open(rt+'\\'+name+'.dat','w')
    fl.write(hdr);fl.write('\n')            #Writing header of *.dat file
    fl.write('4');fl.write('\n')            #Writing the number of variables (x,y,h,z)
    fl.write('Easting');fl.write('\n')      #Writing East coordinates
    fl.write('Northing');fl.write('\n')     #Writing West coordinates
    fl.write('Altitude');fl.write('\n')     #Writing Altitude coordinates
    fl.write(str(var));fl.write('\n')       #Writing name of variable
    
    l=np.size(mtrx,0)   #Layers
    n=np.size(mtrx,1)   #Rows
    m=np.size(mtrx,2)   #Columns    
    
    cont=0
    for k in range(0,l-1,1):
        for i in range(0,n-1,1):
            for j in range(0,m-1,1):
                cont=cont+1
                z=1.*(zia+rz/2)+(1.0*k)*rz
                y=1.*(yia+ry/2)+(1.0*i)*ry
                x=1.*(xia+rx/2)+(1.0*j)*rx
                fl.write(str(round(x,pres))+'   ')
                fl.write(str(round(y,pres))+'   ')
                fl.write(str(round(z,pres))+'   ')
                fl.write(str(round(mtrx[k,i,j],pres))+'   ')
                fl.write('\n')
    fl.close     
    print('*---The file '+name+'.asc was succesfully written---*')

"""
*******************************************************************************
This routine reads a 3D array of data fron a *.dat geoeas file.
""" 
def rd_3Dar(rt,name,tp):
    
    #Read header of *.dat file
    M=np.loadtxt(fname=rt+'\\'+name+'.dat',dtype='str',comments='#',
                 delimiter=None,converters=None,skiprows=0,usecols=None,
                 unpack=False,ndmin=0)
    hdr=np.char.lower(M[0,0])
    nvrbl=int(M[1,0])
    nm1=np.char.lower(M[2,0])
    nm2=np.char.lower(M[3,0])
    nm3=np.char.lower(M[4,0])    
    cols=int(M[2,1])
    rows=int(M[3,1])
    lays=int(M[4,1])
    
    #Read matriz of coordinates and data
    if tp==0:
        M=np.loadtxt(fname=rt+'\\'+name+'.dat',dtype='int32',
                     comments='#',delimiter=None,converters=None,
                     skiprows=nvrbl+2,usecols=None,unpack=False,ndmin=0)        
    elif tp==1:
        M=np.loadtxt(fname=rt+'\\'+name+'.dat',dtype='int64',
                     comments='#',delimiter=None,converters=None,
                     skiprows=nvrbl+2,usecols=None,unpack=False,ndmin=0)
    elif tp==2:
        M=np.loadtxt(fname=rt+'\\'+name+'.dat',dtype='float32',
                     comments='#',delimiter=None,converters=None,
                     skiprows=nvrbl+2,usecols=None,unpack=False,ndmin=0)
    elif tp==3:
        M=np.loadtxt(fname=rt+'\\'+name+'.dat',dtype='float64',
                     comments='#',delimiter=None,converters=None,
                     skiprows=nvrbl+2,usecols=None,unpack=False,ndmin=0)
    else:
        print('*---Invalid type of number ('+str(tp)+')---*')
    
    #Build the 3D array with data
    minar=np.amin(M, axis=0)
    maxar=np.amax(M, axis=0)e
    minx=minar[0]
    print('minx='+str(minx))
    maxx=maxar[0]
    print('maxx='+str(maxx))
    miny=minar[1]
    print('minx='+str(miny))
    maxy=maxar[1]
    print('maxx='+str(maxy))
    
    return hdr,nvrbl,nm1,nm2,nm3
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 10:12:44 2018

@author: cecidip
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_rnSGS(g,E,):
    """
    Plots slices of randomly chosen ensemble members
    If you wanna save the plot call the funtion as:
    fig=plot_rnSGS(g,E)
    fig.savefig("somefile.png")
    """
    x=g.x
    y=g.y
    z=g.z

    N=len(E) #
    Vrn1 = E[np.random.randint(N),:].reshape(g.shape)
    Vrn2 = E[np.random.randint(N),:].reshape(g.shape)
    Vrn3 = E[np.random.randint(N),:].reshape(g.shape)
    Vrn4 = E[np.random.randint(N),:].reshape(g.shape)
    
    fig=plt.figure(figsize=(10,8))
    plt.subplot(221)
    plt.pcolor(x,z,np.squeeze(Vrn1[:,10,:].T), cmap='jet'), plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.colorbar()
    plt.subplot(222)
    plt.pcolor(x,z,np.squeeze(Vrn2[:,10,:].T), cmap='jet'), plt.gca().invert_yaxis()
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.colorbar()
    plt.subplot(223)
    plt.pcolor(x,z,np.squeeze(Vrn3[:,10,:].T), cmap='jet'), plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    plt.subplot(224)
    plt.pcolor(x,z,np.squeeze(Vrn4[:,10,:].T), cmap='jet'), plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    plt.suptitle('rnV_Ee')
    plt.show()
    return fig

def esfera(g,centro,V,pp=0.05, rr=120):
    #centro: centre of the sphere
    #rr:radius of the sphere
    #pp: % of Velocity change 
    #V: velocity model
    corr= np.ones((g.shape))
    mm=-pp/rr
                         
    x0=centro[0]
    y0=centro[1]
    z0=centro[2]             
    #x0=sum(hyp[:,2])/nev
    #y0=sum(hyp[:,3])/nev
    #z0=sum(hyp[:,4])/nev
    for i in range(0,g.x.size):
            for j in range(0,g.y.size):
                for k in range(0,g.z.size):
                    x=g.x[i]
                    y=g.y[j]
                    z=g.z[k]
                    dis=np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2) 
                    
                    if dis<=rr:
                        corr[i,j,k]= (1+pp)#+mm*dis#+4*np.random.randn(1))           
    corrV=corr.flatten()
    #g.toXdmf(corrV, 'esfera2', 'esfera2')
    return V*corrV

def stats(E,Vtrue):
    #Tarrahi2015:"Integration of microseismic monitoring data ..Kalman filter"
    N,m=E.shape #num of esemble members
    tmp=np.sum((E-Vtrue)**2,axis=0)
    RMSE=np.sum(np.sqrt(tmp/N))/m 
    vmean=np.sum(E,axis=0)/N
    tmp1=np.sum((E-vmean)**2,axis=0)
    Sp=np.sum(np.sqrt(tmp1/N))/m
    return (RMSE,Sp)
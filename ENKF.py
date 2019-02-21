#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 16:52:25 2018

@author: cecidip

The EnKF.
  Ref: Evensen, Geir. (2009):
  "Data assimilation: the ensemble Kalman filter}."
  
"""
import numpy as np
import seaborn as sns
sns.set_style("whitegrid")
import hypo #https://github.com/groupeLIAMG/hypopy/blob/master/hypo.py
import multiprocessing as mp
from numpy.random import multivariate_normal
from common import * # https://github.com/nansencenter/DAPPER
from functions import *
from tqdm import tqdm #https://github.com/tqdm/tqdm/tree/master/tqdm


def fore_step(E, hyp, rcv):

    """FORECAST STEP"""
    N,m=E.shape #num of enseble members
    #centro=np.mean(ms[:,2:],axis=0) #basado en los ms
    centro=ms[2:]
    hyp = np.kron(hyp,np.ones((nsta,1)))
    rcv_data = np.kron(np.ones((nev,1)), rcv)
    Ef=np.zeros((E.shape))
    hE=np.zeros((N,len(tt))) # Compute times with inverted hypocentre locations (hE)    
    
    for j in tqdm(range(N)):
        Ef[j,:]=esfera(g,centro,E[j,:],pp=np.random.uniform(-0.1,0.1),rr=np.random.uniform(0,60))
        slowness = 1./Ef[j,:]
        hE[j,:] = g.raytrace(slowness, hyp, rcv_data)-hyp[:,1]
    
    return Ef,hE,N,slowness,rcv_data


#%%
def analysis_step(Ef, hE, R, N):
    
    """ANALYSIS STEP"""  
    # Reference:  P. Raanes 
    # https://github.com/nansencenter/DAPPER/blob/master/da_methods/da_methods.py

    mu = mean(Ef,0) #ensemble mean
    A  = Ef - mu #enesemble anomaly
    hx = mean(hE,0)
    Y  = hE-hx #anomalies of the observed ensemble
    D = multivariate_normal([0]*len(tt), R, N) #https://github.com/rlabbe/Kalman-anYd-Bayesian-Filters-in-Python/blob/master/Appendix-E-Ensemble-Kalman-Filters.ipynb
    C  =  Y.T @ Y +( R*(N-1))
    YC =  mrdiv(Y, C)
    KG = A.T @ YC         
    dE = (KG @ ( tt + D - hE ).T).T 
    Ea  = Ef + dE 
    
    return Ea

#%% 
def tt_error(Ea,slowness,hyp,rcv_data,hE,N):

    """CONTROL INVERSION - compute tt with analysed models to check if they improved"""
    hEa=np.zeros((N,len(tt)))
    for j in tqdm(range(N)):
        slowness = 1./Ea[j,:]
        hEa[j,:] = g.raytrace(slowness, hyp, rcv_data)-hyp[:,1]
            
    err_for = abs(tt- np.mean(hE,axis=0))
    err_ana = abs(tt- np.mean(hEa,axis=0))
    
    deltaTT=plt.figure()
    plt.plot(err_for,'o-')#,label=r'$\|\|\Delta tt fore\|\|$ = {0:6.5f}'.format(np.mean(err_for)))
    plt.plot(err_ana,'r*-')#,label=r'$\|\|\Delta tt analysed\|\|$ = {0:6.5f}'.format(np.mean(err_ana)))
    plt.ylabel('tt',fontsize=18)
    plt.xlabel('Ray num',fontsize=18)
    plt.tick_params(labelsize=15)
    plt.legend()
    return hEa

#%%
if __name__ == '__main__':

    ##  GRID
    xmin = 9600
    xmax = 10440
    ymin = 10036
    ymax = 10416
    zmin = 2626
    zmax = 3206
    dx = 20  # grid cell size, we use cubic cells here
    x = np.arange(xmin, xmax+1, dx)
    y = np.arange(ymin, ymax+1, dx)
    z = np.arange(zmin, zmax+1, dx)
    nthreads= mp.cpu_count() #use all the CPUs for parallel task
    g = hypo.Grid3D(x, y, z, nthreads)
    
    step_i =1 #change for each step
    step = str(step_i)
    step_prev = str(step_i-1)
    
    #%%
    """ imput files """
    ##  Receivers
    rcv=np.loadtxt(open("rcvRND.txt"), skiprows=1)
    ## Observed times  
    tt=np.loadtxt(open("test_ranREC/tt_synth" + step + ".txt"))
    ## Hypos
    hyp=np.loadtxt(open("test_ranREC/h_true" + step + ".txt"))
    ## Previews (or initial) V models, by SGS :
    #E=np.loadtxt(open("E100.txt")).T #ensemble de V
    E=np.loadtxt(open("test_ranREC/Ea"+step_prev+".txt"))
    # true V for RMSE
    Vtrue=np.loadtxt(open("test_ranREC/Vtrue" + step + ".txt"))
    # microseisms coordinates
    ms=np.loadtxt(open("test_ranREC/hinit" + step + ".txt")) 
    
    #%%
    ircv = np.arange(rcv.shape[0]).reshape(-1,1)   # vector of rcv indices
    nsta = rcv.shape[0]
    nev=1
    #nev=len(hyp)
    ## Errors:
    std_noise = 1.e-3
    R=np.eye((len(tt))) * std_noise**2
   
    #%%
    """Perform EnKF"""
    Ef,hE, N , slowness, rcv_data= fore_step(E, hyp, rcv)               #forecast step
    Ea= analysis_step(Ef, hE, R, N)                 #Analysis step
    hEa= tt_error(Ea,slowness, hyp,rcv_data,hE,N)   #Control travel times
    #%% 
    
    """ SAVE RESULTS"""
    
    #slices of E before the forecast
    fig=plot_rnSGS(g,E)
    fig.savefig("test_ranREC/rnEa" + step + ".png")
    #slices de Ef
    fig=plot_rnSGS(g,Ef)
    fig.savefig("test_ranREC/rnEf" + step + ".png")
    #slices E
    fig=plot_rnSGS(g,Ea)
    fig.savefig("test_ranREC/rnEa" + step + ".png")
    #analysed ensembles Ea
    np.savetxt("test_ranREC/Ea" + step + ".txt", Ea)
    #plots Ef Ea and dif
    g.toXdmf(np.mean(Ef,axis=0), "Vp_forcasted" + step ,  "test_ranREC/Vp_forcasted" + step )
    g.toXdmf(np.mean(Ea,axis=0), "Vp_analysed" + step,  "test_ranREC/Vp_analysed" + step)
    g.toXdmf(np.abs(np.mean(Ea,axis=0)-np.mean(Ef,axis=0)), "Vpa-Vpf" + step ,  "test_ranREC/Vpa-Vpf" + step)
    
    # tt analysed hEa and forecasted hE
    np.savetxt("test_ranREC/hE" + step + ".txt", hE)
    np.savetxt("test_ranREC/hEa" + step + ".txt", hEa)
    deltaTT.savefig("test_ranREC/deltaTT" + step )
        

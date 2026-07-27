import sys
import os
import time
import numpy as np
import obspy
from numpy.random import multivariate_normal
from joblib import Parallel, delayed
from tqdm import tqdm

sys.path.append(os.path.abspath("C:/Users/marti/hypopy"))
import hypo
from functions import esfera 

def read_seg(foldername, h):
    """
    Read events per hour.
    Input:
        h = hour of the day (0-24) 
        foldername = location of the files of interest
    Output:
        array with: RCV (x-y-z), SRC (x-y-z), timestamp, P-pick, S-pick
    """
    rcx, rcy, rcz = [], [], []        
    srx, sry, srz = [], [], []
    ttp, tts, timestamp = [], [], []
    
    for filename in os.listdir(foldername):
        if filename.lower().endswith(".sgy"):
            # Modern ObsPy implementation for reading SEGY
            stream = obspy.read(os.path.join(foldername, filename), format='SEGY', unpack_trace_headers=True)
            
            for tr in stream:
                header = tr.stats.segy.trace_header
                if header.hour_of_day == h:
                    timestamp.append(header.hour_of_day * 3600 + header.minute_of_hour * 60 + header.second_of_minute)
                    
                    # XYZ receivers
                    rcx.append(header.group_coordinate_x / 100.0)
                    rcy.append(header.group_coordinate_y / 100.0)
                    rcz.append(header.datum_elevation_at_receiver_group / 100.0)
                    
                    # XYZ source
                    srx.append(header.source_coordinate_x / 100.0)
                    sry.append(header.source_coordinate_y / 100.0)
                    srz.append(header.source_depth_below_surface / 100.0)
                    
                    # P and S pick
                    ttp.append(header.x_coordinate_of_ensemble_position_of_this_trace)   
                    tts.append(header.y_coordinate_of_ensemble_position_of_this_trace)               
    
    ret = np.array([rcx, rcy, rcz, srx, sry, srz, timestamp, ttp, tts]).T
    
    if len(ret) > 0:
        return ret[ret[:, 6].argsort()] # Sort by timestamp
    return ret

def inv_vcte(data_4inv, k, g): 
    """
    Locate hypocenter with constant velocity.
    Inputs:
        data_4inv: file containing EvID#, tt[s], rcvID, rcv X, rcv Y, rcv Z, 0/1 for P/S
        k: Event ID to isolate
        g: Grid object
    """
    DATA = np.copy(data_4inv)
    nev = 1
    
    # Generate initial guess (random) of hypo location
    hinit = np.vstack((
        np.arange(nev),
        np.linspace(0., 50., nev),
        g.x[50] + 10 * np.random.randn(nev),
        g.y[40] + 10 * np.random.randn(nev),
        g.z[40] + 10 * np.random.randn(nev)
    )).T                  

    dat = []
    rcv = []
    
    for i in range(DATA.shape[0]):  
        if DATA[i, 0] == k:
            DATA[i, 0] = 0
            rcv.append(DATA[i, 4:7]) 
            dat.append(DATA[i, 0:4])
            
    # Use standard Python int instead of deprecated np.int
    nsta = int(len(rcv) / 2)
    
    hinit1, res = hypo.hypolocPS(
        np.array(dat), 
        np.array(rcv)[:nsta], 
        V=np.array((6030., 3523.)), 
        hinit=hinit, 
        maxit=15, 
        convh=10., 
        verbose=False
    )                       
    
    return hinit1, np.array(dat)[:nsta, 1], np.array(rcv)[:nsta], nsta    

# --- Parallel Helper Functions ---
def perturb(i, g, centro, Ej, pp_min, pp_max, rr_min, rr_max):
    """ 
    Helper function for parallel spherical perturbation
    """    
    # Calculate dynamic random variables
    if pp_min == pp_max:
        pp_val = pp_min
    else:
        pp_val = np.random.uniform(pp_min, pp_max)
        
    rr_val = np.random.uniform(rr_min, rr_max)
    
    # Apply the sphere (assuming vlim=2500 is only strictly needed for real data, 
    # but safe to leave as a generic limiter or set to None)
    Ef_row = esfera(g, centro, Ej, pp=pp_val, rr=rr_val, vlim=None)
    
    return i, Ef_row

def raytr(j, g, rcv, ms, slowness):
    try:
        # 1. Extract the Origin Time (t0) from column 1
        t0 = ms[:, 1]
        
        # 2. Extract only X, Y, Z (columns 2, 3, 4) and force them to be clean 2D arrays
        source_coords = np.ascontiguousarray(ms[:, 2:5])
        rcv_coords = np.ascontiguousarray(rcv)
        
        # 3. Modern ttcrpy signature: source first, then receivers, then slowness
        tt_travel = g.raytrace(source_coords, rcv_coords, slowness=slowness)  
        
        # 4. The actual Arrival Time is the travel time across the grid + the origin time
        arrival_time = tt_travel + t0
        
        return j, arrival_time
        
    except Exception as e:
        print(f"An exception occurred in iteration {j} : {e}")
        return j, None
# ---------------------------------

def fore_step(E, ms, rcv, nsta, nev, g, pp_min, pp_max, rr_min, rr_max, ms_true=None, n_jobs=-1):
    """
    FORECAST step using modern Joblib parallelization.
    with live tqdm progress bars.
    """
    starttime = time.time()
    N, m = E.shape 
    centro = ms[:, 2:].T
    
    # Raytrace from the TRUE location if provided (Synthetic), otherwise use ms (Real)
    raytrace_source = ms_true if ms_true is not None else ms

    Ef = np.zeros(E.shape) 
    hE = np.zeros((N, len(rcv)))   
    
    ms_kron = np.kron(raytrace_source, np.ones((nsta, 1)))
    rcv_kron = np.kron(np.ones((nev, 1)), rcv)     
    
    # 1. Parallel Slowness Perturbation
    # Generates a quick progress bar for the perturbation step
    jobs_perturb = (delayed(perturb)(i, g, centro, E[i, :], pp_min, pp_max, rr_min, rr_max) for i in range(N))
    Ef_para = []
    
    for result in tqdm(Parallel(n_jobs=n_jobs, return_as="generator")(jobs_perturb), total=N, desc="     Perturbing", ncols=80):
        Ef_para.append(result)
        
    for Efp in Ef_para:
        Ef[Efp[0], :] = Efp[1]  
        
    # 2. Parallel Raytracing
    # Generates the main progress bar for the heavy C++ raytracing step
    jobs_raytr = (delayed(raytr)(i, g, rcv_kron, ms_kron, 1.0 / Ef[i, :]) for i in range(N))
    results = []
    
    for result in tqdm(Parallel(n_jobs=n_jobs, return_as="generator")(jobs_raytr), total=N, desc="     Raytracing", ncols=80):
        results.append(result)
        
    for result in results:
        if result[1] is None:
            raise Exception(f"Raytracing failed at ensemble member {result[0]}")
        hE[result[0], :] = result[1]  

    print(f'     Forecast step took {time.time() - starttime:.2f} seconds')  
    return Ef, hE, rcv_kron, ms_kron

def mrdiv(A, B):
    """ Matrix right division: A / B = A @ inv(B) 
    Computes the matrix product A @ inv(B) in a numerically stable way
    using np.linalg.solve instead of explicit inversion.

    Parameters
    ----------
    A : np.ndarray
        Left-hand matrix.
    B : np.ndarray
        Right-hand matrix.

    Returns
    -------
    X : np.ndarray
        Result of the right matrix division (A / B).

    """
    return np.linalg.solve(B.T, A.T).T

def analysis_step(Ef, hE, R, tt):    
    """
    ANALYSIS STEP 
    Reference: P. Raanes (DAPPER)
    """  
    N, m = Ef.shape
    mu = np.mean(Ef, 0)  # ensemble mean
    A = Ef - mu          # ensemble anomaly
    hx = np.mean(hE, 0)
    Y = hE - hx          # anomalies of the observed ensemble
    
    D = multivariate_normal([0]*len(tt), R, N) 
    
    C = Y.T @ Y + (R * (N - 1))
    YC = mrdiv(Y, C)
    KG = A.T @ YC         
    
    dE = (KG @ (tt + D - hE).T).T 
    Ea = Ef + dE 

    return Ea
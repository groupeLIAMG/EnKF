# main.py
import numpy as np
import pandas as pd
import math
import os
import time
import sys

# Import our custom toolboxes
sys.path.append(os.path.abspath("C:/Users/marti/hypopy"))

import config
import enkf_core as enkf
import hypo

np.linalg.linalg = np.linalg # Forces modern Numpy to accept the deprecated .linalg.linalg calls
'''
def main():
    print("Starting Microseismic EnKF Pipeline...")
    pipeline_start_time = time.time()

    # ==========================================
    # PHASE 1: SETUP & INITIALIZATION
    # ==========================================
    # 1. Build the 3D Grid
    x = np.arange(config.X_MIN, config.X_MAX + 1, config.DX, dtype=float)
    y = np.arange(config.Y_MIN, config.Y_MAX + 1, config.DX, dtype=float)
    z = np.arange(config.Z_MIN, config.Z_MAX + 1, config.DX, dtype=float)
    g = hypo.Grid3d(x, y, z)
   
    print(f"Grid built: {g.shape}")

    # 2. Load the True Velocity Model (for initial hypocenter guessing)
    print(f"Loading true velocity model from {config.TRUE_VELOCITY_FILE}...")
    if config.RUN_MODE == 'REAL':
        # Real data uses a CSV with specific columns
        Vesg = pd.read_csv(config.TRUE_VELOCITY_FILE, skiprows=0)
        Vtrue = np.array(Vesg)[:, 3]     
    else:
        # Synthetic data is a standard text array
        Vtrue = np.loadtxt(config.TRUE_VELOCITY_FILE)    
    

    # 3. Load the Initial Ensemble (E0)
    print(f"Loading initial ensemble from {config.INITIAL_E_FILE}...")
    E = np.loadtxt(config.INITIAL_E_FILE)
    
    # ---------------
    ## TEMPORAL SOLUTION TO FIT DIMENSIONS V*corrV

    E = E[:, :int(np.prod(g.shape))] ## Truncate ensemble parameters to match cell count
    # ---------------
    N, m = E.shape
    print(f"Ensemble loaded: {N} members, {m} parameters.")

    # Ensure output directory exists
    os.makedirs("EnKF_data", exist_ok=True)

    # ===================================================================================
    # PHASE 2: THE TIME LOOP (Hours of the day)
    # ===================================================================================
    # Example: Looping over hour 1 to 5 (adjust range as needed for your data)
    for j in range(12, 15):
        print(f"\n--- Processing Hour {j} ---")
        
        # Read the SEGY data for this hour using our core function
        data = enkf.read_seg(config.FOLDER_NAME, j)
        
        if len(data) == 0:
            print(f"No data found for hour {j}. Skipping.")
            continue
            
        print(f"Data found for hour {j}. Formatting events...")
        
        # Separate data into P and S wave picks per event
        evID = 0
        dataP = []
        dataS = []
        
        for i in range(data.shape[0]):        
            if i == 0:
                rcvID = 0
            elif data[i, 6] != data[i-1, 6]:
                evID += 1
                rcvID = 0 
                
            # If both P and S picks are present (not zero)
            if not math.isclose(data[i, 7], 0) and not math.isclose(data[i, 8], 0):
                # P-wave: indicator 0.0
                dataP.append([np.int64(evID), data[i, 7] * 1e-6, int(rcvID), 0., *data[i, 0:3]])
                # S-wave: indicator 1.0
                dataS.append([np.int64(evID), data[i, 8] * 1e-6, int(rcvID), 1., *data[i, 0:3]])
                rcvID += 1
                
        data_4inv = np.concatenate((np.array(dataP), np.array(dataS)))           
        nms = int(max(data_4inv[:, 0]) + 1) # Number of microseismic events
        print(f"Total events to process in hour {j}: {nms}")

        Eat = np.zeros(E.shape)
        Eft = np.zeros(E.shape)
        valid_events = 0

        # ==========================================
        # PHASE 3: THE EVENT LOOP (Data Assimilation)
        # ==========================================
        for k in range(nms): 
            print(f"  -> Assimilating Event {k+1}/{nms}...")
            
            try:
                # 1. Hypocenter Location
                ms, tt, rcv, nsta = enkf.inv_vcte(data_4inv, k, g) 
            except Exception as e: 
                print(f"     [!] Error locating event {k}: {e}")
                continue
                
            # 2. Check if the event is inside our 3D grid
            if (min(ms[:, 2]) > config.X_MIN and max(ms[:, 2]) < config.X_MAX and 
                min(ms[:, 3]) > config.Y_MIN and max(ms[:, 3]) < config.Y_MAX and 
                min(ms[:, 4]) > config.Z_MIN and max(ms[:, 4]) < config.Z_MAX):
                
                print(f"     Event inside grid at coordinates: {ms[:, 2:]}")
                
                try:
                    # 3. Forecast Step (Parallelized)
                    Ef, hE, rcv_kron, ms_kron = enkf.fore_step(E, ms, rcv, nsta, nev=1, g=g) 
                except Exception as e: 
                    print(f"     [!] Forecast failed: {e}")
                    continue
                
                # 4. Analysis Step
                R = np.eye(len(tt)) * config.STD_NOISE**2
                Ea = enkf.analysis_step(Ef, hE, R, tt) 
                
                # Accumulate results
                Eft += Ef
                Eat += Ea
                valid_events += 1
                
            else:
                print(f"     [!] Event {k} out of grid bounds. Skipping.")
                
        # ==========================================
        # PHASE 4: SAVING RESULTS
        # ==========================================
        if valid_events > 0:
            print(f"Saving averaged results for hour {j}...")
            # Average the ensembles over the valid events
            Eft_avg = Eft / valid_events
            Eat_avg = Eat / valid_events
            
            step_str = f"{j:02d}"
            
            # Save raw txt arrays
            np.savetxt(f"EnKF_data/Ea_{step_str}.txt", Eat_avg)
            
            # Save 3D visualizations using VTK (replaces old toXdmf)
            g.to_vtk({"Vp_forecasted": np.mean(Eft_avg, axis=0)}, f"EnKF_data/Vp_forecasted_{step_str}")
            g.to_vtk({"Vp_analysed": np.mean(Eat_avg, axis=0)}, f"EnKF_data/Vp_analysed_{step_str}")
            g.to_vtk({"Vpa-Vpf": np.abs(np.mean(Eat_avg, axis=0) - np.mean(Eft_avg, axis=0))}, f"EnKF_data/Vpa-Vpf_{step_str}")
            
            # Update E for the next hour's initial state
            E = np.copy(Eat_avg) 

    pipeline_end_time = time.time()
    total_time = pipeline_end_time - pipeline_start_time
    
    # Convert total seconds into hours, minutes, and seconds
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print("\n" + "="*50)
    print(" PIPELINE COMPLETE!")
    print(f" Total Execution Time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")
    print("="*50 + "\n")

    print("Pipeline Complete!")

# This is the trigger that actually runs the code when you type 'python main.py'
if __name__ == "__main__":
    main()
''' 

def main():
    print("Starting Microseismic EnKF Pipeline...")
    pipeline_start_time = time.time() 

    # ==========================================
    # PHASE 1: SETUP & INITIALIZATION
    # ==========================================
    # 1. Build the 3D Grid (Automatically uses the grid from the chosen RUN_MODE)
    x = np.arange(config.X_MIN, config.X_MAX + 1, config.DX, dtype=float)
    y = np.arange(config.Y_MIN, config.Y_MAX + 1, config.DX, dtype=float)
    z = np.arange(config.Z_MIN, config.Z_MAX + 1, config.DX, dtype=float)
    g = hypo.Grid3d(x, y, z)
    
    print(f"Grid built: {g.shape}")

    # 2. Load the True Velocity Model
    print(f"Loading true velocity model from {config.TRUE_VELOCITY_FILE}...")
    if config.RUN_MODE == 'REAL':
        # Real data uses a CSV with specific columns
        Vesg = pd.read_csv(config.TRUE_VELOCITY_FILE, skiprows=0)
        Vtrue = np.array(Vesg)[:, 3] 
    else:
        # Synthetic data is a standard text array
        Vtrue = np.loadtxt(config.TRUE_VELOCITY_FILE)


    # 3. Load the Initial Ensemble (E0)
    print(f"Loading initial ensemble from {config.INITIAL_E_FILE}...")
    E = np.loadtxt(config.INITIAL_E_FILE)
    
    # ----------------------------
    ## TEMPORAL SOLUTION TO FIT DIMENSIONS V*corrV

    E = E[:, :int(np.prod(g.shape))]
    # ----------------------------
    N, m = E.shape
    print(f"Ensemble loaded: {N} members, {m} parameters.")

    # Ensure output directory exists
    os.makedirs("EnKF_data", exist_ok=True)

    # ==========================================
    # PHASE 2 & 3: EXECUTION BRANCH
    # ==========================================
    
    if config.RUN_MODE == 'REAL':
        print("\n--- Running REAL Data Mode (SEGY) ---")
        # Example: Looping over hour 1 to 5 (adjust range as needed for your data)
        for j in range(5, 6):
            print(f"\n--- Processing Hour {j} ---")
            
            # Read the SEGY data for this hour using our core function
            data = enkf.read_seg(config.FOLDER_NAME, j)
            
            if len(data) == 0:
                print(f"No data found for hour {j}. Skipping.")
                continue
                
            print(f"Data found for hour {j}. Formatting events...")
            
            # Separate data into P and S wave picks per event
            evID = 0
            dataP = []
            dataS = []
            
            for i in range(data.shape[0]):        
                if i == 0:
                    rcvID = 0
                elif data[i, 6] != data[i-1, 6]:
                    evID += 1
                    rcvID = 0 

                # If both P and S picks are present (not zero)    
                if not math.isclose(data[i, 7], 0) and not math.isclose(data[i, 8], 0):
                    # P-wave: indicator 0.0
                    dataP.append([np.int64(evID), data[i, 7] * 1e-6, int(rcvID), 0., *data[i, 0:3]])
                    # S-wave: indicator 1.0
                    dataS.append([np.int64(evID), data[i, 8] * 1e-6, int(rcvID), 1., *data[i, 0:3]])
                    rcvID += 1
                    
            data_4inv = np.concatenate((np.array(dataP), np.array(dataS)))           
            nms = int(max(data_4inv[:, 0]) + 1) # Number of microseismic events
            print(f"Total events to process in hour {j}: {nms}")

            Eat = np.zeros(E.shape)
            Eft = np.zeros(E.shape)
            valid_events = 0

            # THE EVENT LOOP (Data Assimilation)
            for k in range(nms): 
                print(f"  -> Assimilating Event {k+1}/{nms}...")
                
                try:
                    # 1. Hypocenter Location
                    ms, tt, rcv, nsta = enkf.inv_vcte(data_4inv, k, g) 
                except Exception as e: 
                    print(f"     [!] Error locating event {k}: {e}")
                    continue
                    
                    # 2. Check if the event is inside our 3D grid
                if (min(ms[:, 2]) > config.X_MIN and max(ms[:, 2]) < config.X_MAX and 
                    min(ms[:, 3]) > config.Y_MIN and max(ms[:, 3]) < config.Y_MAX and 
                    min(ms[:, 4]) > config.Z_MIN and max(ms[:, 4]) < config.Z_MAX):
                    
                    print(f"     Event inside grid at coordinates: {ms[:, 2:]}")
                    
                    try:
                    # 3. Forecast Step (Parallelized)
                        Ef, hE, rcv_kron, ms_kron = enkf.fore_step(
                            E, ms, rcv, nsta, nev=1, g=g, 
                            pp_min=config.PP_MIN, pp_max=config.PP_MAX, 
                            rr_min=config.RR_MIN, rr_max=config.RR_MAX
                        ) 
                    except Exception as e: 
                        print(f"     [!] Forecast failed: {e}")
                        continue
                    
                    # 4. Analysis Step
                    R = np.eye(len(tt)) * config.STD_NOISE**2
                    Ea = enkf.analysis_step(Ef, hE, R, tt) 
                    
                    Eft += Ef
                    Eat += Ea
                    valid_events += 1
                    
                else:
                    print(f"     [!] Event {k} out of grid bounds. Skipping.")
                    
            # ==========================================
            # PHASE 4: SAVING RESULTS (REAL DATA)
            # ==========================================
            if valid_events > 0:
                print(f"Saving averaged results for hour {j}...")
                # Average the ensembles over the valid events
                Eft_avg = Eft / valid_events
                Eat_avg = Eat / valid_events
                
                step_str = f"{j:02d}"
                
                # Save raw txt arrays
                np.savetxt(f"EnKF_data/Ea_{step_str}.txt", Eat_avg)
                # Save 3D visualizations using VTK (replaces old toXdmf)
                g.to_vtk({"Vp_forecasted": np.mean(Eft_avg, axis=0)}, f"EnKF_data/VpA_forecasted_{step_str}")
                g.to_vtk({"Vp_analysed": np.mean(Eat_avg, axis=0)}, f"EnKF_data/VpA_analysed_{step_str}")
                g.to_vtk({"Vpa-Vpf": np.abs(np.mean(Eat_avg, axis=0) - np.mean(Eft_avg, axis=0))}, f"EnKF_data/VpAa-Vpf_{step_str}")
                
                # Update E for the next hour's initial state
                E = np.copy(Eat_avg) 

    elif config.RUN_MODE == 'SYNTHETIC':
        print(f"\n--- Running SYNTHETIC Benchmark Mode (Step {config.SYNTH_STEP}) ---")
        
        step_str = str(config.SYNTH_STEP)
        
        # 1. Load Synthetic Text Files
        print(f"Loading synthetic data from {config.SYNTH_DIR}...")
        rcv = np.loadtxt(os.path.join(config.SYNTH_DIR, "rcvRND.txt"), skiprows=1)
        tt = np.loadtxt(os.path.join(config.SYNTH_DIR, f"tt_synth{step_str}.txt"))
        ms_raw = np.loadtxt(os.path.join(config.SYNTH_DIR, f"hinit{step_str}.txt"))
        hyp_raw = np.loadtxt(os.path.join(config.SYNTH_DIR, f"h_true{step_str}.txt")) # Load true location

        nsta = rcv.shape[0]
        
        # 2. Format 'ms' to match modern enkf_core expectations: [evID, t0, X, Y, Z]
        # Ensure it is a 2D array even if it's just one event
        if ms_raw.ndim == 1: ms_raw = ms_raw.reshape(1, -1)
        if hyp_raw.ndim == 1: hyp_raw = hyp_raw.reshape(1, -1)
        
        ms = np.zeros((1, 5))
        ms[0, 2:5] = ms_raw[0, -3:] # Map the X, Y, Z coordinates to columns 2, 3, 4
        
        ms_true = np.zeros((1, 5))
        ms_true[0, 2:5] = hyp_raw[0, -3:]

        # 3. Data Assimilation
        print("  -> Assimilating Synthetic Event...")
        try:
            Ef, hE, rcv_kron, ms_kron = enkf.fore_step(
                E, ms, rcv, nsta, nev=1, g=g, 
                pp_min=config.PP_MIN, pp_max=config.PP_MAX, 
                rr_min=config.RR_MIN, rr_max=config.RR_MAX,
                ms_true=ms_true
            ) 
            
            R = np.eye(len(tt)) * config.STD_NOISE**2
            Ea = enkf.analysis_step(Ef, hE, R, tt)
            
            # ==========================================
            # PHASE 4: SAVING RESULTS (SYNTHETIC DATA)
            # ==========================================
            print("Saving synthetic results...")
            np.savetxt(f"EnKF_data/SYNTH_Ea_{step_str}.txt", Ea)
            g.to_vtk({"Vp_forecasted": np.mean(Ef, axis=0)}, f"EnKF_data/SYNTH_Vp_forecasted_{step_str}")
            g.to_vtk({"Vp_analysed": np.mean(Ea, axis=0)}, f"EnKF_data/SYNTH_Vp_analysed_{step_str}")
            g.to_vtk({"Vpa-Vpf": np.abs(np.mean(Ea, axis=0) - np.mean(Ef, axis=0))}, f"EnKF_data/SYNTH_Vpa-Vpf_{step_str}")
            
        except Exception as e:
            print(f"     [!] Synthetic Assimilation failed: {e}")

    # ==========================================
    # PHASE 5: PIPELINE TIMING
    # ==========================================
    pipeline_end_time = time.time()
    total_time = pipeline_end_time - pipeline_start_time
    
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print("\n" + "="*50)
    print(" PIPELINE COMPLETE!")
    print(f" Total Execution Time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")
    print("="*50 + "\n")

    print("Pipeline Complete!")

# This is the trigger that actually runs the code when you type 'python main.py'
if __name__ == "__main__":
    main()
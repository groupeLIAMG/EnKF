import os
import numpy as np

# ==========================================
# RAYTRACING ENGINE CONFIGURATION
# ==========================================
# Choose the raytracing algorithm:'SPM'or 'FSM'

RAYTRACE_METHOD = 'FSM'
# Toggle GPU acceleration for FSM (Requires OpenCL support)
USE_GPU = True

# Force PyOpenCL to use Platform 1 (The NVIDIA RTX 5070) and silence warnings
if USE_GPU:
    os.environ["PYOPENCL_CTX"] = "0"
    os.environ["PYOPENCL_COMPILER_OUTPUT"] = "0"

# Enforce Float32 precision
import numpy as np
PRECISION = np.float32 

# ==========================================

# Set to 'REAL' for real SEGY data, or 'SYNTHETIC' to run the benchmark

RUN_MODE = 'REAL'

if RUN_MODE == 'REAL':
    # --- Real Grid ---
    NX, NY, NZ = 121, 121, 88     # Nodes
    DX, DY, DZ = 15.0, 15.0, 15.0 # 15m spacing
    
    X_MIN, X_MAX = 5750.0, 7550.0 
    Y_MIN, Y_MAX = 1900.0, 3700.0
    Z_MIN, Z_MAX = 1450.0, 2755.0

    # --- Real Perturbation Parameters ---
    PP_MIN, PP_MAX = 0.05, 0.05  # Flat 5% perturbation
    RR_MIN, RR_MAX = 50, 100     # Large radius for real data

    # --- Real Paths ---
    FOLDER_NAME = r'C:\Users\marti\FRQ_Project\example\SEGY_FilesForENKF_Picks\2025\02\10'
    TRUE_VELOCITY_FILE = r'C:\Users\marti\FRQ_Project\example\Laronde_3DVM_forINRS\laronde_3dvm_6contacts.vtr'
    INITIAL_E_FILE = r'C:\Users\marti\FRQ_Project\example\Laronde_3DVM_forINRS\larondeVp_ensemble.npy'
    #INITIAL_E_FILE = r'C:\Users\marti\FRQ_Project\EnKF_data\Ea_04.txt'

elif RUN_MODE == 'SYNTHETIC':
    # --- Synthetic Grid (From Dip 2018 script) ---
    X_MIN = 9600
    X_MAX = 10440
    Y_MIN = 10036
    Y_MAX = 10416
    Z_MIN = 2626
    Z_MAX = 3206
    DX = 20.0
    
    # --- Synthetic Perturbation Parameters ---
    PP_MIN, PP_MAX = -0.1, 0.1   # -10% to +10% perturbation
    RR_MIN, RR_MAX = 0, 60       # Small radius for synthetic benchmark

    # --- Synthetic Paths ---
    SYNTH_STEP = 3
    SYNTH_DIR = r'C:\Users\marti\FRQ_Project\example\synthetic'
    
    # 'f-strings' here to automatically update the file names based on the step
    TRUE_VELOCITY_FILE = rf'{SYNTH_DIR}\Vtrue{SYNTH_STEP-1}.txt'
    INITIAL_E_FILE = rf'{SYNTH_DIR}\Ea{SYNTH_STEP-1}.txt'

# --- Algorithm Parameters ---
STD_NOISE = 1.e-3
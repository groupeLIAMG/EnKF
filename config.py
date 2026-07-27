# ==========================================
# Set to 'REAL' for real SEGY data, or 'SYNTHETIC' to run the benchmark
RUN_MODE = 'SYNTHETIC'

if RUN_MODE == 'REAL':
    # --- Real Grid ---
    X_MIN = 9600  
    X_MAX = 10700
    Y_MIN = 9700
    Y_MAX = 10500
    Z_MIN = 2425
    Z_MAX = 3275
    DX = 10.0
    
    # --- Real Perturbation Parameters ---
    PP_MIN, PP_MAX = 0.05, 0.05  # Flat 5% perturbation
    RR_MIN, RR_MAX = 50, 100     # Large radius for real data

    # --- Real Paths ---
    FOLDER_NAME = r'C:\Users\marti\FRQ_Project\example\SEGY_FilesForINRS\SEGY_FilesForINRS\2025\02\17'
    TRUE_VELOCITY_FILE = r'C:\Users\marti\FRQ_Project\example\Glencore_ESG_data\vp_esg.csv'
    INITIAL_E_FILE = r'C:\Users\marti\FRQ_Project\example\Glencore_ESG_data\Ea-1.txt'
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
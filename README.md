# Microseismic Ensemble Kalman Filter (EnKF) Pipeline

A modernized, highly-parallelized Python implementation of the Ensemble Kalman Filter (EnKF) for real-time updating of 3D velocity models using microseismic travel times.

## Repository Structure

`config_enkf.py`: The master configuration file. Contains a dynamic switch `RUN_MODE = 'REAL'` or `RUN_MODE = 'SYNTHETIC'`  to automatically toggle between real mining datasets and benchmark testing.

`enkf_core.py`: Contains the parallelized `fore_step` (Forecast) utilizing `joblib` and the `analysis_step` (Kalman Analysis).

`functions.py`: Contains supporting mathematical and visualization functions (e.g., spherical perturbation algorithms).

`main.py`: The primary execution script. Handles the time-loop (hour-by-hour processing), manages data assimilation flow, and saves `.vtk` and `.txt` outputs.

## Requirements and Installation
This project requires specific custom libraries:

**HypoPy**: Used for 3D grid generation, spatial management, and handling the geometric relationships between seismic sources and receivers.

**ttcrpy**: A high-performance Python/C++ library utilized for computing theoretical travel times through the 3D velocity grid.

Development is made with python version 3.12.13


Install the standard libraries (`ttcrpy` and `hypopy`):

`pip install numpy pandas obspy tqdm joblib ttcrpy`

`pip install git+https://github.com/groupeLIAMG/hypopy.git`

If you have VTK compiled with python on your system, it is possible to save velocity models and raypaths for posterior visualization (e.g. in paraview).

## Usage
Open `config_enkf.py` and set your desired `RUN_MODE` and grid parameters. Ensure paths to your local data directories are correct.

Run the main pipeline:

`python main.py`

References:

@inproceedings{dip2019microseismic,
  title={Microseismic Monitoring of Mines in Real Time with Ensemble Kalman Filter},
  author={Dip, AC and Giroux, B and Gloaguen, E},
  booktitle={81st EAGE Conference and Exhibition 2019},
  year={2019}
}

@Book{evensen2009data,
  title     = {Data assimilation: the ensemble Kalman filter},
  publisher = {Springer Science \& Business Media},
  year      = {2009},
  author    = {Evensen, Geir},
}

@InProceedings{raanes2016intro2,
  author       = {Raanes, Patrick N.},
  title        = {Introduction to Data Assimilation and the Ensemble {K}alman Filter. {S}econd edition.},
  year         = {2016},
  month        = {October},
  organization = {Nansen Environmental and Remote Sensing Center},
}

# EnKF

The Ensemble Kalman filter applied to microseismic data for updating 3D velocity models

! STILL UNDER TEST IN SYNTHETICS !

Requirements

Development is made with python version 3.6

You need to add this to you PYTHONPATH:

https://github.com/groupeLIAMG/hypopy

Note: You need to compile the python wrapper for the C++ raytracing code in https://github.com/groupeLIAMG/ttcr and add it to your PYTHONPATH to be able to run hypo.py

If you have VTK compiled with python on your system, it is possible to save velocity models and raypaths for posterior visualization (e.g. in paraview).

References:
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

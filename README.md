# LaCe Analysis - Jonathan Morrell

This repository contains data and analysis code for the La(p,x) cross-section measurements performed in Sep 2017 at LBNL.  Currently data and code are to be considered as preliminary.

Code has been tested with `python 2.7` on Ubuntu 17.1.

Individual plots with peak fits can be viewed by uncommenting the lines under the __main__ statement (line 255) in code/analysis.py.

## How to run analysis code


To run the analysis:

 1. Change the `_empire_path`, `_mcnp_path`, and `_talys_path` variables in analysis.py and nudat_manager.py to the correct paths to those executables on your local machine.

 2. `make data` - Computes TALYS and EMPIRE cross-sections and scrapes decay data from NNDC.  If EMPIRE doesn't finish successfully re-run `make data` until all energy points have been calculated.

 3. `make analysis` - Runs the analysis code, should take about 45 minutes in total to perform peak fits, calibration, minimize the proton current uncertainty and calculate cross sections.

# ********************************************************************************************************************
# ********************************************************************************************************************
#
# Exoplanetary Transit Reduction
# Author: Walter Martins-Filho
# Contact: walter at on.br
#          waltersmartinsf at gmail.com 
# Created: July, 18, 2016
# Version 1.5
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# Goal: Reduction and photometry of a data set of exoplanetary transits
#
# Code Files:
# - ExoTRed.sh 									Bash shell script to run this code and README instructions
# - ExoTRed.py 									Python script with the step-by-step data reduction and 
#												photometry analysis
# - ExoTRed_input.yaml 							Input file with information of objetcs, observatory, ...
# - login.cl 									IRAF default login.cl file to work with this code.
# - ./pyraf/									PyRAF folder contains PyRAF cache file.
# - ./example/									Example folder. This folder only exist if you ask to author for the
#												big version of ExoTRed.tar.gz file.
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# INPUT files:
#
# data_path/exoplanet*.fits 					Exoplanet science images in this format:
# 												name_of_exoplanet*.fits where the name of the exoplanet 
#												is set on the ExoTRed_input.yaml file
# data_path/bias*.fits 							Bias images files
# data_path/flat*.fits 							Flat images files

# **Notes**: 
#			(1) It is important to make separate folders for separate filters whatever if the nomenclature of the 
#				filter is present in the name of the images.
#			(2) When complete the YAML file, please include all path in data_path and save_path variables.
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# OUTPUT Tree Files:
# 
# save_path/
# save_path/background/         				Folder with *.pik (binaries in python) with sky backgound information.
# save_path/background/bkg_data_*_.pik 			Files with sky background information for subtracted during photometry.
# save_path/phot_results/						Folder contains the results of aperture photometry.
# save_path/phot_results/*_flux_radius_*.csv	*csv files with photometry information in counts. 
#												Set columns in this order:
#												hoststar, err_hostar, refstar1, err_refstar1, refstar2, err_refstar2, ...
# save_path/superbias.fits						Masterbias image in FITS format.
# save_path/Bflat*.fits 						Flat images with subtracted bias.
# save_path/superflat.fits 						Masterflat image in FITS format.
# save_path/ABflat*.fits 						Flat images with subtracted bias and normalized.
# save_path/B*exoplanet*.fits 					Science images with bias subtracted.
# save_path/AB*exoplanet*.fits 					Science images with bias subtracted and flatfielding.
# save_path/results.csv							Table with time information for each science image from Python Routines:
#												image, UTC, JD, ST, ST_isot, RA_SUN, DEC_SUN, D_SUN, HJD
#												where
#												image = image name
#												UTC = Coordinated Universal Time
#												JD = Julian Date
#												ST = Sideral Time
#												ST_isot = Sideral Time in ISOT format
#												RA_SUN = Right ascension of the Sun
#												DEC_SUN = Declination of the Sun
#												D_SUN: Earth distance from the Sun
#												HJD: Heliocentric Julian Date
# save_path/results_iraf_calibrations.csv 		Table with time information for each science image from IRAF Routines
#												and airmass. Set in this order:
#												HJD, JD, ST, Airmass
# save_path/rawflux_StarColumn_*.png            Plots using data from aperture photometry results table with minimum
#												scatter in the flux. 
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# How to run this script:
# 
# To run this script, include it and the input YAML file on the directory that You want to save and type on terminal:
#
# bash ExoTRed.sh
# 
# or, if ExoTRed.sh had executable permission, type:
# 
# ./ExoTRed.sh
#
# ********************************************************************************************************************
# ********************************************************************************************************************
# 
# HOW TO INSTALL DEPENDENCIES:
#
# - Anaconda Python: https://www.continuum.io/downloads
# - Astroconda Environment with Python 2.7 and IRAF legacy : http://astroconda.readthedocs.io/en/latest/installation.html
#
# How to Install Dependencies:
#
# (1) Download the bash script from Anaconda site for your OS system. For OSX (Macbooks) or for Windows there is a 
#     graphical installer. If you choose the bash script, save in a the directory that you want to install Anaconda 
#     Python Environment. The default path is /home/user/anaconda2 or /home/user/anaconda3 (if you choose Python 3.5).
#     After that, You will type on terminal:
#
#     bash Anaconda2-4.1.1-Linux-x86_64.sh
#     
#     or, if you choose Python 3.5:
#
#     bash Anaconda3-4.1.1-Linux-x86_64.sh 
#
#     The script will ask if yopu accpet the terms and if you accpet the path. Following it instructions.
#
# (2) Open a new terminal, this is done for restart your source in .bashrc (or similar). Next, You will need to install 
#     Astroconda Environment with Python 2.7 and IRAF legacy. PyRAF and other packages that works with IRAF only support
#     Python 2.7 You don't need to install Anaconda 2.7 on the STEP 1 because here it will download the python packages 
#     that IRAF needed. On terminal, type:
#     
#     conda config --add channels http://ssb.stsci.edu/astroconda
#     
#     conda create -n astroconda python=2.7 iraf pyraf stsci
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# There is a copy of this in the shell file.
#
# ********************************************************************************************************************
# ********************************************************************************************************************
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
# More in ExoTRed.sh file.
#
# ********************************************************************************************************************
# ********************************************************************************************************************
#
# Goal: Reduction and photometry of a data set of exoplanetary transits
#
# ********************************************************************************************************************
# ********************************************************************************************************************

import sys
import os
import numpy as np
from pandas import DataFrame, read_csv
import glob
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt #plot library
def init_plotting():
    """
    Customize graphical output.
    """
    plt.rcParams['figure.figsize'] = (14.0,8.0)
    plt.rcParams['font.size'] = 20
    #plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 0.65*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1

init_plotting()

#Includinfg sources folder as the main folder for core scripts
sys.path.insert(0, './sources')
# ********************************************************************************************************************
# ********************************************************************************************************************
#import the core script ExoTRed as a internal package
import ExoTRed_core as exotred
import usefull_functions as use

# ********************************************************************************************************************
# ********************************************************************************************************************
# Start to run:
start_time = exotred.time.time() #Start time of this process from the clock of OS.
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Obtain the information in the YAML file: \n'
data_path, save_path, input_file = exotred.input_info('ExoTRed_input.yaml') #include the YAML information
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Create the superbias image. \n'
exotred.masterbias(input_file) #create a superbias image
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Create the superflat image \n'
exotred.masterflat(input_file) #create a superflat image
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Reduce and calibrate science images \n'
exotred.science_reduction(input_file) #reduce and calibrate in flat e bias the science images
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Obtain time information of each science image \n'
exotred.time_info(input_file) #obtain the time information for each science image (based on astropy)
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Obtain airmass and put the time information in the header of each science image'
exotred.time_calibration(input_file) #obtain the airmass, and put all information of time in the header (based on pyraf).
# ********************************************************************************************************************
# ********************************************************************************************************************
time_final = (abs(exotred.time.time()-start_time)/60.) #time in minutes
print 'Duration for reduction = ',time_final,' minutes \n'
# ********************************************************************************************************************
# ********************************************************************************************************************
# print 'Obtain the sky background information for each image: \n'
#exotred.bkg_info(input_file)
# ********************************************************************************************************************
# ********************************************************************************************************************
# print 'Loading sky background data: \n'
# bkg_data, bkg_rms = exotred.bkg_info(input_file)
# ********************************************************************************************************************
# ********************************************************************************************************************
print 'Making the aperture photometry: \n'
exotred.phot_aperture(input_file)
# ********************************************************************************************************************
# ********************************************************************************************************************

original_path = os.getcwd()
os.chdir(input_file['save_path'])
print 'Photometry Analysis \n'
print 'Reading '+input_file['exoplanet']+'*.csv files ....\n'
files_csv = np.sort(glob.glob(input_file['save_path']+'/phot_results/'+input_file['exoplanet']+'*.csv'))
scatter = np.zeros(len(files_csv))
for i in range(len(files_csv)):
    phot_data = read_csv(files_csv[i])
    scatter[i] = np.std(phot_data[str(0)])
    use.update_progress((i+1.)/len(files_csv))
hjd = read_csv(input_file['save_path']+'/results_iraf_calibrations.csv')
airmass = hjd.Airmass
hjd = hjd.HJD.values
print '... done!'
id_min = scatter.argmin() #index of the min scatter file
id_max = scatter.argmax() #index for the maximum scatter file
print 'The smallest scatter is: '+str(files_csv[id_min])
print 'Which is file: '+files_csv[id_min]
print('Working @'+files_csv[id_min]+' that is the min scatter')
print('... Read '+files_csv[id_min]+' ...')
data_min_scatter = read_csv(files_csv[id_min])
print('... done.')
#Setting the table with minimum scatter
rawdata, err_rawdata = data_min_scatter[str(0)].values,data_min_scatter[str(1)].values
# ********************************************************************************************************************
# ********************************************************************************************************************

print 'Make a beautifull plot from flux results: \n'
for i in range(len(data_min_scatter.columns)-2):
    if (i % 2 == 0) and (i != int(data_min_scatter.columns[-1])-2):
        refstar = data_min_scatter[str(i+2)].values
        err_refstar = data_min_scatter[str(i+3)].values
        rawflux = rawdata/refstar
        eflux = rawflux*np.sqrt((err_rawdata/rawdata)**2 + (err_refstar/refstar)**2)
        f = plt.figure()
        plt.suptitle('Minimum scatter of the flux \n @'+files_csv[id_min]+'\ncolumns = '+str(i+3)+' Reference Star/ '+str(i+4)+' Err_flux',fontsize=12)
        gs1 = GridSpec(2, 2, width_ratios=[1,2],height_ratios=[4,1])
        gs1.update(wspace=0.5)
        ax1 = plt.subplot(gs1[:-1, :])
        ax2 = plt.subplot(gs1[-1, :])
        ax1.grid()
        ax1.errorbar(hjd,rawflux,yerr=eflux,ecolor='g')
        ax1.set_xticklabels([])
        ax1.set_ylabel('Relative Flux')
        ax2.grid()
        ax2.plot(hjd,airmass,color='green')
        plt.yticks(np.array([airmass.min(), (airmass.min()+airmass.max())/2., airmass.max()]))
        ax2.set_xlabel('JD')
        ax2.set_ylabel('airmass')
        plt.savefig(save_path+'/rawflux_StarColumn_'+str(i+3)+'.png')

os.chdir(original_path)
print 'Total duration = ',abs(exotred.time.time()-start_time)/60.,' minutes'
plt.show()
# ********************************************************************************************************************
# ********************************************************************************************************************

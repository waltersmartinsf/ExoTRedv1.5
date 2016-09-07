"""
This file contains all the functions for ExoTRed Package

@author: Walter Martins-Filho
e-mail: walter at on.br
        waltersmartinsf at gmail.com

"""


"""
Loading packges and python useful scritps
"""
from pyraf import iraf #loading iraf package
from login import * #loading login.cl parameters for iraf
from ExoSetupTaskParameters import * #loading setup from PyExoDRPL
import usefull_functions as use
import glob #package for list files
import os #package for control bash commands
import yaml #input data without any trouble
import string
import numpy as np
from numpy import arange,array
from astropy.io import fits
from astropy.time import Time #control time in fits images
from astropy.coordinates import SkyCoord,get_sun,ICRS #get the Sun position from a instant of time
from photutils import CircularAperture, aperture_photometry,CircularAnnulus,Background #Photometry routines 
from astropy.table import hstack
from astropy.io.misc import fnpickle, fnunpickle #create binary files
from pandas import DataFrame, read_csv
from string import split #use to unconcanated a string in parts
import time
#import string #transform a list in a string of caracters
#from pandas import HDFStore #save data in a database


def input_info(path_input):
    """
    Obtain information about the files that will be reduced or analised.
    ___
    Input:

    path_input: string, path of your yaml file directory with the information of the data

    Return:

    save_path: string, path where is your data to be reduced.
    data_path: string, path where the data reduced will be saved.
    input_file: dictionary, information from YAML file.

    """
    #path for your data directory, path for your data save, and names for the lists
    #Import with yaml file: input path and prefix information for files
    input_file = glob.glob('ExoTRed_input.yaml')
    if input_file: #if exist the input file, the code will obtain all information We need to run all tasks
        if len(input_file) == 1: #if there is only one yaml file, obtain data and save paths, and return that with a dictionary with information
            print 'reading input file ... \n'
            file = yaml.load(open(input_file[0])) #creating our dictionary of input variables
            data_path = file['data_path']
            save_path = file['save_path']
            print '....  done! \n'
            if len(input_file) > 1: #if are more than one yaml file,. the code will ask to you remove the others.
                print 'reading input file ... \n'
                print '.... there is more than 1 input_path*.yaml.\n \nPlease, remove the others files that you do not need. \n'
                raise SystemExit
    else:
        #if aren't a yaml file, the code ask for you to put a valid yaml file path.
        print 'There is no input_path*.yaml. \nPlease, create a input file describe in INPUT_PARAMETERS.'
        raise SystemExit
    input_file = file #creating a better name to our dictionary info
    return data_path, save_path, input_file

def masterbias(input_file):
    """
    Obtain the masterbias.fits image.
    ___
    Input:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    Output:
    It is possible that the function return some of these values:

    0. Create the masterbias image on the save_path.
    1. It do not create the masterbias image, because of some error
    ___
    """
    #Set original directory
    original_path = os.getcwd()
    save_path = input_file['save_path']
    data_path = input_file['data_path']
    #Change your directory to data diretory
    os.chdir(data_path)
    #list all bias images
    bias = glob.glob('bias*.fits')
    print 'Loading bias images \nTotal of bias files = ',len(bias),'\nFiles = \n'
    print bias
    print '\nCreating superbias \n'
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #copy bias images to save_path
    os.system('cp bias*.fits '+save_path)
    #change to sabe_path
    os.chdir(save_path)
    #verify if previous superbias exist
    if os.path.isfile('superbias.fits') == True:
        os.system('rm superbias.fits')
    #create the list of bias images
    bias_list = string.join(bias,',')
    #combine the bias image and create the superbias
    iraf.imcombine(bias_list,'superbias.fits')
    iraf.imstat('superbias.fits')
    #clean previos bias files
    print '\n Cleaning bias*.fits images ....\n'
    os.system('rm bias*.fits')
    print '\n.... done.'
    #print output
    #test of outpu value
    #os.remove('superbias.fits')
    #Verify if the image was created:
    output = glob.glob('superbias*.fits')
    if len(output) != 0:
        output = 0
    else:
        output = 1
    #Return to original directory
    os.chdir(original_path)
    #END of the masterbias reduction messsage
    print '\nsuperbias.fits created!\n'
    print '\nEND of superbias reduction!\n'
    #obtain the value of return
    if output == 1:
        print '!!! ERROR/WARNING !!!'
        print 'Check if the superbias was created or if there is more than one superbias image.'
    return output

def masterflat(input_file):
    """
    Obtain the masterflat image for calibration.
    ___
    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    OUTPUT:
    It is possible that the function return some of these values:

    0. Create the masterflat image on the save_path.
    1. It do not create the masterflat image, because of some erros.
    """
    #set original directory
    original_path = os.getcwd()
    data_path = input_file['data_path']
    save_path = input_file['save_path']
    #Change your directory to data diretory
    os.chdir(data_path)
    #list all flat images
    flat = glob.glob('flat*.fits')
    print 'Loading flat images \nTotal of flat files = ',len(flat),'\nFiles = \n'
    print flat
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #create a list of bias images and copy images to save_path
    os.system('cp flat*.fits '+save_path)
    #creating the names of flat with bias subctracted
    bflat = []
    for i in flat:
        bflat.append('B'+i)
    print '\n Names os flat images with bias subtracted: \n \n',bflat
    #change for save_path directory
    os.chdir(save_path)
    #verify if previous superbias exist
    if os.path.isfile('superflat.fits') == True:
        os.system('rm superflat.fits')
    #verify if exits previous bflat*.fits files and remove then.
    for i in bflat:
        if os.path.isfile(i) == True:
            os.system('rm -f '+i)
    print '\nCreating superflat .... \n'
    #create the list of flat images  and bflat images
    #flat = string.join(flat,',')
    #bflat = string.join(bflat,',')
    print '\n Subtracting bias from flat images and creating bflat images.... \n'
    #iraf.imarith()
    for i in range(len(flat)):
        iraf.imarith(flat[i],'-','superbias.fits',bflat[i])
        #print statistics from bflat*.fits images
        iraf.imstat(bflat[i])
    print '\n .... done \n'
    #clean previos flat*.fits files
    print '\n Clean flat*.fits images .... \n'
    os.system('rm flat*.fits')
    print '\n .... done. \n'
    #normalizing each flat
    print '\nNormalizing each flat ....\n'
    #checking if mean from numpy is the same from your bflat images using imstat
    #take the mean of each bflat image
    bflat_mean = np.zeros(len(bflat))
    for i in range(len(bflat)):
        image = fits.getdata(bflat[i])
        image = np.array(image,dtype='Float64')
        bflat_mean[i] = round(np.mean(image))
    image = 0 #clean image allocate to this variable
    print 'The mean of each bflat image, respectivaly ...'
    print bflat_mean
    #creating the names of bflat images after the normalization:
    abflat = []
    for i in bflat:
        abflat.append('A'+i)
    print '\n Names os bflat images with bias subtracted and normalizad: \n \n',abflat
    #verify if exist previous ABflat*.fits images and remove then.
    for i in abflat:
        if os.path.isfile(i) == True:
            os.system('rm -f '+i)
    for i in range(len(abflat)):
        iraf.imarith(bflat[i],'/',bflat_mean[i],abflat[i])
    print '\n.... done!\n'
    # print '\n Cleaning bflat*.fits images ....\n'
    # os.system('rm Bflat*.fits')
    print '\n.... done.\n'
    print 'Statistics of the abflat*.fits images .... \n'
    for i in range(len(abflat)):
        iraf.imstat(abflat[i])
    print '\n Combining abflat images ....\n'

    # ablist = string.join(abflat,',')
    # iraf.imcombine(ablist,'superflat.fits')
    #change how import flat files
    #usning the abflat list of flat files We will create a pandas python dataframe
    ablist = DataFrame(ablist)
    ablist.columns=['flat_files']
    ablist.to_csv('flat_list',index_label=False,index=False,header=False)
    #combine all flat images
    iraf.imcombine('@flat_list','superflat.fits')
    iraf.imstat('superflat.fits')
    print '\n .... done. \n'
    # print '\nCleaning ABflat*.fits images ....\n'
    # os.system('rm ABflat*.fits')
    print '\n.... done!'
    #Verify if the image was created:
    output = glob.glob('superflat*.fits')
    if len(output) != 0:
        output = 0
    else:
        output = 1
    #Return to original directory
    os.chdir(original_path)
    #last mensage
    print '\n MASTERFLAT.FITS created! \n'
    print '\n END of Data Reduction for create a masterflat.fits file. \n'
    #obtain the value of return
    if output == 1:
        print '!!! ERROR/WARNING !!!'
        print 'Check if the superbias was created or if there is more than one superbias image.'
    return output

def science_reduction(input_file):
    """
    Calibrate science images with masterflat (or superflat) and masterbias (or superbias) images.
    ___
    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    OUTPUT:
    It is possible that the function return some of these values:

    0. Create the masterflat image on the save_path.
    1. It do not create the masterflat image, because of some erros.
    """
    #name of the planet
    planet = input_file['exoplanet']
    #set original directory
    original_path = os.getcwd()
    save_path = input_file['save_path']
    data_path = input_file['data_path']
    #Change your directory to data diretory
    os.chdir(data_path)
    #list all flat images
    exoplanet = glob.glob(planet+'*.fits')
    print '\nLoading exoplanet images \nTotal of '+planet+'*.fits  files = ',len(exoplanet),'\nFiles = \n'
    print exoplanet
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #create a list of bias images and copy images to save_path
    print '\nCopy science images to save_path directory to main reduction: ....'
    os.system('cp '+planet+'*.fits '+save_path)
    print '\n .... done. \n'
    #change to save_path
    os.chdir(save_path)
    #create the names for exoplanet science mages with bias subtracted
    bexoplanet = []
    for i in exoplanet:
        bexoplanet.append('B'+i)
        #verify if previous superbias exist
        if os.path.isfile('B'+i) == True:
            os.system('rm B'+i)
    print '\n Will be create this images: \n'
    print bexoplanet
    #exoplanet = string.join(exoplanet,',') #create the list string of exoplanet science images
    #bexoplanet = string.join(bexoplanet,',')#create the list string of bexoplanet science images
    print '\nSubtracting superbias.fits from all '+planet+'*.fits images ....\n'
    for i in range(len(exoplanet)):
        iraf.imarith(exoplanet[i],'-','superbias.fits',bexoplanet[i])
        use.update_progress((i+1.)/len(bexoplanet))
    print '\n.... cleaning '+planet+'*.fits images\n'
    os.system('rm '+planet+'*.fits')
    print '\n Statistics of B'+planet+'*.fits images: \n'
    for i in range(len(bexoplanet)):
        iraf.imstat(bexoplanet[i])
    print '\nFlatfielding the B'+planet+'*.fits ....\n'
    #create the names for exoplanet science images with bias subtracted and flatfielding
    abexoplanet = []
    for i in bexoplanet:
        abexoplanet.append('A'+i)
        #verify if previous superbias exist
        if os.path.isfile('A'+i) == True:
            os.system('rm A'+i)
    print '\n Will be create this images: \n'
    print abexoplanet
    #flatifielding images
    for i in range(len(abexoplanet)):
        iraf.imarith(bexoplanet[i],'/','superflat.fits',abexoplanet[i])
        use.update_progress((i+1.)/len(abexoplanet))
    # print '\n.... cleaning B'+planet+'*.fits images\n'
    # os.system('rm B'+planet+'*.fits')
    print '\n Statistics of AB'+planet+'*.fits images: \n'
    for i in range(len(abexoplanet)):
        iraf.imstat(abexoplanet[i])
    os.chdir(original_path) #change to save_path
    return

def time_info(input_file):
    """
    Obtain the Sideral Time and the Heliocentric Jullian Date from the header of the images.
    ___
    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    OUTPUT:
    It is possible that the function return some of these values:

    0. Create the masterflat image on the save_path.
    1. It do not create the masterflat image, because of some erros.
    """
    original_path = os.getcwd() #set original directory
    save_path = input_file['save_path']
    planet = input_file['exoplanet'] #set exoplanet name
    print '\nObtain the images .... \n'
    print 'Change to ', save_path
    os.chdir(save_path) #change to save directory where is our scvience images
    images = sorted(glob.glob('AB'+input_file['exoplanet']+'*.fits'))
    print '\nImages = \n',images
    tempo_loc = [] #time object
    SUN = [] #Sun coordinate object
    ra_sun, dec_sun, dsun = np.zeros(len(images)),np.zeros(len(images)),np.zeros(len(images)) #sun coordinates
    JD = np.zeros(len(images)) #julian date from time object
    ST = np.zeros(len(images))
    HJD = np.zeros(len(images))
    #create the exoplanet object coordianate
    exoplanet = SkyCoord(dec=input_file['DEC'],ra=input_file['RA'],unit=('deg','deg'),frame=input_file['frame'])
    print '\nObtain data info from header ....\n'
    for i in range(len(images)):
        hdr = fits.getheader(images[i])
        UTC = hdr['date-obs']+'T'+hdr['UT'] #string that contain the time in UTC in isot format
        tempo_loc.append(Time(UTC,scale=input_file['scale-time'],format='isot',location=(input_file['lon-obs'],input_file['lat-obs'])))#,input_data['altitude'])))
        JD[i] = tempo_loc[i].jd
        ST[i] = tempo_loc[i].sidereal_time('apparent').hour
        SUN.append(get_sun(tempo_loc[i]))
        ra_sun[i],dec_sun[i] = SUN[i].ra.deg, SUN[i].dec.deg
        dsun[i] = SUN[i].distance.value
        HJD[i] = use.hjd_date(JD[i],dsun[i],dec_sun[i],ra_sun[i],exoplanet.dec.deg,exoplanet.ra.deg,circular_orbit=input_file['circular_orbit'])
        use.update_progress((i+1.)/len(images))
    print '\n.... done.\n'
    print '\n Time from header = \n'
    #print '\nImages ** UTC (YYYY-MM-DDTHH:MM:SS) ** JD (7d.5d) ** ST (hours) ** ST (HH:MM:SS) ** Sun Coordinate (epoch,RA,DEC,Distance) (deg,deg,AU) \n'
    ST_string = []
    for i in range(len(images)):
        ST1 = int(ST[i])
        ST2 = int((ST[i]-ST1)*60.)
        ST3 = (((ST[i]-ST1)*60.)-ST2)*60
        ST_string.append(str(ST1)+':'+str(ST2)+':'+str(ST3))
        tempo_loc[i] = tempo_loc[i].value
        use.update_progress((i+1.)/len(images))
        #print images[i], ' ** ',tempo_loc[i], ' ** ', JD[i], ' ** ', ST[i],' ** ',ST_string[i],' ** ',sun_loc[i],' ** ',HJD[i]
    print '\nSave data file ... \n'
    data = DataFrame([images,tempo_loc,list(JD),list(ST),list(ST_string),list(ra_sun),list(dec_sun),list(dsun),list(HJD)]).T
    data.columns=['images','UTC','JD','ST','ST_isot','RA_SUN','DEC_SUN','D_SUN','HJD']
    print data
    data.to_csv('results.csv')
    os.chdir(original_path)
    return

def time_calibration(input_file):
    """
    Obtain the calibration for time (hjd) by pyraf and the airmass for each image. Include in the header all information.
    """
    original_path = os.getcwd()
    save_path = input_file['save_path']
    #change to save data reduction directory
    os.chdir(save_path)
    print '\n Reading the list of images ....\n'
    planet = input_file['exoplanet'] #set exoplanet name
    images = sorted(glob.glob('AB'+planet+'*.fits'))
    print images
    #include de RA,DEC and epoch of the exoplanet
    RA,DEC,epoch = input_file['RA'],input_file['DEC'],input_file['epoch']
    #obtain ST JD using iraf task and introduce in the header
    for i in range(len(images)):
        hdr = fits.getheader(images[i])
        if int(split(hdr['UT'],':')[0]) < int(hdr['timezone']):
            new_date = use.yesterday(hdr['date-obs'])
            #print images[i], new_date
        else:
            new_date = hdr['date-obs']
        year,month,day = split(new_date,'-')
        iraf.asttimes(year=year,month=month,day=day,time=hdr['loctime'],obs=input_file['observatory'])
        JD = iraf.asttimes.jd #obtain julian date
        LMST = iraf.asttimes.lmst #obtain the sideral time
        LMST = use.sexagesimal_format(LMST) #convert sideral time in sexagesimal format
        iraf.hedit(images[i],'ST',LMST,add='yes',verify='no',show='no',update='yes') #create the ST keyword in the header
        iraf.ccdhedit(images[i],'LMST',LMST,type='string') #include the mean sideral time in the header
        iraf.ccdhedit(images[i],'JD',JD,type='string') #include de julian date in the header
        #include RA, and DEC of the object in your header
        iraf.ccdhedit(images[i],"RA",RA,type="string") #include right ascention in the header
        iraf.ccdhedit(images[i],"DEC",DEC,type="string")  #include declination in the header
        iraf.ccdhedit(images[i],"epoch",epoch,type="string") #include epoch in the header
        # use.update_progress((i+1.)/len(images))
    print '\n Setting airmass ....\n'
    for i in range(len(images)):
        print '# ',images[i]
        #iraf.hedit(images[i],'airmass',airmass,add='yes')
        #iraf.hedit(images[i],'HJD',HJD,add='yes')
        iraf.setairmass.observatory = input_file['observatory']
        iraf.setairmass(images[i])
        iraf.setjd.time = 'ut'
        iraf.setjd(images[i])
    print '\n.... done.\n'
    #export information
    hjd, jd, airmass, st = [],[],[],[]
    for i in range(len(images)):
        hdr = fits.getheader(images[i])
        hjd.append(hdr['HJD'])
        jd.append(hdr['JD'])
        airmass.append(hdr['airmass'])
        st.append(hdr['st'])
    #saving the data
    data = DataFrame([list(hjd),list(jd),list(st),list(airmass)]).T
    data.columns = ['HJD','JD','ST','Airmass']
    data.to_csv('results_iraf_calibrations.csv')
    #change to workings directory
    os.chdir(original_path)
    return

def bkg_info(input_file):
    """
    Obtain the sky backgound for each science image.

    More in: http://photutils.readthedocs.io/en/latest/api/photutils.background.Background.html#photutils.background.Background

    WARNING:

    This routine only need to be run one time for the same set of images.
    ___


    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.
    """
    #set the original directory
    original_path = os.getcwd()
    save_path = input_file['save_path']
    os.chdir(save_path)
    planet = input_file['exoplanet']
    tempo = time.time()
    print 'Obtain background data for each image ... \n'
    if not os.path.exists('background'): #if background does not exist, create!
        os.makedirs('background')

    images = sorted(glob.glob('AB'+planet+'*.fits'))
    #if background exist, check if files bkg_data_image_name_.pik exist
    #if not exist, then create, else: get out of here! XD
    if os.path.exists('background') == True :
        value = []
        for i in images:
            value.append(os.path.isfile('./background/'+'bkg_data_'+i+'_.pik'))
    if (False in value) == True:
        print 'Does not exist all files to all images in the sample.'
        print 'Calculating ...'
        print 'This will take some time... go drink some coffe'
        print ' while you wait for the routine finish \n'
        for i in images:
            im = fits.getdata(i,header=False)
            im = np.array(im,dtype='Float64')
            bkg = Background(im,tuple(input_file['skysection'])) #estimating the background using a boxpixel
            fnpickle(bkg,'./background/'+'bkg_data_'+i+'_.pik')
            use.update_progress((float(images.index(i))+1.)/len(images))
    else:
        print 'The sky background files *.pik exist. \n'
    print 'Sky backgound obtained.'
    print 'Total time = ',abs(time.time()-tempo)/60.,' minutes'
    os.chdir(original_path)
    return

def bkg_read(input_file):
    """
    Read the sky backgound files obtained by bkg_info routine and return the background mask and noise data.
    """
    #set the original directory
    tempo = time.time()
    original_path = os.getcwd()
    save_path = input_file['save_path']
    planet = input_file['exoplanet']
    #change to save data reduction directory
    os.chdir(save_path)
    print 'Reading sky backgound files ...'
    images = sorted(glob.glob('AB'+planet+'*.fits'))
    bkg_data = []
    bkg_rms = []
    for i in range(len(images)):
        bkg = fnunpickle('./background/'+'bkg_data_'+images[i]+'_.pik')
        bkg_data.append(bkg.background)
        bkg_rms.append(bkg.background_rms)
        use.update_progress((i+1.)/len(images))
    os.chdir(original_path)
    print 'total time = ',abs(time.time()-tempo)/60.,' minutes'
    return bkg_data, bkg_rms

def phot_aperture(input_file,bkg_data,bkg_rms):
    """
    Obtain the aperture photometry to the list of apertures in the input_file dictionary.
    ___
    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.
    bkg_data: sky background mask from bkg_info routine
    bkg_rms: sky background noise from bkg_info routine

    """
    #set the original directory
    original_path = os.getcwd()
    save_path = input_file['save_path']
    planet = input_file['exoplanet']
    #radii = np.arange(input_file['apertures'][0],input_file['apertures'][1],0.1)
    radii = np.array(input_file['apertures'])
    #change to save data reduction directory
    os.chdir(save_path)
    if not os.path.exists('phot_results'):
        os.makedirs('phot_results')
    tempo = time.time()
    print 'Starting aperture photometry'
    print 'Saving results on: '+save_path+'/phot_results/'
    #check the number of objects to make the photometry
    N_obj = len(input_file['pxpositions'])/2.
    print 'Number of objects = ',N_obj
    positions = [] #create the positions variable (X,Y) in pixels unit on the CCD
    for i in range(len(input_file['pxpositions'])):
        if i % 2 == 0: #if the number is a even (or not a odd), the turple is created
            positions.append((input_file['pxpositions'][i],input_file['pxpositions'][i+1]))
    print 'Radius from ',radii[0],' to ',radii[-1],'\n'
    images = sorted(glob.glob('AB'+planet+'*.fits'))
    for radius in radii:
        flux_data = []
        for i in range(len(images)):
            im = fits.getdata(images[i],header=False)
            im = array(im,dtype='Float64')
            phot_table = aperture_photometry(im - bkg_data[i], CircularAperture(positions, radius),
                                             error=bkg_rms[i], effective_gain=float(input_file['gain']))
            phot_table_flux = np.array([]) #saving results of aperture photometry
            for j in range(len(phot_table['aperture_sum'])):
                phot_table_flux = np.concatenate((phot_table_flux,np.array([phot_table['aperture_sum'][j]])),axis=0)
                phot_table_flux = np.concatenate((phot_table_flux,np.array([phot_table['aperture_sum_err'][j]])),axis=0)
            flux = np.concatenate((phot_table_flux,np.array([images[i]])),axis=0)
            # flux = [phot_table['aperture_sum'][0], phot_table['aperture_sum'][1],phot_table['aperture_sum_err'][0],
            #         phot_table['aperture_sum_err'][1],images[i]]
            flux_data.append(flux)
        flux_data = DataFrame(flux_data)#,columns=['hoststar','refstar','hoststar_err','refstar_err','image'])
        flux_data.to_csv('./phot_results/'+planet+'_flux_radius_'+str(radius)+'.csv',index=False)
        use.update_progress((float(np.where(radii == radius)[0])+1.)/len(radii))
    print 'Time total = ',abs(time.time()-tempo)/60.,' minutes'
    os.chdir(original_path)

def phot_readData(input_file):
    """
    Read the aperture photometry files and return the normalized flux, rawflux, of the exoplanet, 
    with the error, eflux, the heliocentric julian date, hjd, and teh airmass.
    ___
    
    INPUT
    
    input_file: dict, dictionary of information for you data images sample.
    """
    original_path = os.getcwd()
    os.chdir(input_file['save_path'])
    print 'Reading '+input_file['exoplanet']+'*.csv files ....'
    files_csv = np.sort(glob.glob(input_file['save_path']+'/phot_results/'+input_file['exoplanet']+'*.csv'))
    scatter = np.zeros(len(files_csv))
    for i in range(len(files_csv)):
        phot_data = read_csv(files_csv[i])
        scatter[i] = np.std(phot_data.hoststar)
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
    rawflux = data_min_scatter.hoststar.values/data_min_scatter.refstar.values
    eflux = rawflux*np.sqrt((data_min_scatter.hoststar_err.values/data_min_scatter.hoststar.values)**2 + (data_min_scatter.refstar_err.values/data_min_scatter.refstar.values)**2)
    os.chdir(original_path)
    return rawflux,eflux,hjd,airmass
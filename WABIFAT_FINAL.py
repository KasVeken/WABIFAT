#WABIFAT: Wsclean, Aegean and Bane Interaction to create fluxdensity versus Frequency And Time plots.


#If you want to make a frequency vs fluxdensity plot, run lets_go_freq(). Follow the DO_STEPS in the function!
#If you want to make a time vs fluxdensity plot, run lets_go_time(). Follow the DO_STEPS in the function!

print(" ")
print(" ")
print("WABIFAT: Wsclean, Aegean and Bane Interaction to create fluxdensity versus Frequency And Time plots")
print("WABIFAT version 2020.05.31 ")
print("Author: Kas Veken (veken@strw.leidenuniv.nl)")
print
print("                                        ")
print("                                        ")
print("                           XXXxxx       ")
print("                 .XXXXX:.  XX    X.     ")
print("             .xXXX      XXX  X.   ;XX   ")
print("            :XX.           Xx  X..  XX  ")
print("          .XXx   ++       ++ Xx  X: XX  ")
print("         :Xx.      ++    ++    X  x XX  ")
print("         :X  X   ++        ++  XX: Xx   ")
print("          XXXX   xxxxxxxxxxxx  XX xX.   ")
print("           XX    XX        XX   :.Xx    ")
print("           XX    XX        Xx    XX     ")
print("          .XX    XX===     X:  ;XX      ")
print("           XX.   XX=====  XX  .;X       ")
print("           :X: :x X======XX  .XX.       ")
print("            X X. X XXXXXXX:  XX.        ")
print("       xx   .X X: X         XX.         ")
print("    XXX  XX. 'X X. X       XX           ")
print("   XX  == XXX :X X  X     XX.           ")
print("   xXX  =   .XXx X.  X   XX             ")
print("     ;XXXXx    ;X X  X  :X.             ")
print("         .xxXXXXX  XX   XX:xx;          ")
print("               XXX       ;XXXXXX.       ")
print("             XXXXx            XX        ")
print("           xX      :XXXXx.    xX        ")
print("           ;XX   ;XX;  :XXXXXXX;        ")
print("             XXXXX;                     ")
print("                                        ")
print("                                        ")

import numpy as np
from astropy.io import fits
import os
from matplotlib import pyplot as plt
import os.path
import time
import copy
from astropy import units as u
from astropy.io import fits
import math
import sys
import scipy.ndimage
import re
import sp_module as sp
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
from astropy.coordinates import ICRS
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance

start = time.time()

def paths(lofar_number, starfolder, makedir=True): 
	if makedir == True:
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder)
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/table_folder")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/wscleaned_channels_folder")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/wscleaned_intervals_folder")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/rmsANDbkg_folder")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/results")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badfits")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodfits")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/MFSfits")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/misfits")	
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badtables")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodtables")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badrmsANDbkg")
		os.system("mkdir "+"/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodrmsANDbkg")
		

	star_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder
	table_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/table_folder/"
	wscleaned_channels_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/wscleaned_channels_folder/"
	wscleaned_intervals_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/wscleaned_intervals_folder/"
	rmsANDbkg_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/rmsANDbkg_folder/"
	msfiles_path = "/net/lofar3/data2/veken/msfiles/" 
	results_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/results/"
	badfits_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badfits/"
	goodfits_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodfits/"
	MFSfits_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/MFSfits/"
	misfits_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/misfits/"
	ds9regions_path = "/net/lofar9/data2/veken/ds9regions/" 
	wabifat_folder_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"
	badtables_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badtables/"
	goodtables_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodtables/"
	badrmsANDbkg_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/badrmsANDbkg/"
	goodrmsANDbkg_path = "/net/"+lofar_number+"/data2/veken/wabifat_folder/"+starfolder+"/goodrmsANDbkg/"
	
	return starfolder, table_path, wscleaned_channels_path, wscleaned_intervals_path, rmsANDbkg_path, msfiles_path, results_path, badfits_path, goodfits_path, MFSfits_path, misfits_path, ds9regions_path, wabifat_folder_path, badtables_path, goodtables_path, badrmsANDbkg_path, goodrmsANDbkg_path
starname, table_path, wscleaned_channels_path, wscleaned_intervals_path, rmsANDbkg_path, msfiles_path, results_path, badfits_path, goodfits_path, MFSfits_path, misfits_path, ds9regions_path, wabifat_folder_path, badtables_path, goodtables_path, badrmsANDbkg_path, goodrmsANDbkg_path = paths("lofar9", "GJ3729_time_stokesV_29july", True)	

def wsclean_for_freq(sizeX, sizeY, autothreshold, chan_lower, chan_upper, chan_out, n_iter, stokes, outname, msfile):
	os.system('wsclean  -minuv-l 90 -size '+str(sizeX)+' '+str(sizeY)+' -reorder -channel-range '+str(chan_lower)+' '+str(chan_upper)+' -channels-out '+str(chan_out)+' -weight briggs -0.5  -weighting-rank-filter 3 -clean-border 1 -padding 1.4 -mgain 0.8 -fit-beam -auto-mask 2.5 -auto-threshold '+str(autothreshold)+' -pol '+stokes+' -name '+wscleaned_channels_path+outname+' -scale 1.5arcsec -niter '+str(n_iter)+' '+msfile) 
	#In the wsclean command it returns 5 .fits files: image, model, psf, dirty, residual. Because the last 4 of these are not used in this program anymore, we'd better delete them to save space.
	#We added "-no-dirty" to the wsclean command, so that we dont have to delete it anymore. In case you want to delete, remove "-no-dirty" and add below: +wscleaned_channels_path+'*dirty.fits'
	  

def wsclean_for_time(sizeX, sizeY, autothreshold, int_lower, int_upper, int_out, n_iter, stokes, outname, msfile):
	os.system('wsclean  -minuv-l 90 -size '+str(sizeX)+' '+str(sizeY)+' -reorder  -interval '+str(int_lower)+' '+str(int_upper)+' -intervals-out '+str(int_out)+' -weight briggs -0.5  -weighting-rank-filter 3 -clean-border 1 -padding 1.4 -mgain 0.8 -fit-beam -auto-mask 2.5 -auto-threshold '+str(autothreshold)+' -pol '+stokes+' -name '+wscleaned_intervals_path+outname+' -scale 1.5arcsec -niter '+str(n_iter)+' '+msfile) 	
	#In the wsclean command it returns 5 .fits files: image, model, psf, dirty, residual. Because the last 4 of these are not used in this program anymore, we'd better delete them to save space.
	#We added "-no-dirty" to the wsclean command, so that we dont have to delete it anymore. In case you want to delete, remove "-no-dirty" and add below: +wscleaned_channels_path+'*dirty.fits'
	#os.system('rm '+wscleaned_intervals_path+'*model.fits '+wscleaned_intervals_path+'*residual.fits '+wscleaned_intervals_path+'*psf.fits ')   
#-channel-range '+str(chan_lower)+' '+str(chan_upper)+'



def BANE_and_aegean(FREQ_or_TIME, sigma, ra, dec, coord_margin, outname, wscleaned_path, starname, SNR_accept, misfitrun = False, npixel = 1000, CENTRAL = 999999.99999, negative = '--negative'):
	
	print(FREQ_or_TIME)
	print(misfitrun)
	wscleaned_name = outname+'-image.fits'
	if os.path.isfile(wscleaned_path+wscleaned_name):
		print(wscleaned_path+wscleaned_name, "exists")
	else: 
		print(wscleaned_path+wscleaned_name, "does not exist!")
		print("For some reason, wsclean did probably not make this image.")
		print("We now return peakflux = 0, errpeakflux = 0, rms = 0, central = 0, table_comp_name = failed_table_comp_+outname bkg_name2 = failed_bkg_+outname, rms_name2 = failed_rms_+outname")
		return 0, 0, 0, 0, "failed_table_comp_"+outname , "failed_bkg_"+outname , "failed_rms_"+outname
		
	#These 4 upper and lower bounds define a square around the target coordinates, with a "radius" of 3arcseconds.
	lower_ra = ra - coord_margin  
	upper_ra = ra + coord_margin 
	lower_dec = dec - coord_margin
	upper_dec = dec + coord_margin
	
	#Here we define a bunch of names, which are needed in the BANE and Aegean runs. Make sure you have created the necessary directories
	
	if FREQ_or_TIME == "FREQ":
		hdulist = fits.open(wscleaned_path+wscleaned_name)
		centralfreq_in_Hz = hdulist[0].header['CRVAL3']
		central = centralfreq_in_Hz/(10**6) #Gives back the central frequency of the image in MHz
		fitsname =  wscleaned_name.rstrip('.fits') #Yes, the naming is stupid.
		central = CENTRAL
		os.system(' BANE ' + wscleaned_path+wscleaned_name + ' --out='+rmsANDbkg_path+str(central)+'MHz'+fitsname)
		rms_name = rmsANDbkg_path+str(central)+'MHz'+fitsname+'_rms.fits' #this is the name BANE gives to your rms file.
		bkg_name = rmsANDbkg_path+str(central)+'MHz'+fitsname+'_bkg.fits' #this is the name BANE gives to your bkg file.
		table_name = 'table_'+str(central)+'MHz'+fitsname #this is in the name that Aegean will give to the fits table.

		#The following names are returned by this function, dont pay attention to it for now.
		bkg_name2 = str(central)+'MHz'+fitsname+'_bkg.fits'
		rms_name2 = str(central)+'MHz'+fitsname+'_rms.fits'
		table_comp_name = table_name+'_comp.fits'
	
	if FREQ_or_TIME == "TIME":
		hdulist = fits.open(wscleaned_path+wscleaned_name)
		central_time_in_ISOT = hdulist[0].header['DATE-OBS']
		t = Time(central_time_in_ISOT, format='isot')
		central = (t.jd)*24.0 #Gives back the central time of the image in Julian Date, in hours!
		print('central time for this slice/image is:', central, CENTRAL)
		fitsname =  wscleaned_name.rstrip('.fits')
		print('fitsname:',fitsname)
		central = CENTRAL
		os.system(' BANE ' + wscleaned_path+wscleaned_name + ' --out='+rmsANDbkg_path+str(central)+'hrs'+fitsname)
		rms_name = rmsANDbkg_path+str(central)+'hrs'+fitsname+'_rms.fits' #this is the name BANE gives to your rms file.
		bkg_name = rmsANDbkg_path+str(central)+'hrs'+fitsname+'_bkg.fits' #this is the name BANE gives to your bkg file.
		table_name = 'table_'+str(central)+'hrs'+fitsname #this is in the name that Aegean will give to the fits table.
		
		#The following names are returned by this function, dont pay attention to it for now. The only difference is that it doesnt have paths, and the table has _comp.fits in it.
		bkg_name2 = str(central)+'hrs'+fitsname+'_bkg.fits'
		rms_name2 = str(central)+'hrs'+fitsname+'_rms.fits'
		table_comp_name = table_name+'_comp.fits'
	
	
	#Here we define peakflux as 0, if we detect our targets flux, it will not be zero anymore and we wont return zero. 
	peakflux = 0 
	local_rms = 0
	err_peakflux = 0
	#Here it will run aegean on the current file. After that it will check if it has a peak flux at the necessary coordinates, if so: note it down as peakflux. 
	os.system('aegean ' ' --noise ' + rms_name + ' --background ' + bkg_name + ' ' + wscleaned_path+wscleaned_name + ' --table '+ table_path+table_name+".fits "+negative+' --seedclip '+str(sigma)) 
	flux_within_region=[]
	rms_within_region=[]
	err_flux_within_region=[]
	if os.path.isfile(table_path+table_name+'_comp.fits') and os.path.isfile(rms_name) and os.path.isfile(bkg_name):
		print("all the files exist created by BANE and Aegean")
		fitstable = fits.open(table_path+table_name+'_comp.fits') #Here it reads in the table made by Aegean.
		fits_RA = fitstable[1].data['ra']
		#print('fits_RA', fits_RA)
		fits_DEC = fitstable[1].data['dec']
		for i in range(len(fits_RA)):
			if fits_RA[i] >= lower_ra and fits_RA[i] <= upper_ra and fits_DEC[i] >= lower_dec and fits_DEC[i] <= upper_dec:
				tmp_flux = fitstable[1].data['peak_flux'][i]
				tmp_rms = fitstable[1].data['local_rms'][i]
				tmp_err_flux = fitstable[1].data['err_peak_flux'][i]
				flux_within_region.append(tmp_flux) 
				rms_within_region.append(tmp_rms)
				err_flux_within_region.append(tmp_err_flux)
		if len(flux_within_region) != 0:
			print('len(flux_within_region != 0')
			maxvalue = max(flux_within_region)
			print('maxvalue:', maxvalue)
			minvalue = min(flux_within_region)
			print('minvalue:', minvalue)
			if abs(minvalue) > abs(maxvalue): #Finding the peakflux, which can be positive or negative. 
				print('abs(min)>abs(max)')
				peakflux = minvalue
			else:
				peakflux = maxvalue
			index = flux_within_region.index(peakflux)
			if peakflux >= 0: #this is just to make the rms negative as well if the peakflux is negative. In my opinion then it's easier to compare signal and noise to eachother.
				local_rms = rms_within_region[index]
				print('peakflux >=0, so local rms:',local_rms)
			else: 
				local_rms = rms_within_region[index]*(-1)
				print('peakflux<0, so local rms:', local_rms)
			err_aegean_flux = err_flux_within_region[index]
			print('err_aegean_flux', err_aegean_flux)
			if np.isnan(err_aegean_flux) == True or abs(err_aegean_flux) == 1:
				print("err_aegaen_flux = nan or -1.0")
				err_aegean_flux = 0
			err_peakflux = ((err_aegean_flux**2)+((0.1*abs(peakflux))**2))**(0.5)
			if err_peakflux >= 30.0/1000.0:
				print("err_peakflux >= 30mJy:", err_peakflux)
				err_peakflux = 0.1*abs(peakflux)

			print("abs(5xlocal)" , abs(SNR_accept*local_rms))
			if abs(peakflux) <= abs(SNR_accept*local_rms):
				print("peakflux <= "+str(SNR_accept)+"xlocal_rms")
				peakflux = 0 
				local_rms = 0
				err_peakflux = 0
			
		else: 
			peakflux = 0 
			local_rms = 0
			err_peakflux = 0
			print('I think len(within_region)=0:', len(flux_within_region))
			
		print('achteraf peakflux:', peakflux)
		print('acteraf local_rms:', local_rms)
		print('acterafl err_peakflux:', err_peakflux)
	
		#print(within_region)
		#break #Break here, because in principle this should be the source. 
	else:
		print("CAREFULLLLLL!!!!!!! File: "+table_path+table_name+"_comp.fits does not exist. Apparently there is not even one bright enough source there?") 
	
	

	#If this is a misfit-run and it didnt find a peakflux with high enough SNR, then we want to get the rms from the _rms.fits files created in BANE.
	if misfitrun == True and peakflux == 0 and os.path.isfile(rms_name): 
		x, y, r = sp.regextract(ds9regions_path+starname+"_rms_region_"+str(npixel)+"pixels.reg")
		image = fits.getdata(rms_name)
		rms_flux, rms_eflux, sky, skyerr = sp.aper(image, xc=[x[0]], yc=[y[0]], apr=[r[0]], phpadu=1., skyrad=[-1], exact=True, flux=True, setskyval=0., silent=True)
		print("misfit total flux in rms region:", rms_flux[0])
		misfit_rms = rms_flux[0]/npixel 
		print("misfit rms:", misfit_rms)
		local_rms = misfit_rms[0]



	#Here we remove the table.fits, rms.fits and the bkg.fits files because we never really look at them after running this program. 
	#os.system('rm '+rms_name+' '+bkg_name+' '+table_path+table_name+'_comp.fits') 
	#if peakflux != 0 and local_rms != 0:
	#	SNR = abs(peakflux)/abs(local_rms)
	#	print('SNR this round =', SNR)
	#	print('SNR_accept*local_rms=',SNR_accept*local_rms)
	#	if SNR <= abs(SNR_accept*local_rms):
	#		peakflux = 0 
	#		local_rms = 0
	#		err_peakflux = 0
	
	print('Peakflux for this round =', peakflux)
	print('Local_rms for this round =', local_rms)
	print('Err_peakflux for this round =', err_peakflux)
	print('central =', central)
	

	return peakflux, local_rms, err_peakflux, central, table_comp_name, bkg_name2, rms_name2


def frequency_run(plot = False):
	#DO_STEP 1: Make the paths in the beginning of the program correct. If necessary create new folders.
	
	#DO_STEP 1.5: Make a "large" region (r ~ 50-70") around the coordinates corresponding to your target in ds9. Then make a region anywhere else (I usually pick a point relatively far away from the source) with a small size. This last region has to be created in order to not iterate over a 0-d array in the sp_module we import at the beginning. Save these regions in the ds9regions under the name "[starname]_rms_region_[#pixels in region]pixels.reg", in the ds9regions folder like you can see in the paths-definition, change the regions coordinates from fk5 to physical, this is important! Before you close down ds9, save the number of pixels in the first region and note them down below at 'pixels_rms_region'. It's usefull to save the number of pixels in the name, so that you dont have to check it again every time. 

	#DO_STEP 1.6: Also save the observation date from the header of the previously made fits file. And save it below at 'date_obs'. Save it as 'yyyy-mm-dd' 
	
	#DO_STEP 2: Fill in the following variables. These are constant throughout the whole programm.
	
	#Variables for whole programm:
	#starname = 'CR_Dra0' #Name in final file name. So dont use spaces.
	#starname_coord = 'CR Dra' #Name needed for SkyCoord below. 
	#starfield = 'ELAIS-N1' 


	inputfits = fits.open(wabifat_folder_path+'WABIFAT_input_data_16july.fits')
	input_data = inputfits[1].data[36]
	
	starname = str(input_data['star']) #Name in final file name. So dont use spaces.
	starfield = str(input_data['field'])
	date_obs = str(input_data['lotss_obs_date']) #Get this from the header of your test images. It is the time of observation by LOFAR
	ra = input_data['RA_gaia_to_lotss_epoch']
	dec = input_data['DEC_gaia_to_lotss_epoch']
	msfile = msfiles_path+str(input_data['msfile'])

	


	SNR_accept = 5.0
	
	#Variables for wsclean:
	stokes = 'V'
	negative = '--negative ' #If stokes V: '--negative', if stokes I: ' ' 
	sizeX = 1024
	sizeY = 1024
	n_iter = 10000
	autothreshold = 2
	
	#Variables for BANE and aegean: 
	seedclip = 4.0
	coord_margin = 0.5*0.3*0.002777 #3arcsec in degrees
	ds9region_starname = 'YY_Gem' #This is CR_Dra0 for all the CR_Dra's
	pixels_rms_region = input_data['npixel_for_rms_region']#CR_Dra0: 4939 # LP212-62: 20197


	
	#Initial starting values for particular measurement set, this is true regardless of what range you want to look at:
	chanMHz = input_data['MHz_per_chan'] #[MHz], frequency width of one channel. (Can be calculated by: (endMHz-startMHz)/totalchans = chanMHz, but this went wrong for some reason so just type 0.31 or 0.4 yourself (0.39).
	totalchans = input_data['nr_freq_chans_for_msfile'] #Total amount of channels spaced between startMHz and endMHz in the measurement set. Usually this should be 120, but for e.g. CR Dra it is 200. 
	startMHz = input_data['bandwidth_low'] #[MHz], left boundary of measurement set in MHz.
	endMHz = input_data['bandwidth_up'] #[MHz], right boundary of measurement set in MHz.
	
	 
	
	#Initial starting values for this run:
	lower_chan = 0 #First channel you want to include in the frequency vs flux plot. 
	upper_chan = totalchans #Last channel is exclusive because "np.arange" takes out 1.
	initial_check_width = 5
	filename = 'IC'+str(lower_chan)+'to'+str(upper_chan)+'_n'+str(n_iter)+'size'+str(sizeX)+'x'+str(sizeY)+'star'+starname+'field'+starfield+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes
	#IC = Initial Chan, AT = AutoThreshold, SC = SeedClip. The rest should be obvious. 
	misfitnumber = 2	

	numberlist = [] 
	for i in range(totalchans): #This for loop creates a list full of string type numbers, with 4 indices in total: eg. ['0000', '0001', '0002',..., '0119'].  								
								#This is needed later to read in the wscleaned files 
		if i < 10: 
			numberlist.append('000'+str(i))
		if i >= 10 and i < 100:
			numberlist.append('00'+str(i))
		if i >= 100 and i < 1000:
			numberlist.append('0'+str(i))
		if i >= 1000:
			numberlist.append(str(i))
	print(numberlist)

	print('You are looking at star:', starname, 'field:', starfield)
	print('with sizeX x SizeY:', sizeX, 'x', sizeY)
	print('n_iter =', n_iter)
	print('autothreshold =', autothreshold)
	print('seedclip =',seedclip)
	print('coords (ra,dec) in degrees:', ra, dec)

	#Lists and variables for frequency_run():
	all_chans = np.arange(lower_chan, upper_chan, 1)
	final_flux_list = []
	final_freq_list = []
	misfit_chans = []
	new_all_chans = []
	misfit_check_chans = []
	misfit_cons_series = 1
	final_misfits_chans = []
	final_slice_radius = []
	final_rms_list = []
	final_misfit_radius = []
	final_misfit_rms = []
	final_misfits_freq_list = []
	final_err_flux_list = []
	failed_list = []
	
	#x0 and y0 are made to plot a horizontal line at y=0 later in the code. But they are defined here because we need the the initial channel boundaries.
	x0=np.linspace(startMHz+(lower_chan*chanMHz),startMHz+((upper_chan-1)*chanMHz)+(chanMHz/2.0),2)  #+(chanMHz/2.0)
	y0=[0,0]
	print('x0',x0)
	
	#In the first for-loop WABIFAT runs through the necessary check_widths from the initial given check_width above, to width of the full range. 
	#In the second for-loop WABIFAT runs through the channels, by obeying certain if/else statements it tries to find consecutive series of channels (e.g. 24,25,26,27...). 
	#It will then go through this consecutive series in check_width steps (called Slices), and apply Wsclean, BANE and Aegean on it (by using the third for-loop and more if/else statements) to find the 		#corresponding peakfluxes. 
	#If the peakflux of a Slice has its SNR too low, it will go through the channels again in the next round with a larger check_width. 
	for check_width in range(initial_check_width, (upper_chan-lower_chan+2)): #+2 makes the maximum checkwidth 1 higher than upperscan, you can change the upper limit to limit ur maximum checkwidth.
		series_check_chans = []
		for i in range(len(all_chans)-1): 
			print('all_chans[i]:',all_chans[i],'i:',i,'all_chans[i+1]:',all_chans[i+1])
			if all_chans[i] != all_chans[-2] and all_chans[i] + 1 == all_chans[i+1] : 
				series_check_chans.append(all_chans[i])
				print('first if', i)
				continue
			elif all_chans[i] == all_chans[-2]  or all_chans[i] + 1 != all_chans[i+1]: 
				print('first elif', i)
				series_check_chans.append(all_chans[i])
				if all_chans[i] == all_chans[-2] and all_chans[-2] + 1 == all_chans[-1]:
					series_check_chans.append(all_chans[-1])
				#if all_chans[i] == all_chans[-2] and all_chans[-2] + 1 != all_chans[-1]:
				#	misfit_chans.append(all_chans[-1])
				print('series_check_chans', series_check_chans)
 				cons_series = len(series_check_chans)
				
				if check_width > cons_series: #If this happens we want to add it to the misfits, because it's in principle not gonna give anything good, because we already checked that width. 
											  #Also the modulo will always be 1 for 1%2, 1%3,... and 2 for 2%3... and so on. So it will not work.
					misfit_chans += series_check_chans
					print('checkwidth>cons_series:', check_width, cons_series, i)
					print('misfit_chans 1',misfit_chans)
					series_check_chans = []
					cons_series = 1
					continue 
				else:
					print('normal wsclean else:',i)
					modulo = cons_series % check_width #Modulo is the amount of channels that are misfits
					print('modulo',modulo)
					cons_series -= modulo
					print('cons_series',cons_series)
					print('2nd series_check_chans', series_check_chans)
					tmp_series_array = np.asarray(series_check_chans)
					print('tmp_series_array', tmp_series_array)
					if modulo == 0:
						series_array = tmp_series_array
						cons_series = len(series_array)
					else: 
						series_array = tmp_series_array[:-modulo] #Cuts out the misfits.
						cons_series = len(series_array)
						tmp_misfits_array = tmp_series_array[(len(tmp_series_array)-modulo):]
						tmp_list = tmp_misfits_array.tolist()
						if modulo >= misfitnumber: #If we have more than 1 misfits we add it to the normal cycle of chans to check again, so that we try more combinations. Otherwise the misfits will fill up quickly.
							new_all_chans += tmp_list
						else:
							misfit_chans += tmp_list
						print('misfit_chans 2',misfit_chans)
					print('series_array',series_array)
					series_lower_chan = series_array[0] #This is the first channel of this consecutive series
					series_upper_chan = series_array[-1] #This is the last channel of this consecutive series
					series_chan_out = cons_series/check_width 
					series_outname = 'cw'+str(check_width)+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename
					wsclean_for_freq(sizeX, sizeY, autothreshold, series_lower_chan, series_upper_chan+1, series_chan_out, n_iter, stokes, series_outname, msfile)
					

					os.system('mv '+wscleaned_channels_path+series_outname+'-MFS-image.fits '+MFSfits_path+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename+'-MFS-image.fits')
					os.system('mv '+wscleaned_channels_path+series_outname+'-MFS-dirty.fits '+MFSfits_path+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename+'-MFS-dirty.fits')
					os.system('mv '+wscleaned_channels_path+series_outname+'-MFS-psf.fits '+MFSfits_path+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename+'-MFS-psf.fits')
					os.system('mv '+wscleaned_channels_path+series_outname+'-MFS-residual.fits '+MFSfits_path+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename+'-MFS-residual.fits')
					os.system('mv '+wscleaned_channels_path+series_outname+'-MFS-model.fits '+MFSfits_path+'chan'+str(series_lower_chan)+'to'+str(series_upper_chan)+filename+'-MFS-model.fits')

					series_sliced = np.split(series_array, series_chan_out) #Creates the Slices Wsclean already looked in for a certain check_width. E.g. check_width = 2 so series_chan_out = 6/2=3:[0,1,2,3,4,5]>>[[0,1],[2,3],[4,5]]
					print('series_sliced',series_sliced)
					for j in range(series_chan_out): 
						#If Wsclean only creates one image within a certain range, the image file will NOT have a 4 digit number in its name (e.g. name-image.fits instead of: name-'0001'-image.fits).
						if cons_series == check_width:
							tmp_outname = series_outname
						else:
							tmp_outname = series_outname+'-'+numberlist[j]
						print('tmp_outname',tmp_outname)
						
						Slice = series_sliced[j] 
						tmp_chan_up = Slice[-1] #Last channel in this particular slice.
						tmp_chan_low = Slice[0] #First channel in this particular slice. 
						centralfreq_for_Slice = startMHz+(Slice[0]*chanMHz)+(check_width*(chanMHz/2.0)) #in MHz
						slice_radius = centralfreq_for_Slice-((Slice[0]*chanMHz)+startMHz)
						slice_radius_notgood = check_width*(chanMHz/2.0)

						peakflux, local_rms, err_peakflux, centralfreq_for_Slice_notgood, table_comp_name, bkg_name, rms_name = BANE_and_aegean("FREQ",seedclip, ra, dec, coord_margin, tmp_outname, wscleaned_channels_path, ds9region_starname, SNR_accept, CENTRAL = centralfreq_for_Slice, negative = negative)
						

						
						if peakflux == 0:
							print("Back in frequencyrun(): peakflux == 0 and your slice is:", Slice)
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-image.fits '+badfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-image.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-psf.fits '+badfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-psf.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-model.fits '+badfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-model.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-residual.fits '+badfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-residual.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-dirty.fits '+badfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-dirty.fits')	
							os.system('mv '+table_path+table_comp_name+' '+badtables_path+table_comp_name)	
							os.system('mv '+rmsANDbkg_path+bkg_name+' '+badrmsANDbkg_path+bkg_name)
							os.system('mv '+rmsANDbkg_path+rms_name+' '+badrmsANDbkg_path+rms_name)
							for chan in Slice:
								new_all_chans.append(chan)
							if centralfreq_for_Slice == 0:
								print("for some reason wsclean failed to make an image.")
								listslice = Slice.tolist()
								failed_list += listslice
								failed_list.append(-9999) 
						else:	
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-image.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-image.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-psf.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-psf.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-model.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-model.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-residual.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-residual.fits')
							os.system('mv '+wscleaned_channels_path+tmp_outname+'-dirty.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_Slice)+'MHz'+'chan'+str(tmp_chan_low)+'to'+str(tmp_chan_up)+filename+'-dirty.fits')
							os.system('mv '+table_path+table_comp_name+' '+goodtables_path+table_comp_name)	
							os.system('mv '+rmsANDbkg_path+bkg_name+' '+goodrmsANDbkg_path+bkg_name)
							os.system('mv '+rmsANDbkg_path+rms_name+' '+goodrmsANDbkg_path+rms_name)
							final_flux_list.append(peakflux)
							final_rms_list.append(local_rms)
							final_freq_list.append(centralfreq_for_Slice)
							final_slice_radius.append(slice_radius)
							final_err_flux_list.append(err_peakflux)
					print('new_all_chans', new_all_chans, i)
					series_check_chans = []
					cons_series = 1 
		os.system('rm '+badfits_path+'*.fits')
		os.system('rm '+badtables_path+'*.fits')
		os.system('rm '+badrmsANDbkg_path+'*.fits')
		chaos_new_all_chans = copy.copy(new_all_chans)
		print('chaos new_all_chans', chaos_new_all_chans)
		new_all_chans.sort()
		print('sorted new_all_chans', new_all_chans)
		if len(new_all_chans) == 0: #If this statement is true, then we want to break the for loop because the only remaining channels are in the misfits list.
			break
		else:
			all_chans = new_all_chans
			print('second all_chans',all_chans, check_width)
			new_all_chans = []
			chaos_new_all_chans = []
	 	

	#Here we will run the last check on the misfits to see if adjacent misfits exist, and if so, if they give a significant result.
	#First we have to sort the misfits_chans, because they are not per se appended in increasing order. We copy the unsorted way to chaos_misfits, for safety.
	chaos_misfits = copy.copy(misfit_chans)
	print('chaos_misfits:',chaos_misfits)
	#np.savetxt(results_path+'FREQ_chaos_list'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'.txt', final_flux_list, delimiter=',')
	misfit_chans.sort() 
	print('misfits before last check:',misfit_chans)
	np.savetxt(results_path+'FREQ_misfits_before_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'.txt', misfit_chans, delimiter=',')
	if len(misfit_chans) == 0:	
		print('There are no misfits')
		final_misfits_chans = []	
	else:
		if len(misfit_chans) == 1: #If there is exactly 1 misfit (which will occur almost never, only if you check a small range), then we add another random number such that the following for loop will work
			misfit_chans.append(-1000) #-1000 is a number that will never be the number of a channel because its negative. 	
		print("misfit_chans met -1000:", misfit_chans)
		print('there are 1 or more misfits')
		for y in range(len(misfit_chans)-1):
			print('round', y, misfit_chans[y])
			if misfit_chans[y] != misfit_chans[-2] and misfit_chans[y] + 1 == misfit_chans[y+1]:
				print('first if',y)
				misfit_check_chans.append(misfit_chans[y])
				continue
			elif misfit_chans[y] == misfit_chans[-2] or misfit_chans[y] + 1 != misfit_chans[y+1]:
				print('first elif',y)
				misfit_check_chans.append(misfit_chans[y])
				if misfit_chans[y] == misfit_chans[-2] and misfit_chans[y] + 1 == misfit_chans[-1]:
					misfit_check_chans.append(misfit_chans[-1])
				misfit_cons_series = len(misfit_check_chans)
				if misfit_cons_series == 10000:
					print("something is really wrong with your code mate")
				else: 
					print('misfit_check_chans',misfit_check_chans,y)
					misfit_lower = misfit_check_chans[0]
					misfit_upper = misfit_check_chans[-1]
					misfit_check_width = misfit_cons_series #Because we are just going to try one combination of adjacent misfits, it will be the length of this series.
					misfit_outname = 'misfit_cw'+str(misfit_check_width)+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename
					wsclean_for_freq(sizeX, sizeY, autothreshold, misfit_lower, misfit_upper+1, 1, n_iter, stokes, misfit_outname, msfile) #Note the '1'. 
					centralfreq_for_misfit = startMHz+(misfit_lower*chanMHz)+(misfit_check_width*(chanMHz/2.0))
					misfit_radius = centralfreq_for_misfit-((misfit_lower*chanMHz)+startMHz)
					misfit_radius_notgood = misfit_check_width*(chanMHz/2.0)
					misfit_peakflux, misfit_local_rms, misfit_err_peakflux, centralfreq_for_misfit_notgood, misfit_table_comp_name, misfit_bkg_name, misfit_rms_name = BANE_and_aegean("FREQ", seedclip, ra, dec, coord_margin, misfit_outname, wscleaned_channels_path, ds9region_starname, SNR_accept, True, pixels_rms_region, CENTRAL = centralfreq_for_misfit, negative = negative)
					if misfit_peakflux != 0: #If this is true, then the misfit is not a misfit anymore, and can be incorporated in the 'goodfits'.
						print('misfit_peakflux',misfit_peakflux)
						
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-image.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-psf.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-model.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-residual.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-dirty.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+goodtables_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+goodrmsANDbkg_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+goodrmsANDbkg_path+misfit_rms_name)
						final_flux_list.append(misfit_peakflux)
						final_rms_list.append(misfit_local_rms)
						final_freq_list.append(centralfreq_for_misfit)
						final_slice_radius.append(centralfreq_for_misfit-(misfit_lower*chanMHz)-startMHz)
						final_err_flux_list.append(misfit_err_peakflux)
					else:
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-image.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-psf.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-model.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-residual.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-dirty.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+misfits_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+misfits_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+misfits_path+misfit_rms_name)
						final_misfits_freq_list.append(centralfreq_for_misfit)
						final_misfit_radius.append(misfit_radius)
						final_misfit_rms.append(misfit_local_rms)
						final_misfits_chans += misfit_check_chans
						print('tmp final_misfits_chans:', final_misfits_chans)
						if centralfreq_for_misfit == 0:
							print("for some reason wsclean failed to make an image.")
							failed_list += misfit_check_chans
							failed_list.append(-9999)
					misfit_check_chans = []
				if misfit_chans[y] == misfit_chans[-2] and misfit_chans[y] + 1 != misfit_chans[y+1] and misfit_chans[y+1] != -1000:
					print("The last lonely misfit")
					misfit_check_chans = []
					misfit_check_chans.append(misfit_chans[-1])
					misfit_lower = misfit_check_chans[0]
					misfit_upper = misfit_check_chans[-1]
					misfit_check_width = misfit_cons_series #Because we are just going to try one combination of adjacent misfits, it will be the length of this series.
					misfit_outname = 'misfit_cw'+str(misfit_check_width)+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename
					wsclean_for_freq(sizeX, sizeY, autothreshold, misfit_lower, misfit_upper+1, 1, n_iter, stokes, misfit_outname, msfile) #Note the '1'. 
					centralfreq_for_misfit = startMHz+(misfit_lower*chanMHz)+(misfit_check_width*(chanMHz/2.0))
					misfit_radius = centralfreq_for_misfit-((misfit_lower*chanMHz)+startMHz)
					misfit_radius_notgood = misfit_check_width*(chanMHz/2.0)
					misfit_peakflux, misfit_local_rms, misfit_err_peakflux, centralfreq_for_misfit_notgood, misfit_table_comp_name, misfit_bkg_name, misfit_rms_name = BANE_and_aegean("FREQ", seedclip, ra, dec, coord_margin, misfit_outname, wscleaned_channels_path, ds9region_starname, SNR_accept, True, pixels_rms_region, CENTRAL = centralfreq_for_misfit, negative = negative)
					if misfit_peakflux != 0: #If this is true, then the misfit is not a misfit anymore, and can be incorporated in the 'goodfits'.
						print('misfit_peakflux',misfit_peakflux)
						
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-image.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-psf.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-model.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-residual.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-dirty.fits '+goodfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+goodtables_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+goodrmsANDbkg_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+goodrmsANDbkg_path+misfit_rms_name)
						final_flux_list.append(misfit_peakflux)
						final_rms_list.append(misfit_local_rms)
						final_freq_list.append(centralfreq_for_misfit)
						final_slice_radius.append(centralfreq_for_misfit-(misfit_lower*chanMHz)-startMHz)
						final_err_flux_list.append(misfit_err_peakflux)
					else:
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-image.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-psf.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-model.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-residual.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_channels_path+misfit_outname+'-dirty.fits '+misfits_path+'centralfreq'+str(centralfreq_for_misfit)+'MHz'+'chan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+misfits_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+misfits_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+misfits_path+misfit_rms_name)
						final_misfits_freq_list.append(centralfreq_for_misfit)
						final_misfit_radius.append(misfit_radius)
						final_misfit_rms.append(misfit_local_rms)
						final_misfits_chans += misfit_check_chans
						if centralfreq_for_misfit == 0:
							print("for some reason wsclean failed to make an image.")
							failed_list += misfit_check_scans
							failed_list.append(-9999)
						print('tmp final_misfits_chans:', final_misfits_chans)
						
						
					
	np.savetxt(results_path+'FREQ_failed_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', failed_list, delimiter=',')
	
	#Because WABIFAT goes through the channels from smaller to larger check_width, the found peakfluxes will not be appended to the list in an increasing order of frequency. 
	#In this for-loop the final lists are sorted such that the frequency is increasing from lowest to largest, and the corresponding fluxlist and radiuslist is sorted with it. 
	for z in range(len(final_freq_list)): 
		swap = z + np.argmin(final_freq_list[z:])
		(final_freq_list[z], final_freq_list[swap]) = (final_freq_list[swap], final_freq_list[z])
		(final_flux_list[z], final_flux_list[swap]) = (final_flux_list[swap], final_flux_list[z])
		(final_slice_radius[z], final_slice_radius[swap]) = (final_slice_radius[swap], final_slice_radius[z])
		(final_rms_list[z], final_rms_list[swap]) = (final_rms_list[swap], final_rms_list[z])
		(final_err_flux_list[z], final_err_flux_list[swap]) = (final_err_flux_list[swap], final_err_flux_list[z])
	
	np.savetxt(results_path+'FREQ_slice_radius_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_slice_radius, delimiter=',')
	np.savetxt(results_path+'FREQ_flux_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_flux_list, delimiter=',')
	np.savetxt(results_path+'FREQ_freq_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_freq_list, delimiter=',')
	np.savetxt(results_path+'FREQ_rms_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_rms_list, delimiter=',')
	np.savetxt(results_path+'FREQ_err_flux_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_err_flux_list, delimiter=',')
	print('final_slice_radius:', final_slice_radius)
	print('final_freq_list:', final_freq_list)
	print('final_flux_list:', final_flux_list)
	print('final_rms_list:', final_rms_list)
	print('final_err_flux_list:', final_err_flux_list)

	
	#The misfits are set to zero and plotted as red "x'es" 
	misfits_freq = []
	misfits_zero = []
	if len(final_misfits_chans) != 0:
		np.savetxt(results_path+'FREQ_misfits_chan_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfits_chans, delimiter=',')
		np.savetxt(results_path+'FREQ_final_misfits_centralfreq_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfits_freq_list, delimiter=',')
		np.savetxt(results_path+'FREQ_final_misfit_radius_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfit_radius, delimiter=',')
		np.savetxt(results_path+'FREQ_final_misfit_rms_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfit_rms, delimiter=',')
		
	print('final_misfits_chan_list:', final_misfits_chans)
	print('final_misfits_centralfreq_list:', final_misfits_freq_list)
	print('final_misfit_radius',final_misfit_radius)
	print('final_misfit_rms', final_misfit_rms)

	if plot == True: #Better just plot yourself afterwards in a/the separate plot.py file.
		fluxarray = np.asarray(final_flux_list)
		freqarray = np.asarray(final_freq_list)
		rmsarray = np.asarray(final_rms_list)
		errfluxarray = np.asarray(final_err_flux_list)
		misfit_rmsarray = np.asarray(final_misfit_rms)
		misfit_centralfreqarray = np.asarray(final_misfits_freq_list)


		
		"""
		#In the following code notice that flux*1000.0 is to put it in mJy/beam, and the rmsarray*100 is to create an arrow 10% the size of the value of the rms, also in mJy/beam (*1000*0.1)
		plt.figure()
		plt.plot(x0,y0,c='green')
		plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, uplims=True, lolims=True, c = 'blue', label = str(len(fluxarray))+' data points', fmt='.') 
		plt.errorbar(freqarray,rmsarray*1000.0, xerr = final_slice_radius, yerr=rmsarray*100, uplims=True, c = 'aqua', label = 'detection local RMS',fmt='.') 
		plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0, xerr = final_misfit_radius, yerr=misfit_rmsarray*100, uplims=True, c = 'red', label = 'non-detection RMS',fmt='.') 
		plt.scatter(misfits_freq, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_freq))+' misfit channels')
		plt.title('Frequency VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		#plt.savefig(results_path+'freq_vs_flux_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip))))
		"""
		plt.plot(x0,y0,c='green')
		plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(freqarray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_freq, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_freq))+' misfit channels')
		plt.title('Frequency VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		
	
	end = time.time()	
	print('Time occured:'+str((end - start)/3600.0)+'hrs')	
	print
	print("Thank you for using WABIFAT, have a nice day!")
	print("                                        ")
	print("                                        ")
	print("                           XXXxxx       ")
	print("                 .XXXXX:.  XX    X.     ")
	print("             .xXXX      XXX  X.   ;XX   ")
	print("            :XX.           Xx  X..  XX  ")
	print("          .XXx   ++       ++ Xx  X: XX  ")
	print("         :Xx.      ++    ++    X  x XX  ")
	print("         :X  X   ++        ++  XX: Xx   ")
	print("          XXXX   xxxxxxxxxxxx  XX xX.   ")
	print("           XX    XX        XX   :.Xx    ")
	print("           XX    XX        Xx    XX     ")
	print("          .XX    XX===     X:  ;XX      ")
	print("           XX.   XX=====  XX  .;X       ")
	print("           :X: :x X======XX  .XX.       ")
	print("            X X. X XXXXXXX:  XX.        ")
	print("       xx   .X X: X         XX.         ")
	print("    XXX  XX. 'X X. X       XX           ")
	print("   XX  == XXX :X X  X     XX.           ")
	print("   xXX  =   .XXx X.  X   XX             ")
	print("     ;XXXXx    ;X X  X  :X.             ")
	print("         .xxXXXXX  XX   XX:xx;          ")
	print("               XXX       ;XXXXXX.       ")
	print("             XXXXx            XX        ")
	print("           xX      :XXXXx.    xX        ")
	print("           ;XX   ;XX;  :XXXXXXX;        ")
	print("             XXXXX;                     ")
	print("                                        ")
	print("                                        ")
	



def time_run(plot = False):
	#DO_STEP 1: Make the paths in the beginning of the program correct. If necessary create new folders.
	
	#DO_STEP 1.5: Make a "large" region (r ~ 50-70") around the coordinates corresponding to your target in ds9. Then make a region anywhere else (I usually pick a point relatively far away from the source) with a small size. This last region has to be created in order to not iterate over a 0-d array in the sp_module we import at the beginning. Save these regions in the ds9regions under the name "[starname]_rms_region_[#pixels in region]pixels.reg", in the ds9regions folder like you can see in the paths-definition, change the regions coordinates from fk5 to physical, this is important! Before you close down ds9, save the number of pixels in the first region and note them down below at 'pixels_rms_region'. It's usefull to save the number of pixels in the name, so that you dont have to check it again every time.

	#DO_STEP 1.6: Also save the observation date from the header of the previously made fits file. And save it below at 'date_obs'. Save it as 'yyyy-mm-dd' 
	
	#DO_STEP 2: Fill in the following variables. These are constant throughout the whole programm.
	inputfits = fits.open(wabifat_folder_path+'WABIFAT_input_data_16july.fits')
	input_data = inputfits[1].data[41]
	
	starname = str(input_data['star']) #Name in final file name. So dont use spaces.
	starfield = str(input_data['field'])
	#date_obs = str(input_data['lotss_obs_date']) #Get this from the header of your test images. It is the time of observation by LOFAR
	ra = input_data['RA_gaia_to_lotss_epoch']
	dec = input_data['DEC_gaia_to_lotss_epoch']
	msfile = msfiles_path+str(input_data['msfile'])
	
	#Variables for wsclean:
	SNR_accept = 5.0
	stokes = 'V'
	negative = '--negative ' #if stokes V: ' ', if stokes I: '--negative'
	sizeX = 1600
	sizeY = 1600 
	n_iter = 10000
	autothreshold = 2
	
	#Variables for BANE and aegean: 
	seedclip = 4.0
	coord_margin = 0.5*0.3*0.002777 #3arcsec in degrees
	ds9region_starname = 'GJ3729' #THIS IS 'CR_Dra0' for all the CR_Dra's, so be careful!
	pixels_rms_region = input_data['npixel_for_rms_region']#CR_Dra0: 6843 # LP212-62: 20197

	#Initial starting values for particular measurement set, this is true regardless of what range you want to look at:
	scansec = input_data['sec_per_chan']/3600.0 #[sec], Time in one scan in sec. (Can be calculated by: (endsec-startsec)/totalscans = scansec, but it should for now allways be 16s. 
	print("scansec", scansec)
	totalscans = input_data['nr_time_chans_for_msfile'] #Total amount of scans spaced between startsec and endsec in the measurement set. Usually this should be 1798, but for e.g. CR Dra it is 1789. 1723 for CR_Dra10to20 i think.
	misfitnumber = 2

	#This is actually not needed anymore.
	startsec = 0.0/3600.0 #[hr], left boundary of measurement set in sec/3600 so in hours.
	endsec = (scansec*totalscans)/3600.0 #CRDRA0:28624  #LP212-62:28768 #[sec], right boundary of measurement set in sec.
	

	
	##Initial starting values for this run:
	lower_scan = 0 #First scan you want to include in the time vs flux plot. 
	upper_scan = totalscans #Last scan is exclusive because "np.arange" takes out 1.
	initial_check_width = 325
	filename = 'IS'+str(lower_scan)+'to'+str(upper_scan)+'_n'+str(n_iter)+'size'+str(sizeX)+'x'+str(sizeY)+'star'+starname+'field'+starfield+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes
	#IS = Initial Scan, AT = AutoThreshold, SC = SeedClip. The rest should be obvious. 

	numberlist = [] 
	for i in range(totalscans): #This for loop creates a list full of string type numbers, with 4 indices in total: eg. ['0000', '0001', '0002',..., '0119'].  
								#This is needed later to read in the wscleaned files 
		if i < 10: 
			numberlist.append('000'+str(i))
		if i >= 10 and i < 100:
			numberlist.append('00'+str(i))
		if i >= 100 and i < 1000:
			numberlist.append('0'+str(i))
		if i >= 1000:
			numberlist.append(str(i))
	#print(numberlist)

	print('You are looking at star:', starname, 'field:', starfield)
	print('with sizeX x SizeY:', sizeX, 'x', sizeY)
	print('n_iter =', n_iter)
	print('autothreshold =', autothreshold)
	print('seedclip =',seedclip)
	print('coords (ra,dec) in degrees:', ra, dec)

	#Lists and variables for time_run():
	all_scans = np.arange(lower_scan, upper_scan, 1) #Creates array consisting of all the scans you want (e.g.[0,1,2,3,...,1797,1798])
	final_flux_list = []
	final_time_list = []
	misfit_scans = []
	new_all_scans = []
	misfit_check_scans = []
	misfit_cons_series = 1
	final_misfits_scans = []
	final_slice_radius = []
	final_rms_list = []
	final_misfit_radius = []
	final_misfit_rms = []
	final_misfits_time_list = []
	final_err_flux_list = []
	failed_list=[]
	
	#x0 and y0 are made to plot a horizontal line at y=0 later in the code. But they are defined here because we need the the initial scan boundaries.
	x0=np.linspace(startsec+(lower_scan*scansec),startsec+((upper_scan-1)*scansec)+(scansec/2.0),2)  #+(scansec/2.0)
	y0=[0,0]
	print('x0',x0)
	
	#In the first for-loop WABIFAT runs through the necessary check_widths from the initial given check_width above, to width of the full range. 
	#In the second for-loop WABIFAT runs through the scans. By obeying certain if/else statements it tries to find consecutive series of scans (e.g. 24,25,26,27...). 
	#It will then go through this consecutive series in check_width steps (called Slices), and apply Wsclean, BANE and Aegean on it (by using the third for-loop and more if/else statements) to find the 		#corresponding peakfluxes. 
	#If the peakflux of a Slice has its SNR too low, it will go through the scans again in the next round with a larger check_width. 
	for check_width in range(initial_check_width, (upper_scan-lower_scan+2)): #+2 makes the maximum checkwidth 1 higher than upperscan, you can change the upper limit to limit ur maximum checkwidth.
		series_check_scans = []
		for i in range(len(all_scans)-1): #all_scans is e.g.: [0 1 2 ... 1796 1797] for upper_scan = 1798. So this for-loop goes through the numbers 0, 1, 2, ..., 1795, 1796
			print('all_scans[i]:',all_scans[i],'i:',i,'all_scans[i+1]:',all_scans[i+1])
			if all_scans[i] != all_scans[-2] and all_scans[i] + 1 == all_scans[i+1]: #if true, this scan's right-adjacent is itself + 1, and NOT the last scan in all_scans.
				series_check_scans.append(all_scans[i])
				print('first if', i)
				continue
			elif all_scans[i] == all_scans[-2] or all_scans[i] + 1 != all_scans[i+1]: 
				print('first elif', i)
				series_check_scans.append(all_scans[i])
				if all_scans[i] == all_scans[-2] and all_scans[-2] + 1 == all_scans[-1]:
					series_check_scans.append(all_scans[-1])
				print('series_check_scans', series_check_scans)
 				cons_series = len(series_check_scans)
				
				if check_width > cons_series: #If this happens we want to add it to the misfits, because it's in principle not gonna give anything good, because we probalby already checked that width. 
											  #Also the modulo will always be 1 for 1%2, 1%3,... and 2 for 2%3... and so on. So it will not work.
					misfit_scans += series_check_scans
					print('checkwidth>cons_series:', check_width, cons_series, i)
					print('misfit_scans 1',misfit_scans)
					series_check_scans = []
					cons_series = 1
					continue 
				else:
					print('normal wsclean else:',i)
					modulo = cons_series % check_width #Modulo is the amount of scans that are misfits
					print('modulo',modulo)
					cons_series -= modulo
					print('cons_series',cons_series)
					print('2nd series_check_scans', series_check_scans)
					tmp_series_array = np.asarray(series_check_scans)
					print('tmp_series_array', tmp_series_array)
					if modulo == 0:
						series_array = tmp_series_array
						cons_series = len(series_array)
					else: 
						series_array = tmp_series_array[:-modulo] #Cuts out the misfits.
						cons_series = len(series_array)
						tmp_misfits_array = tmp_series_array[(len(tmp_series_array)-modulo):]
						tmp_list = tmp_misfits_array.tolist()
						if modulo >= misfitnumber: #If we have more than 1 misfits we add it to the normal cycle of scans to check again, so that we try more combinations. Otherwise the misfits will fill up quickly.
							new_all_scans += tmp_list
						else:
							misfit_scans += tmp_list
						print('misfit_scans 2',misfit_scans)
					print('series_array',series_array)
					series_lower_scan = series_array[0] #This is the first scan of this consecutive series
					series_upper_scan = series_array[-1] #This is the last scan of this consecutive series
					series_scan_out = cons_series/check_width 
					series_outname = 'cw'+str(check_width)+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename
					wsclean_for_time(sizeX, sizeY, autothreshold, series_lower_scan, series_upper_scan+1, series_scan_out, n_iter, stokes, series_outname, msfile)
					

					os.system('mv '+wscleaned_intervals_path+series_outname+'-MFS-image.fits '+MFSfits_path+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename+'-MFS-image.fits')
					os.system('mv '+wscleaned_intervals_path+series_outname+'-MFS-dirty.fits '+MFSfits_path+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename+'-MFS-dirty.fits')
					os.system('mv '+wscleaned_intervals_path+series_outname+'-MFS-psf.fits '+MFSfits_path+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename+'-MFS-psf.fits')
					os.system('mv '+wscleaned_intervals_path+series_outname+'-MFS-residual.fits '+MFSfits_path+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename+'-MFS-residual.fits')
					os.system('mv '+wscleaned_intervals_path+series_outname+'-MFS-model.fits '+MFSfits_path+'scan'+str(series_lower_scan)+'to'+str(series_upper_scan)+filename+'-MFS-model.fits')

					series_sliced = np.split(series_array, series_scan_out) #Creates the Slices Wsclean already looked in for a certain check_width. E.g. check_width = 2:[0,1,2,3,4,5]>>[[0,1],[2,3],[4,5]]
					print('series_sliced',series_sliced)
					for j in range(series_scan_out): 
						#If Wsclean only creates one image within a certain range, the image file will NOT have a 4 digit number in its name (e.g. name-image.fits instead of: name-'t0001'-image.fits).
						if cons_series == check_width:
							tmp_outname = series_outname
						else:
							tmp_outname = series_outname+'-t'+numberlist[j]
						print('tmp_outname',tmp_outname)
						
						Slice = series_sliced[j] 
						tmp_scan_up = Slice[-1] #Last scan in this particular slice.
						tmp_scan_low = Slice[0] #First scan in this particular slice. 
						centraltime_for_Slice = startsec+(Slice[0]*scansec)+(check_width*(scansec/2.0)) #in hours
						slice_radius = centraltime_for_Slice-((Slice[0]*scansec)+startsec)
						slice_radius_notgood = check_width*(scansec/2.0)
						print("slice_radius for slice:", Slice, "is: ", slice_radius)
						peakflux, local_rms, err_peakflux, centraltime_for_Slice_notgood, table_comp_name, bkg_name, rms_name = BANE_and_aegean("TIME",seedclip, ra, dec, coord_margin, tmp_outname, wscleaned_intervals_path, ds9region_starname, SNR_accept, CENTRAL = centraltime_for_Slice, negative = negative)
						
						if peakflux == 0:
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-image.fits '+badfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-image.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-psf.fits '+badfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-psf.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-model.fits '+badfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-model.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-residual.fits '+badfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-residual.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-dirty.fits '+badfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-dirty.fits')
							os.system('mv '+table_path+table_comp_name+' '+badtables_path+table_comp_name)	
							os.system('mv '+rmsANDbkg_path+bkg_name+' '+badrmsANDbkg_path+bkg_name)
							os.system('mv '+rmsANDbkg_path+rms_name+' '+badrmsANDbkg_path+rms_name)
							for scan in Slice:
								new_all_scans.append(scan)
							if centraltime_for_Slice == 0:
								print("for some reason wsclean failed to make an image.")
								listslice = Slice.tolist()
								failed_list += listslice
								failed_list.append(-9999)
						else:	
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-image.fits '+goodfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-image.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-psf.fits '+goodfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-psf.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-model.fits '+goodfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-model.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-residual.fits '+goodfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-residual.fits')
							os.system('mv '+wscleaned_intervals_path+tmp_outname+'-dirty.fits '+goodfits_path+'centraltime'+str(centraltime_for_Slice)+'hrs'+'scan'+str(tmp_scan_low)+'to'+str(tmp_scan_up)+filename+'-dirty.fits')	
							os.system('mv '+table_path+table_comp_name+' '+goodtables_path+table_comp_name)	
							os.system('mv '+rmsANDbkg_path+bkg_name+' '+goodrmsANDbkg_path+bkg_name)
							os.system('mv '+rmsANDbkg_path+rms_name+' '+goodrmsANDbkg_path+rms_name)
							final_flux_list.append(peakflux)
							final_rms_list.append(local_rms)
							final_time_list.append(centraltime_for_Slice)
							print("final_slice_radius now:", final_slice_radius)
							final_slice_radius.append(slice_radius)
							final_err_flux_list.append(err_peakflux)
					print('new_all_scans', new_all_scans, i)
					series_check_scans = []
					cons_series = 1 
		os.system('rm '+badfits_path+'*.fits')
		os.system('rm '+badtables_path+'*.fits')
		os.system('rm '+badrmsANDbkg_path+'*.fits')
		chaos_new_all_scans = copy.copy(new_all_scans)
		print('chaos new_all_scans', chaos_new_all_scans)
		new_all_scans.sort()
		print('sorted new_all_scans', new_all_scans)
		if len(new_all_scans) == 0: #If this statement is true, then we want to break the checkwidth-for-loop because the only remaining scans are in the misfits list.
			break
		else:
			all_scans = new_all_scans
			print('second all_scans',all_scans, check_width)
			new_all_scans = []
			chaos_new_all_scans = []
	
	#os.system('rm '+wscleaned_intervals_path+'*model.fits '+wscleaned_intervals_path+'*residual.fits '+wscleaned_intervals_path+'*psf.fits ') 	

	#Here we will run the last check on the misfits to see if adjacent misfits exist, and if so, if they give a significant result.
	#First we have to sort the misfits_scans, because they are not per se appended in increasing order. We copy the unsorted way to chaos_misfits, for safety.
	chaos_misfits = copy.copy(misfit_scans)
	print('chaos_misfits:',chaos_misfits)
	#np.savetxt(results_path+'TIME_chaos_list'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'.txt', chaos_misfits, delimiter=',')
	misfit_scans.sort() 
	print('misfits before last check:',misfit_scans)
	np.savetxt(results_path+'TIME_misfits_before_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'.txt', misfit_scans, delimiter=',')
	#if len(misfit_scans)  < initial_check_width + 1:
	if len(misfit_scans) == 0:	
		print('There are no misfits')
		final_misfits_scans = []
	else:
		if len(misfit_scans) == 1: #If there is exactly 1 misfit (which will occur almost never, only if you check a small range), then we add another random number such that the following for loop will work.
			misfit_scans.append(-1000) #-1000 is a number that will never be the number of a channel because its negative. 
		print('there are 1 or more than 1 misfits')
		for y in range(len(misfit_scans)-1):
			print('round', y, misfit_scans[y])
			if misfit_scans[y] != misfit_scans[-2] and misfit_scans[y] + 1 == misfit_scans[y+1]:
				print('first if',y)
				misfit_check_scans.append(misfit_scans[y])
				continue
			elif misfit_scans[y] == misfit_scans[-2] or misfit_scans[y] + 1 != misfit_scans[y+1]:
				print('first elif',y)
				misfit_check_scans.append(misfit_scans[y])
				if misfit_scans[y] == misfit_scans[-2] and misfit_scans[y] + 1 == misfit_scans[-1]:
					misfit_check_scans.append(misfit_scans[-1])
				misfit_cons_series = len(misfit_check_scans)
				if misfit_cons_series == 100000:
					print("something is really wrong with your code mate")
				else: 
					print('misfit_check_scans',misfit_check_scans,y)
					misfit_lower = misfit_check_scans[0]
					misfit_upper = misfit_check_scans[-1]
					misfit_check_width = misfit_cons_series #Because we are just going to try one combination of adjacent misfits, it will be the length of this series.
					misfit_outname = 'misfit_cw'+str(misfit_check_width)+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename
					wsclean_for_time(sizeX, sizeY, autothreshold, misfit_lower, misfit_upper+1, 1, n_iter, stokes, misfit_outname, msfile) #Note the '1'. 
					centraltime_for_misfit = startsec+(misfit_lower*scansec)+(misfit_check_width*(scansec/2.0))
					misfit_radius = centraltime_for_misfit-((misfit_lower*scansec)+startsec)
					misfit_radius_notgood = misfit_check_width*(scansec/2.0)
					print("misfit_radius for misfit:", misfit_check_scans, "is: ", misfit_radius)
					misfit_peakflux, misfit_local_rms, misfit_err_peakflux, centraltime_for_misfit_notgood, misfit_table_comp_name, misfit_bkg_name, misfit_rms_name = BANE_and_aegean("TIME", seedclip, ra, dec, coord_margin, misfit_outname, wscleaned_intervals_path, ds9region_starname, SNR_accept,  True, pixels_rms_region, CENTRAL = centraltime_for_misfit, negative = negative)
					
					if misfit_peakflux != 0: #If this is true, then the misfit is not a misfit anymore, and can be incorporated in the 'goodfits'.	
						print('misfit_peakflux',misfit_peakflux)
						
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-image.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-psf.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-model.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-residual.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-dirty.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+goodtables_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+goodrmsANDbkg_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+goodrmsANDbkg_path+misfit_rms_name)
						final_flux_list.append(misfit_peakflux)
						final_rms_list.append(misfit_local_rms)
						final_time_list.append(centraltime_for_misfit)
						final_slice_radius.append(misfit_radius)
						print("final_slice_radius now:", final_slice_radius)
						final_err_flux_list.append(misfit_err_peakflux)
					else:
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-image.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-psf.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-model.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-residual.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-dirty.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+misfits_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+misfits_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+misfits_path+misfit_rms_name)
						final_misfits_time_list.append(centraltime_for_misfit)
						final_misfit_radius.append(misfit_radius)
						print("final_misfit_radius now:", final_misfit_radius)
						final_misfit_rms.append(misfit_local_rms)
						final_misfits_scans += misfit_check_scans
						if centraltime_for_misfit == 0:
							print("for some reason wsclean failed to make an image.")
							failed_list += misfit_check_scans
							failed_list.append(-9999)

						print('tmp final_misfits_scans:', final_misfits_scans)
						#os.system('rm '+wscleaned_intervals_path+misfit_outname+'-image.fits ') 
					misfit_check_scans = []
				if misfit_scans[y] == misfit_scans[-2] and misfit_scans[y] + 1 != misfit_scans[-1] and misfit_scans[y+1] != -1000:
					print("The last lonely misfit")
					misfit_check_scans = []
					misfit_check_scans.append(misfit_scans[-1])
					misfit_cons_series = len(misfit_check_scans) #is just 1.
					misfit_lower = misfit_check_scans[0] #is the same as the last element because there is just 1 scan left.
					misfit_upper = misfit_check_scans[-1] #is the same as the first element because there is just 1 scan left.
					misfit_check_width = misfit_cons_series #Because we are just going to try one combination of adjacent misfits, it will be the length of this series.
					misfit_outname = 'misfit_cw'+str(misfit_check_width)+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename
					wsclean_for_time(sizeX, sizeY, autothreshold, misfit_lower, misfit_upper+1, 1, n_iter, stokes, misfit_outname, msfile) #Note the '1'. 
					centraltime_for_misfit = startsec+(misfit_lower*scansec)+(misfit_check_width*(scansec/2.0))
					misfit_radius = centraltime_for_misfit-((misfit_lower*scansec)+startsec)
					misfit_radius_notgood = misfit_check_width*(scansec/2.0)
					misfit_peakflux, misfit_local_rms, misfit_err_peakflux, centraltime_for_misfit_notgood, misfit_table_comp_name, misfit_bkg_name, misfit_rms_name = BANE_and_aegean("TIME", seedclip, ra, dec, coord_margin, misfit_outname, wscleaned_intervals_path, ds9region_starname, SNR_accept,  True, pixels_rms_region, CENTRAL = centraltime_for_misfit, negative = negative)
					
					if misfit_peakflux != 0: #If this is true, then the misfit is not a misfit anymore, and can be incorporated in the 'goodfits'.	
						print('misfit_peakflux',misfit_peakflux)
						
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-image.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-psf.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-model.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-residual.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-dirty.fits '+goodfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+goodtables_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+goodrmsANDbkg_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+goodrmsANDbkg_path+misfit_rms_name)
						final_flux_list.append(misfit_peakflux)
						final_rms_list.append(misfit_local_rms)
						final_time_list.append(centraltime_for_misfit)
						final_slice_radius.append(misfit_radius)
						final_err_flux_list.append(misfit_err_peakflux)
					else:
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-image.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-image.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-psf.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-psf.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-model.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-model.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-residual.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-residual.fits')
						os.system('mv '+wscleaned_intervals_path+misfit_outname+'-dirty.fits '+misfits_path+'centraltime'+str(centraltime_for_misfit)+'hrs'+'scan'+str(misfit_lower)+'to'+str(misfit_upper)+filename+'-dirty.fits')
						os.system('mv '+table_path+misfit_table_comp_name+' '+misfits_path+misfit_table_comp_name)	
						os.system('mv '+rmsANDbkg_path+misfit_bkg_name+' '+misfits_path+misfit_bkg_name)
						os.system('mv '+rmsANDbkg_path+misfit_rms_name+' '+misfits_path+misfit_rms_name)
						final_misfits_time_list.append(centraltime_for_misfit)
						final_misfit_radius.append(misfit_radius)
						final_misfit_rms.append(misfit_local_rms)
						final_misfits_scans += misfit_check_scans
						if centraltime_for_misfit == 0:
							print("for some reason wsclean failed to make an image.")
							failed_list += misfit_check_scans
							failed_list.append(-9999)
						print('tmp final_misfits_scans:', final_misfits_scans)
						


	np.savetxt(results_path+'TIME_failed_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', failed_list, delimiter=',')

	#Because WABIFAT goes through the scans from smaller to larger check_width, the found peakfluxes will not be appended to the list in a increasing order of time. 
	#In this for-loop the final lists are sorted such that the time is increasing from early to late, and the corresponding fluxlist and radiuslist are sorted with it. 
	for z in range(len(final_time_list)): 
		swap = z + np.argmin(final_time_list[z:])
		(final_time_list[z], final_time_list[swap]) = (final_time_list[swap], final_time_list[z])
		(final_flux_list[z], final_flux_list[swap]) = (final_flux_list[swap], final_flux_list[z])
		(final_slice_radius[z], final_slice_radius[swap]) = (final_slice_radius[swap], final_slice_radius[z])
		(final_rms_list[z], final_rms_list[swap]) = (final_rms_list[swap], final_rms_list[z])
		(final_err_flux_list[z], final_err_flux_list[swap]) = (final_err_flux_list[swap], final_err_flux_list[z])
	
	np.savetxt(results_path+'TIME_slice_radius_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_slice_radius, delimiter=',')
	np.savetxt(results_path+'TIME_flux_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_flux_list, delimiter=',')
	np.savetxt(results_path+'TIME_time_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_time_list, delimiter=',')
	np.savetxt(results_path+'TIME_rms_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_rms_list, delimiter=',')
	np.savetxt(results_path+'TIME_err_flux_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_err_flux_list, delimiter=',')
	print('final_slice_radius:', final_slice_radius)
	print('final_time_list:', final_time_list)
	print('final_flux_list:', final_flux_list)
	print('final_rms_list:', final_rms_list)
	print('final_err_flux_list:', final_err_flux_list)

	#The misfits are set to zero and plotted as red "x'es" 
	misfits_time = []
	misfits_zero = []
	if len(final_misfits_scans) != 0:
		np.savetxt(results_path+'TIME_misfits_scan_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfits_scans, delimiter=',')	
		np.savetxt(results_path+'TIME_final_misfits_centraltime_list_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfits_time_list, delimiter=',')
		np.savetxt(results_path+'TIME_final_misfit_radius_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfit_radius, delimiter=',')
		np.savetxt(results_path+'TIME_final_misfit_rms_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip)))+'pol'+stokes+'.txt', final_misfit_rms, delimiter=',')


	print('final_misfits_scanlist:', final_misfits_scans)
	print('final_misfits_centraltime_list:', final_misfits_time_list)
	print('final_misfit_radius',final_misfit_radius)
	print('final_misfit_rms', final_misfit_rms)


	if plot == True: #Better just plot yourself afterwards in a/the separate plot.py file.
		fluxarray = np.asarray(final_flux_list)
		timearray = np.asarray(final_time_list)
		rmsarray = np.asarray(final_rms_list)
		errfluxarray = np.asarray(final_err_flux_list)
		misfit_rmsarray = np.asarray(final_misfit_rms)
		misfit_centraltimearray = np.asarray(final_misfits_time_list)


		end = time.time()	
		print('Time occured:'+str((end - start)/3600)+'hrs')	

		plt.figure()
		plt.plot(x0,y0,c='green')
		plt.errorbar(timearray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(timearray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_time, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_time))+' misfit channels')
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		
		"""
		#In the following code notice that flux*1000.0 is to put it in mJy/beam, and the rmsarray*100 is to create an arrow 10% the size of the value of the rms, also in mJy/beam (*1000*0.1)
		plt.figure()
		plt.plot(x0,y0,c='green')
		plt.errorbar(timearray,fluxarray*1000.0, xerr = final_slice_radius, yerr=errfluxarray*1000, uplims=True, lolims=True, c = 'blue', label = str(len(fluxarray))+' data points', fmt='.') 
		plt.errorbar(timearray,rmsarray*1000.0*5, xerr = final_slice_radius, yerr=rmsarray*100, uplims=True, c = 'aqua', label = 'detection 5x local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray,misfit_rmsarray*1000.0*5, xerr = final_misfit_radius, yerr=misfit_rmsarray*100, uplims=True, c = 'red', label = 'non-detection 5x RMS',fmt='.') 
		plt.scatter(misfits_time, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_time))+' misfit channels')
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Time (hours)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		#plt.savefig(results_path+'time_vs_flux_'+starname+'_field'+starfield+'_n'+str(n_iter)+'AT'+str(autothreshold)+'SC'+str(int(round(seedclip))))
		plt.show()		
		"""
	end = time.time()	
	print('Time occured:'+str((end - start)/3600.0)+'hrs')
	print
	print("Thank you for using WABIFAT, have a nice day!")
	print("                                        ")
	print("                                        ")
	print("                           XXXxxx       ")
	print("                 .XXXXX:.  XX    X.     ")
	print("             .xXXX      XXX  X.   ;XX   ")
	print("            :XX.           Xx  X..  XX  ")
	print("          .XXx   ++       ++ Xx  X: XX  ")
	print("         :Xx.      ++    ++    X  x XX  ")
	print("         :X  X   ++        ++  XX: Xx   ")
	print("          XXXX   xxxxxxxxxxxx  XX xX.   ")
	print("           XX    XX        XX   :.Xx    ")
	print("           XX    XX        Xx    XX     ")
	print("          .XX    XX===     X:  ;XX      ")
	print("           XX.   XX=====  XX  .;X       ")
	print("           :X: :x X======XX  .XX.       ")
	print("            X X. X XXXXXXX:  XX.        ")
	print("       xx   .X X: X         XX.         ")
	print("    XXX  XX. 'X X. X       XX           ")
	print("   XX  == XXX :X X  X     XX.           ")
	print("   xXX  =   .XXx X.  X   XX             ")
	print("     ;XXXXx    ;X X  X  :X.             ")
	print("         .xxXXXXX  XX   XX:xx;          ")
	print("               XXX       ;XXXXXX.       ")
	print("             XXXXx            XX        ")
	print("           xX      :XXXXx.    xX        ")
	print("           ;XX   ;XX;  :XXXXXXX;        ")
	print("             XXXXX;                     ")
	print("                                        ")
	print("                                        ")	





#frequency_run()
time_run()














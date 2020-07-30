import numpy as np
import matplotlib.pyplot as plt
import time
from astropy import units as u
import sp_module as sp
from astropy.io import fits
import math
import numpy
import sys
import scipy.ndimage
import re
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
from astropy.table import Table



def propermotioncorrector():
	
	#gaia values, get these values from gaia website: http://gaia.ari.uni-heidelberg.de/singlesource.html, and plug them in at
	gaia_ra_list = [155.967154663, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166,244.2729637317166,244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 244.2729637317166, 167.964677326, 167.964677326, 336.995585368, 154.898887344, 202.694752018, 202.694752018, 357.420661272, 357.420661272, 357.420661272, 177.730486605, 218.380392642, 169.879477025, 177.779128394, 206.506241420, 113.655009889, 147.066983127, 166.352920955, 156.056289291, 256.866412680, 187.261081061, 135.924782345, 163.422229118]
	gaia_dec_list = [43.892552057, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 55.26724796887516, 33.536954552, 33.536954552, 57.697009626, 19.869815836, 24.23262253, 24.23262253, 36.425075968, 36.425075968, 36.425075968, 48.373233442, 34.296337410, 46.692694805, 35.273104447, 42.018804661, 31.869076649, 51.247355392, 43.525774120, 39.043043423, 64.625481119, 41.729454287, 34.805205380, 52.884641569] 
	gaia_pmra_list = [178.47, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, 94.55898207375307, -176.10, -176.10, -3.07, -498.60, -51.97, -51.97, -0.86, -0.86, -0.86, -1545.70, -90.54, 308.23, -273.16, 11.16, -201.49, -39.17, -4339.89, -92.15, -124.56, -187.45, 133.59, 26.75]
	gaia_pmdec_list = [3.81, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, -433.4149307694429, 116.43, 116.43, -2.37, -43.68, -22.26, -22.26, -47.40, -47.40, -47.40, -962.82, -35.88, -608.70, 253.36, -14.71, -97.10, -88.40, 960.78, 100.60, 224.21, -216.05, -218.24, 40.03]
	gaia_parallax_list = [54.93, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 49.352056624799985, 74.86, 74.86, 0.24, 201.37, 4.61, 4.61, 7.17, 7.17, 7.17, 124.41, 20.90, 95.32, 114.14, 4.37, 66.23, 27.65, 204.06, 55.97, 36.05, 42.44, 27.05, 22.13]
	gaia_ref_epoch_list = [2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5, 2015.5]
	
	

	obs_date_list = ["2015-06-17 20:11:23.0", "2014-05-19 19:49:27.0", "2014-05-20 19:46:31.0", "2014-05-22 19:30:08.0", "2014-05-26 19:30:08.0", "2014-06-02 19:30:08.0","2014-06-03 19:30:08.0", "2014-06-05 19:30:08.0", "2014-06-10 19:50:08.0", "2014-06-12 19:50:08.0", "2014-06-27 20:06:06.0", "2014-07-07 00:33:22.8", "2015-06-07 20:11:08.0", "2015-06-12 20:11:08.0", "2015-06-17 20:11:23.0", "2015-06-19 17:58:08.0", "2015-06-26 20:11:08.0", "2015-06-29 20:11:08.0", "2015-07-01 20:11:08.0", "2015-08-07 18:11:08.0", "2015-08-22 16:11:08.0", "2015-08-21 16:11:08.0", "2019-02-13 21:08:18.0", "2019-09-23 06:00:09.0", "2019-02-07 08:42:33.0", "2019-07-13 11:11:08.0", "2019-04-09 20:11:08.0", "2019-04-04 21:00:08.0", "2016-09-15 19:52:48.0", "2016-10-04 18:23:55.0", "2019-02-28 08:43:56.0", "2014-06-15 15:00:08.0", "2015-10-10 08:46:30.0", "2014-06-08 16:00:08.0", "2017-09-21 08:00:08.0", "2015-05-07 17:56:02.0", "2017-10-15 03:00:08.0", "2018-08-18 07:13:08.0", "2016-03-22 18:41:32.0", "2015-07-08 11:05:42.0", "2018-05-02 22:00:08.0", "2015-01-21 22:58:08.0", "2018-06-11 10:58:15.0", "2015-02-08 20:11:08.0"]
	starname_list = ["LP212-62", "CR_Dra0","CR_Dra1", "CR_Dra2", "CR_Dra3", "CR_Dra4", "CR_Dra5", "CR_Dra6", "CR_Dra7", "CR_Dra8", "CR_Dra9", "CR_Dra10", "CR_Dra11", "CR_Dra12", "CR_Dra13", "CR_Dra14", "CR_Dra15", "CR_Dra16", "CR_Dra17", "CR_Dra18", "CR_Dra19", "CR_Dra20", "CW_Uma1", "CW_Uma2", "DO_Cep", "AD_Leo", "FK_Com1", "FK_Com2", "OU_And1", "OU_And2", "OU_And3", "GJ1151", "2MASS_J14333139+3417472", "LP169-22", "GJ450", "BD+42_2437", "YY_Gem", "2MASS_J09481615+511451", "GJ412", "HAT_182-00605", "G_240-45", "GJ3729", "LP259-39", "J1054129+5253040"]
	starfield_list = ["P156+42", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "P166+35", "P167+32", "P334+58", "P155+19", "P201+25", "P204+25", "P000+36", "P356+36", "P357+38", "P15Hetdex13", "P217+35", "P9Hetdex01", "P169+35", "P205+42", "P115+32", "P147+52", "P167+42", "P154+40", "P260+65", "P184+42", "P135+34", "P164+55"]
	npixel_list = [20197, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 4939, 6861, 6861, 5761, 1581, 7000, 7000, 2855, 2855, 2855, 4934, 4036, 4272, 4238, 3502, 3003, 3593, 3521, 2808, 3149, 4638, 3286, 3543]
	freq_chan_list = [120, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120]
	freq_chanwidth_list = [0.4, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
	time_chan_list = [1798, 1789, 1789, 1789, 1789, 1789, 1789, 1789, 1789, 1789, 1789, 1789, 1723, 1723, 1723, 1723, 1723, 1723, 1723, 1723, 1723, 1723, 1798, 1798, 1798, 1610, 1798, 1632, 1798, 1798, 1798, 1789, 1792, 1512, 1798, 1948, 1620, 1798, 1798, 1798, 1798, 1948, 1798, 1948]
	time_chanwidth_list = [16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0] 
	chan_low_list = [120, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115,  115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120,120, 120, 120, 120, 120]
	chan_up_list = [168, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168]
	msfile_list = ['LP212_P156+42_Subtracted_phaseshifted.ms', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive1.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive2.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive3.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive4.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive5.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive6.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive7.calibrated',  'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive8.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive9.calibrated',  'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive10.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive11.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive12.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive13.calibrated',   'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive14.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive15.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive16.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive17.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive18.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive19.calibrated', 'ELAIS-N1_CR_DRa.dysco.sub.shift.avg.weights.ms.archive20.calibrated', 'P166+35_CW_UMa.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P167+32_CW_UMa.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P334+58_DO_Cep.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P155+19_AD_Leo.dysco.sub.shift.avg.weights.ms.archive.calibrated', 
'P201+25_FK_Com.dysco.sub.shift.avg.weights.ms.archive.calibrated', 
'P204+25_FK_Com.dysco.sub.shift.avg.weights.ms.archive.calibrated', 
'P000+36_OU_And.dysco.sub.shift.avg.weights.ms.archive.calibrated',
'P356+36_OU_And.dysco.sub.shift.avg.weights.ms.archive.calibrated',
'P357+38_OU_And.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P15Hetdex13_GJ1151.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P217+35_2MASS_J14333139+3417472.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P9Hetdex01_LP_169-22.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P179+35_GJ450.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P205+42_BD+42_2437.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P115+32_YY_Gem.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P147+52_2MASS_J09481615+5114518.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P167+42_GJ412.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P154+40_HAT_182-00605.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P260+65_G240-45.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P184+42_GJ3729.dysco.sub.shift.avg.weights.ms.archive.calibrated', 'P135+34_LP_259-39.dysco.sub.shift.avg.weights.ms.archive0.calibrated', 'P164+55_2MASS_J1054129+5253040.dysco.sub.shift.avg.weights.ms.archive0.calibrated']




	lofar_ra_list = []
	lofar_dec_list = []
	for i in range(len(obs_date_list)):

		#lofar obersvation date:
		date_obs = obs_date_list[i]
		
		c = SkyCoord(ra=gaia_ra_list[i] * u.deg,
		       	dec=gaia_dec_list[i] * u.deg,
		        distance=Distance(parallax=gaia_parallax_list[i] * u.mas,allow_negative=False),
		        pm_ra_cosdec=gaia_pmra_list[i] * u.mas/u.yr,
		        pm_dec=gaia_pmdec_list[i] * u.mas/u.yr,
		        obstime=Time(gaia_ref_epoch_list[i], format='decimalyear'))
		print('hier onder komt c')
		print(c)
		
		epoch_lotss = Time(date_obs, format='iso')

		c_gaia_to_lotss_epoch = c.apply_space_motion(epoch_lotss)
		print('hier onder komt c_gaia_to_lotss_epoch')
		print(c_gaia_to_lotss_epoch)


		newra= c_gaia_to_lotss_epoch.ra.value
		newdec = c_gaia_to_lotss_epoch.dec.value
		print('newra:',newra)
		print('newdec:',newdec)
		lofar_ra_list.append(newra)
		lofar_dec_list.append(newdec) 
	
	print(len(starname_list))
	print(lofar_ra_list)
	print(lofar_dec_list)
	print(len(obs_date_list))

	t = Table([starname_list, starfield_list, lofar_ra_list, lofar_dec_list, obs_date_list, gaia_ra_list, gaia_dec_list, gaia_pmra_list, gaia_pmdec_list, gaia_parallax_list, gaia_ref_epoch_list, npixel_list, freq_chan_list, freq_chanwidth_list, time_chan_list, time_chanwidth_list, msfile_list, chan_low_list, chan_up_list], names=('star', 'field', 'RA_gaia_to_lotss_epoch', 'DEC_gaia_to_lotss_epoch', 'lotss_obs_date', 'RA_gaia', 'DEC_gaia', 'pmra_gaia', 'pmdec_gaia', 'parallax_gaia', 'gaia_obs_date', 'npixel_for_rms_region', 'nr_freq_chans_for_msfile', 'MHz_per_chan', 'nr_time_chans_for_msfile', 'sec_per_chan', 'msfile', 'bandwidth_low', 'bandwidth_up'))
	t.write('WABIFAT_input_data_16july.fits', format='fits', overwrite='True')
	print(t)

#propermotioncorrector()	



































############################################################################ WABIFAT RESULTS #######################################################################################################

def yerr_bar_fixer(fluxarray, errfluxarray, boundary):
	errfluxlist = [] 
	for i in range(len(fluxarray)):
		if errfluxarray[i] > boundary:
			errfluxlist.append(0.1*abs(fluxarray[i]))
		else:	
			errfluxlist.append(errfluxarray[i])
	errfluxarray = np.asarray(errfluxlist)
	return errfluxarray
	
def negative_stokes_I_deleter(fluxarray, timearray, errfluxarray, final_slice_radius, rmsarray):
	time_list = []
	flux_list = []
	final_slice_radius_list = []
	errflux_list = []
	rms_list = []
	for i in range(len(fluxarray)):
		if fluxarray[i] >= 0:
			time_list.append(timearray[i])
			flux_list.append(fluxarray[i])
			errflux_list.append(errfluxarray[i])
			final_slice_radius_list.append(final_slice_radius[i])
			rms_list.append(rmsarray[i])
	return np.asarray(time_list), np.asarray(flux_list), np.asarray(errflux_list), np.asarray(final_slice_radius_list), np.asarray(rms_list)


def solo_plot():
	ForT = "TIME"
	starname = "GJ3729"#"2MASS_J14333139+3417472" #Don't use spaces
	starfield = "P184+42" 
	n_iter = 10000
	seedclip = 4.0
	autothreshold = 2
	SNR_accept = 5
	stokes = 'V'
	path = "/home/kas/Documents/MRP1/results/"+starname+"/time_stokesV_15july/results/"
	
	#FREQ
	startMHz = 120 #[MHz], left boundary of measurement set in MHz.
	endMHz = 168 #[MHz], right boundary of measurement set in MHz.
	totalchans = 120 #Total amount of channels spaced between startMHz and endMHz in the measurement set. Usually this should be 120, but for e.g. CR Dra it is 200. 
	chanMHz = 0.4 #[MHz], frequency width of one channel. (Can be calculated by: (endMHz-startMHz)/totalchans = chanMHz, but this went wrong for some reason so just type 0.31 or 0.4 yourself (0.39). 
	lower_chan = 0 #First channel you want to include in the frequency vs flux plot. 
	upper_chan = 120 #Last channel is exclusive because "np.arange" takes out 1.
	freqx0=np.linspace(startMHz+(lower_chan*chanMHz),startMHz+((upper_chan-1)*chanMHz)+(chanMHz/2.0),2)
	freqy0=[0,0]

	
	#TIME
	
	totalscans = 172300 #Total amount of scans spaced between startsec and endsec in the measurement set. Usually this should be 1798, but for e.g. CR Dra it is 1789. 1723 for CR_Dra10to20 i think.
	scansec = 16.0/3600.0 #[sec], Time in one scan in sec. (Can be calculated by: (endsec-startsec)/totalscans = scansec, but it should for now allways be 16s. 
	lower_scan = 0 #First scan you want to include in the time vs flux plot. 
	upper_scan = 172300 #Last scan is exclusive because "np.arange" takes out 1.
	startsec = 0.0/3600.0 #[sec], left boundary of measurement set in sec.
	endsec = scansec*totalscans#CRDRA0:28624  #LP212-62:28768 #[sec], right boundary of measurement set in sec.
	


	if ForT == "FREQ":
		#stokes I
		freqarray = np.loadtxt(path+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		fluxarray = np.loadtxt(path+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		rmsarray = np.loadtxt(path+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		errfluxarray = np.loadtxt(path+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_centralfreqarray = np.loadtxt(path+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(path+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		"""
		#This part will fix that the time is not plotted in Julian Date, but just from 0 to 8 hours. 
		if freqarray[0] >= misfit_centralfreqarray[0]:
			lowest_freq = misfit_centralfreqarray[0]-final_misfit_radius[0]
		else: 
			lowest_freq = freqarray[0]-final_slice_radius[0]
		
		if freqarray[-1] >= misfit_centralfreqarray[-1]:
			highest_freq = freqarray[-1]+final_slice_radius[-1]
		else: 
			highest_freq = misfit_centralfreqarray[-1]+final_misfit_radius[-1]
		freqx0 = [lowest_freq, highest_freq] 
		freqy0=[0,0]
		"""
		
		#random arguments: uplims=True, lolims=True  ####str(len(fluxarray))+
		plt.figure()
		plt.plot(freqx0,freqy0,c='green')
		#plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = '1 x '+str(SNR_accept)+r'$\sigma$ detection', fmt='.') 
		plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(freqarray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		#plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*-1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.title('Frequency VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.title('Star: '+starname+', field: '+starfield)
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.savefig(path+starname+'_'+ForT+'_pol'+stokes+'.eps',format='eps')
		plt.show()		
	

	if ForT == "TIME": 
		timearray = np.loadtxt(path+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		fluxarray = np.loadtxt(path+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		rmsarray = np.loadtxt(path+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		errfluxarray = np.loadtxt(path+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_centraltimearray = np.loadtxt(path+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(path+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		
	
		radiuslist=final_slice_radius.tolist()
		misfitradiuslist = final_misfit_radius.tolist()
		#print("sum of radii x 2 =", sum(radiuslist)*2.0 + sum(misfitradiuslist)*2.0)
		#This part will replace the huge y_error bars (above boundary) with 0.1*peakflux
		#errfluxarray = yerr_bar_fixer(fluxarray, errfluxarray, 20*10**(-3))
		#timearray, fluxarray, errfluxarray, final_slice_radius, rmsarray = negative_stokes_I_deleter(fluxarray, timearray, errfluxarray, final_slice_radius, rmsarray)
	
		timelist = timearray.tolist()
		fluxlist = fluxarray.tolist()
		slicelist = final_slice_radius.tolist()
		rmslist = rmsarray.tolist()
		errfluxlist = errfluxarray.tolist() 
		if isinstance(timelist, float):
			timelist = [timelist]
			fluxlist = [fluxlist]
			slicelist = [slicelist]
			rmslist = [rmslist]
			errfluxlist = [errfluxlist]
		timearray = np.asarray(timelist)
		fluxarray = np.asarray(fluxlist)
		final_slice_radius = np.asarray(slicelist)
		rmsarray = np.asarray(rmslist)
		errfluxarray = np.asarray(errfluxlist)

		misfit_timelist = misfit_centraltimearray.tolist()
		misfit_rmslist = misfit_rmsarray.tolist()
		misfit_radiuslist = final_misfit_radius.tolist()
		if isinstance(misfit_timelist, float):
			misfit_timelist = [misfit_timelist]
			misfit_rmslist = [misfit_rmslist]
			misfit_radiuslist = [misfit_radiuslist]
		misfit_centraltimearray = np.asarray(misfit_timelist)
		misfit_rmsarray = np.asarray(misfit_rmslist)
		final_misfit_radius = np.asarray(misfit_radiuslist)
		print(misfit_centraltimearray)
		print(misfit_rmsarray)
		print(final_misfit_radius)
		#print("misfit_timearray[0]:", misfit_timelist[0])
		#print("misfit_radius[0]:", final_misfit_radius[0])
		
		print(len(timearray))
		if len(timearray) == 0:
			
			lowest_time = misfit_centraltimearray[0]-final_misfit_radius[0]
			highest_time = misfit_centraltimearray[-1]+final_misfit_radius[-1]
		else: 
			#This part will fix that the time is not plotted in Julian Date, but just from 0 to 8 hours. 
			if timearray[0] >= misfit_centraltimearray[0]:
				lowest_time = misfit_centraltimearray[0]-final_misfit_radius[0]
			else: 
				lowest_time = timearray[0]-final_slice_radius[0]
			
			if timearray[-1] >= misfit_centraltimearray[-1]:
				highest_time = timearray[-1]+final_slice_radius[-1]
			else: 
				highest_time = misfit_centraltimearray[-1]+final_misfit_radius[-1]
		timex0 = [lowest_time-lowest_time, highest_time-lowest_time] 
		timey0=[0,0]
		
		#timex0=[0, 2*3.36]
		#timey0 = [0,0]
	
		plt.figure()
		plt.plot(timex0,timey0,c='green')
		plt.errorbar(timearray-lowest_time, fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(timearray-lowest_time, rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray-lowest_time, misfit_rmsarray*-1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		plt.title('Star: '+starname+', field: '+starfield)
		#plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.savefig(path+starname+'_'+ForT+'_pol'+stokes+'.eps',format='eps')
		plt.show()		

solo_plot()






def duo_plot():
	ForT = "FREQ"
	starname = "BD+42_2437" #Don't use spaces
	starfield = "P205+42" 
	n_iter = 10000
	seedclip = 4.0
	autothreshold = 2
	SNR_accept = 5
	stokesI = 'I'
	stokesV = 'V'
	
	pathI = "/home/kas/Documents/MRP1/results/BD+42_2437/freq_stokesI_22june/results/" 
	pathV = "/home/kas/Documents/MRP1/results/BD+42_2437/freq_stokesV_22june/results/"
	
	#FREQ
	startMHz = 120 #[MHz], left boundary of measurement set in MHz.
	endMHz = 168 #[MHz], right boundary of measurement set in MHz.
	totalchans = 120 #Total amount of channels spaced between startMHz and endMHz in the measurement set. Usually this should be 120, but for e.g. CR Dra it is 200. 
	chanMHz = 0.4 #[MHz], frequency width of one channel. (Can be calculated by: (endMHz-startMHz)/totalchans = chanMHz, but this went wrong for some reason so just type 0.31 or 0.4 yourself (0.39). 
	lower_chan = 0 #First channel you want to include in the frequency vs flux plot. 
	upper_chan = 120 #Last channel is exclusive because "np.arange" takes out 1.
	freqx0=np.linspace(startMHz+(lower_chan*chanMHz),startMHz+((upper_chan-1)*chanMHz)+(chanMHz/2.0),2)
	freqy0=[0,0]

	
	#TIME
	
	totalscans = 1798 #Total amount of scans spaced between startsec and endsec in the measurement set. Usually this should be 1798, but for e.g. CR Dra it is 1789. 1723 for CR_Dra10to20 i think.
	scansec = 16.0/3600.0 #[sec], Time in one scan in sec. (Can be calculated by: (endsec-startsec)/totalscans = scansec, but it should for now allways be 16s. 
	lower_scan = 0 #First scan you want to include in the time vs flux plot. 
	upper_scan = 1789 #Last scan is exclusive because "np.arange" takes out 1.
	startsec = 0.0/3600.0 #[sec], left boundary of measurement set in sec.
	endsec = scansec*totalscans#CRDRA0:28624  #LP212-62:28768 #[sec], right boundary of measurement set in sec.
	


	if ForT == "FREQ":
		#stokes I:
		freqarrayI = np.loadtxt(pathI+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		fluxarrayI = np.loadtxt(pathI+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_slice_radiusI = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		rmsarrayI = np.loadtxt(pathI+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		errfluxarrayI = np.loadtxt(pathI+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_centralfreqarrayI = np.loadtxt(pathI+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_rmsarrayI = np.loadtxt(pathI+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_misfit_radiusI = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	
		#stokes V:
		freqarrayV = np.loadtxt(pathV+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		fluxarrayV = np.loadtxt(pathV+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		rmsarrayV = np.loadtxt(pathV+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		errfluxarrayV = np.loadtxt(pathV+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)

		misfit_centralfreqarrayV = np.loadtxt(pathV+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfit_rmsarrayV = np.loadtxt(pathV+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	
		
		#This part will replace the huge y_error bars (above boundary in Jy) with 0.1*peakflux
		#errfluxarrayI = yerr_bar_fixer(fluxarrayI, errfluxarrayI, 10*10**(-3))
		#errfluxarrayV = yerr_bar_fixer(fluxarrayV, errfluxarrayV, 10*10**(-3))
		
		#Random arguments: uplims=True, lolims=True  ####str(len(fluxarray))+
		plt.figure()
		plt.plot(freqx0,freqy0,c='green')
		#stokes I:
		plt.errorbar(freqarrayI,fluxarrayI*1000.0, xerr = final_slice_radiusI, yerr = errfluxarrayI*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarrayI))+'x '+str(SNR_accept)+r'$\sigma$ detections stokes I', fmt='.') 
		plt.errorbar(freqarrayI,rmsarrayI*1000.0*SNR_accept, xerr = final_slice_radiusI, yerr=rmsarrayI*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS stokes I',fmt='.') 
		plt.errorbar(misfit_centralfreqarrayI,misfit_rmsarrayI*1000.0*SNR_accept, xerr = final_misfit_radiusI, yerr=misfit_rmsarrayI*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection stokes I, '+str(SNR_accept)+'*RMS',fmt='.') 
		#stokes V:
		plt.errorbar(freqarrayV,fluxarrayV*1000.0, xerr = final_slice_radiusV, yerr = errfluxarrayV*1000, solid_capstyle='projecting', capsize=4, c = 'magenta', label = str(len(fluxarrayV))+'x '+str(SNR_accept)+r'$\sigma$ detections stokes V', fmt='.') 
		plt.errorbar(freqarrayV,rmsarrayV*1000.0*SNR_accept, xerr = final_slice_radiusV, yerr=rmsarrayV*100, c = 'pink', label = 'detection, '+str(SNR_accept)+'*local RMS stokes V',fmt='.') 
		plt.errorbar(misfit_centralfreqarrayV,misfit_rmsarrayV*-1000.0*SNR_accept, xerr = final_misfit_radiusV, yerr=misfit_rmsarrayV*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'black', label = 'non-detection stokes V, '+str(SNR_accept)+'*RMS',fmt='.') 

		plt.title('Frequency VS Fluxdensity in Stokes I and V for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		
	

	if ForT == "TIME": 
		#stokes I:
		timearrayI = np.loadtxt(pathI+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		fluxarrayI = np.loadtxt(pathI+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_slice_radiusI = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		rmsarrayI = np.loadtxt(pathI+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		errfluxarrayI = np.loadtxt(pathI+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_centraltimearrayI = np.loadtxt(pathI+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_rmsarrayI = np.loadtxt(pathI+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_misfit_radiusI = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	
		#stokes V:
		timearrayV = np.loadtxt(pathV+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		fluxarrayV = np.loadtxt(pathV+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		rmsarrayV = np.loadtxt(pathV+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		errfluxarrayV = np.loadtxt(pathV+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfit_centraltimearrayV = np.loadtxt(pathV+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfit_rmsarrayV = np.loadtxt(pathV+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	

		
		#This part will replace the huge y_error bars (above boundary in Jy) with 0.1*peakflux
		errfluxarrayI = yerr_bar_fixer(fluxarrayI, errfluxarrayI, 10*10**(-3))
		errfluxarrayV = yerr_bar_fixer(fluxarrayV, errfluxarrayV, 10*10**(-3))
		timearrayI, fluxarrayI, errfluxarrayI, final_slice_radiusI, rmsarrayI = negative_stokes_I_deleter(fluxarrayI, timearrayI, errfluxarrayI, final_slice_radiusI, rmsarrayI)


		#This part will set your lower and upper time limit, (which is the first time value and the last (plus and minus their channelwidths))
		if timearrayI[0] >= misfit_centraltimearrayI[0]:
			lowest_timeI = misfit_centraltimearrayI[0]-final_misfit_radiusI[0]
		else: 
			lowest_timeI = timearrayI[0]-final_slice_radiusI[0]
		
		if timearrayI[-1] >= misfit_centraltimearrayI[-1]:
			highest_timeI = timearrayI[-1]+final_slice_radiusI[-1]
		else: 
			highest_timeI = misfit_centraltimearrayI[-1]+final_misfit_radiusI[-1]

		if timearrayV[0] >= misfit_centraltimearrayV[0]:
			lowest_timeV = misfit_centraltimearrayV[0]-final_misfit_radiusV[0]
		else: 
			lowest_timeV = timearrayV[0]-final_slice_radiusV[0] 
		
		if timearrayV[-1] >= misfit_centraltimearrayV[-1]:
			highest_timeV = timearrayV[-1]+final_slice_radiusV[-1]
		else: 
			highest_timeV = misfit_centraltimearrayV[-1]+final_misfit_radiusV[-1]


		if lowest_timeI >= lowest_timeV:
			lowest_time = lowest_timeV
		else:
			lowest_time = lowest_timeI
		if highest_timeI >= highest_timeV:
			highest_time = highest_timeI
		else:
			highest_time = highest_timeV 

		timex0 = [lowest_time-lowest_time, highest_time-lowest_time] 
		timey0=[0,0]

	
		plt.figure()
		plt.plot(timex0,timey0,c='green')
		#stokes I:
		plt.errorbar(timearrayI-lowest_time,fluxarrayI*1000.0, xerr = final_slice_radiusI, yerr = errfluxarrayI*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(SNR_accept)+r'$\sigma$ detections (I)', fmt='.') 
		plt.errorbar(timearrayI-lowest_time,rmsarrayI*1000.0*SNR_accept, xerr = final_slice_radiusI, yerr=rmsarrayI*100, c = 'aqua', label = 'detections, '+str(SNR_accept)+'x local RMS  (I)',fmt='.') 
		plt.errorbar(misfit_centraltimearrayI-lowest_time,misfit_rmsarrayI*1000.0*SNR_accept, xerr = final_misfit_radiusI, yerr=misfit_rmsarrayI*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detections, '+str(SNR_accept)+'x RMS  (I)',fmt='.') 
		#stokes V:
		plt.errorbar(timearrayV-lowest_time,fluxarrayV*1000.0, xerr = final_slice_radiusV, yerr = errfluxarrayV*1000, solid_capstyle='projecting', capsize=4, c = 'magenta', label = str(SNR_accept)+r'$\sigma$ detections (V)', fmt='.') 
		plt.errorbar(timearrayV-lowest_time,rmsarrayV*1000.0*SNR_accept, xerr = final_slice_radiusV, yerr=rmsarrayV*100, c = 'pink', label = 'detections, '+str(SNR_accept)+'x local RMS (V)',fmt='.') 
		plt.errorbar(misfit_centraltimearrayV-lowest_time,misfit_rmsarrayV*-1000.0*SNR_accept, xerr = final_misfit_radiusV, yerr=misfit_rmsarrayV*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'black', label = 'non-detections, '+str(SNR_accept)+'x RMS (V)',fmt='.') 

		plt.text(0.32, 28.3, 'Stokes I', color = 'blue', fontsize = 20)
		plt.text(0.32, -25.9, 'Stokes V', color = 'magenta', fontsize = 20)

		plt.title('Time VS Fluxdensity in Stokes I and V for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		




#duo_plot()
####################################################################################################################################################################################################



















############################################################################ Forced Fit Results #######################################################################################################
def forced_plot():
	ForT = "TIME"
	starname = "CR_Dra0" #Don't use spaces
	starname_forced = "CR_Dra0_time_stokesI_31may"
	starfield = "ELAIS-N1" 
	n_iter = 10000
	seedclip = 4.0
	autothreshold = 2
	SNR_accept = 5
	stokes = 'I'
	path = "/home/kas/Documents/MRP1/results/CR_Dra0/time_stokesI_31may/results/"
	forced_fitter_path = '/home/kas/Documents/MRP1/random/forced_fitter/'
	
	#Frequency Run:
	""" 
	if ForT == "FREQ":
		#stokes I
		freqarray = np.loadtxt(path+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		fluxarray = np.loadtxt(path+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		rmsarray = np.loadtxt(path+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		errfluxarray = np.loadtxt(path+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_centralfreqarray = np.loadtxt(path+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(path+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		
		#random arguments: uplims=True, lolims=True  ####str(len(fluxarray))+
		plt.figure()
		plt.plot(freqx0,freqy0,c='green')
		#plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = '1 x '+str(SNR_accept)+r'$\sigma$ detection', fmt='.') 
		plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections stokes I', fmt='.') 
		plt.errorbar(freqarray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		plt.title('Frequency VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.legend(loc='center right')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		
	"""
	#Time run:
	if ForT == "TIME": 
		#reads in the three files made by the forced_fitter definition.
		fluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/"+starname_forced+"_forced_peak_flux.txt", unpack = True)
		rmsarray = np.loadtxt(forced_fitter_path+starname_forced+"/"+starname_forced+"_forced_local_rms.txt" ,unpack=True)
		errfluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/"+starname_forced+"_forced_err_peak_flux.txt",unpack=True)
		print(len(fluxarray))
		print(len(rmsarray))
		print(len(errfluxarray))

		#reads in the files made by WABIFAT, for the time and channelwidths, these of course are not changed by the forced_fitter. 
		timearray = np.loadtxt(path+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_centraltimearray = np.loadtxt(path+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		#this is a slightly overkill way to sort add the times and slices of the detections and the nondetections together, and consequently sort them such that the times and slices correspond to the fluxes, errors, and local_rms's of the forced fitter. 
		det_time = timearray.tolist()
		nondet_time = misfit_centraltimearray.tolist()
		det_slice = final_slice_radius.tolist()
		nondet_slice = final_misfit_radius.tolist()
		totaltime_not_sorted = det_time+nondet_time
		totalslice_not_sorted = det_slice+nondet_slice

		for z in range(len(totaltime_not_sorted)): 
			swap = z + np.argmin(totaltime_not_sorted[z:])
			(totaltime_not_sorted[z], totaltime_not_sorted[swap]) = (totaltime_not_sorted[swap], totaltime_not_sorted[z])
			(totalslice_not_sorted[z], totalslice_not_sorted[swap]) = (totalslice_not_sorted[swap], totalslice_not_sorted[z])

		totaltime_sorted = np.asarray(totaltime_not_sorted)
		totalslice_sorted = np.asarray(totalslice_not_sorted)

		#determining the lowest time value and highest time value such that we can plot the green y=0 line I always plot. 
		lowest_time = totaltime_sorted[0]-totalslice_sorted[0]
		highest_time = totaltime_sorted[-1]+totalslice_sorted[-1]
		timex0 = [lowest_time, highest_time] 
		timey0=[0,0]

		#plot
		plt.figure()
		plt.plot(timex0,timey0,c='green')
		plt.errorbar(totaltime_sorted-lowest_time, fluxarray*1000.0, xerr = totalslice_sorted, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(totaltime_sorted-lowest_time, rmsarray*1000.0, xerr = totalslice_sorted, c = 'aqua', label = 'detection, local RMS',fmt='.') 
		#plt.errorbar(misfit_centraltimearray-lowest_time, misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		

#forced_plot() 




#########################################################################################################################################################################################################
















def polfrac_plot():
	ForT = "FREQ"
	starname = "LP212-62" #Don't use spaces
	starfield = "P156+42" 
	n_iter = 100
	seedclip = 4.0
	autothreshold = 2
	SNR_accept = 5
	stokesI = 'I'
	stokesV = 'V'
	pathI = "/home/kas/Documents/MRP1/results/LP212-62/freq_stokesI/results/"
	pathV = "/home/kas/Documents/MRP1/resultsLP212-62/results_freq_3april/" 
	
	#FREQ
	startMHz = 120 #[MHz], left boundary of measurement set in MHz.
	endMHz = 168 #[MHz], right boundary of measurement set in MHz.
	totalchans = 120 #Total amount of channels spaced between startMHz and endMHz in the measurement set. Usually this should be 120, but for e.g. CR Dra it is 200. 
	chanMHz = 0.4 #[MHz], frequency width of one channel. (Can be calculated by: (endMHz-startMHz)/totalchans = chanMHz, but this went wrong for some reason so just type 0.31 or 0.4 yourself (0.39). 
	lower_chan = 0 #First channel you want to include in the frequency vs flux plot. 
	upper_chan = 120 #Last channel is exclusive because "np.arange" takes out 1.
	freqx0=np.linspace(startMHz+(lower_chan*chanMHz),startMHz+((upper_chan-1)*chanMHz)+(chanMHz/2.0),2)
	freqy0=[0,0]

	
	#TIME
	startsec = 0.0/3600.0 #[sec], left boundary of measurement set in sec.
	endsec = 28768.0/3600.0#CRDRA0:28624  #LP212-62:28768 #[sec], right boundary of measurement set in sec.
	totalscans = 1798 #Total amount of scans spaced between startsec and endsec in the measurement set. Usually this should be 1798, but for e.g. CR Dra it is 1789. 1723 for CR_Dra10to20 i think.
	scansec = 16.0/3600.0 #[sec], Time in one scan in sec. (Can be calculated by: (endsec-startsec)/totalscans = scansec, but it should for now allways be 16s. 
	lower_scan = 0 #First scan you want to include in the time vs flux plot. 
	upper_scan = 1798 #Last scan is exclusive because "np.arange" takes out 1.
	timex0=np.linspace(startsec+(lower_scan*scansec),startsec+((upper_scan-1)*scansec)+(scansec/2.0),2)  #+(scansec/2.0)
	timey0=[0,0]


	if ForT == "FREQ":
		#stokes I
		freqarray = np.loadtxt(pathI+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		fluxarray = np.loadtxt(pathI+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_slice_radius = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		rmsarray = np.loadtxt(pathI+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		errfluxarray = np.loadtxt(pathI+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)

		misfit_centralfreqarray = np.loadtxt(pathI+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(pathI+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	


		#stokes V
		freqarrayV = np.loadtxt(pathV+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		fluxarrayV = np.loadtxt(pathV+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		rmsarrayV = np.loadtxt(pathV+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		errfluxarrayV = np.loadtxt(pathV+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)

		misfit_centralfreqarrayV = np.loadtxt(pathV+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfit_rmsarrayV = np.loadtxt(pathV+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	

		misfits_freqV = np.loadtxt(pathV+ForT+"_misfits_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfits_zeroV = np.loadtxt(pathV+ForT+"_misfits_zero_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		

		frac_fluxlist = []
		frac_freqlist = []
		frac_radiuslist = []
		for i in range(len(freqarray)):
			for j in range(len(freqarrayV)):
				if freqarray[i] == freqarrayV[j]:
					frac_fluxlist.append(abs(fluxarrayV[j]/fluxarray[i]))
					frac_freqlist.append(freqarray[i])
					if final_slice_radius[i] >= final_slice_radiusV[j]:
						frac_radiuslist.append(final_slice_radiusV[j])
					else:
						frac_radiuslist.append(final_slice_radius[i])
				#else: 
				#	not_same_freqlist.append(freqarray[i])
				#	not_same_freq
				
		frac_fluxarray = np.asarray(frac_fluxlist)
		frac_freqarray = np.asarray(frac_freqlist)
		frac_radiusarray = np.asarray(frac_radiuslist)
		



		#, uplims=True, lolims=True
		plt.figure()
		plt.plot(freqx0,freqy0,c='green')
		plt.errorbar(frac_freqarray,frac_fluxarray, xerr = frac_radiusarray, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(frac_fluxarray))+"x stokes(V/I)", fmt='.') 
		#plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection stokes I, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_freq, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_freq))+' misfit channels')
		
		

		#plt.title('Frequency VS Fluxdensity in Stokes '+stokesI+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.title('Stokes(V/I), LP212-62, SC4, AT2, SNRA5')
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('abs(V/I)')

		plt.show()		
	

	if ForT == "TIME": 
		timearray = np.loadtxt(path+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		fluxarray = np.loadtxt(path+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		rmsarray = np.loadtxt(path+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		errfluxarray = np.loadtxt(path+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)

		misfit_centraltimearray = np.loadtxt(path+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(path+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		misfits_time = np.loadtxt(path+ForT+"_misfits_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfits_zero = np.loadtxt(path+ForT+"_misfits_zero_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		
		# This part will fix the SNR acceptance of 5sigma. It also makes the huge errorbars smaller if above a certain boundary.
		print(len(fluxarray))
		print(len(timearray))
		print(len(rmsarray))
		print(len(final_slice_radius))
		print(len(errfluxarray))


		fluxlist = []
		rmslist = [] 
		timelist = []
		radiuslist = [] 
		errfluxlist = [] 
		for i in range(len(fluxarray)):
			if abs(fluxarray[i])/abs(rmsarray[i]) >= 5:
				fluxlist.append(fluxarray[i])
				print("fluxlist", fluxlist)
				rmslist.append(rmsarray[i])
				timelist.append(timearray[i])
				print("timelist", timelist)
				radiuslist.append(final_slice_radius[i])
				if errfluxarray[i] > 10*10**(-3):
					errfluxlist.append(0.1*fluxarray[i])
				else:	
					errfluxlist.append(errfluxarray[i])
		fluxarray = np.asarray(fluxlist)
		rmsarray = np.asarray(rmslist)
		timearray = np.asarray(timelist)
		final_slice_radius = np.asarray(radiuslist)
		errfluxarray = np.asarray(errfluxlist)

		print(len(fluxarray))
		print(len(timearray))
		print(len(rmsarray))
		print(len(final_slice_radius))
		print(len(errfluxarray))
		

		plt.figure()
		plt.plot(timex0,timey0,c='green')
		plt.errorbar(timearray,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(timearray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray,misfit_rmsarray*-1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_time, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_time))+' misfit channels')
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		




#polfrac_plot()


def hallo_plot():
	ForT = "TIME"
	starname = "CR_Dra" #Don't use spaces
	starfield = "P156+42" 
	n_iter = 100
	seedclip = 4.0
	autothreshold = 2
	SNR_accept = 5
	stokes = 'V'
	
	path = "/home/kas/Documents/MRP1/results/LP212-62/time_stokesV_9april/results/"
	pathV = "/home/kas/Documents/MRP1/resultsLP212-62/results_freq_3april/" 
	
	#FREQ
	startMHz = 115 #[MHz], left boundary of measurement set in MHz.
	endMHz = 177 #[MHz], right boundary of measurement set in MHz.
	totalchans = 160 #Total amount of channels spaced between startMHz and endMHz in the measurement set. Usually this should be 120, but for e.g. CR Dra it is 200. 
	chanMHz = 0.39 #[MHz], frequency width of one channel. (Can be calculated by: (endMHz-startMHz)/totalchans = chanMHz, but this went wrong for some reason so just type 0.31 or 0.4 yourself (0.39). 
	lower_chan = 0 #First channel you want to include in the frequency vs flux plot. 
	upper_chan = 160 #Last channel is exclusive because "np.arange" takes out 1.
	freqx0=np.linspace(startMHz+(lower_chan*chanMHz),startMHz+((upper_chan-1)*chanMHz)+(chanMHz/2.0),2)
	freqy0=[0,0]

	
	#TIME
	
	totalscans = 1789 #Total amount of scans spaced between startsec and endsec in the measurement set. Usually this should be 1798, but for e.g. CR Dra it is 1789. 1723 for CR_Dra10to20 i think.
	scansec = 16.0/3600.0 #[sec], Time in one scan in sec. (Can be calculated by: (endsec-startsec)/totalscans = scansec, but it should for now allways be 16s. 
	lower_scan = 0 #First scan you want to include in the time vs flux plot. 
	upper_scan = 1789 #Last scan is exclusive because "np.arange" takes out 1.
	startsec = 0.0/3600.0 #[sec], left boundary of measurement set in sec.
	endsec = scansec*totalscans#CRDRA0:28624  #LP212-62:28768 #[sec], right boundary of measurement set in sec.
	


	if ForT == "FREQ":
		#stokes I
		freqarray = np.loadtxt(pathI+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		fluxarray = np.loadtxt(pathI+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_slice_radius = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		rmsarray = np.loadtxt(pathI+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		errfluxarray = np.loadtxt(pathI+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)

		misfit_centralfreqarray = np.loadtxt(pathI+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(pathI+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	

		misfits_freq = np.loadtxt(pathI+ForT+"_misfits_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
		misfits_zero = np.loadtxt(pathI+ForT+"_misfits_zero_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)

		#stokes V
		"""
		freqarrayV = np.loadtxt(pathV+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		fluxarrayV = np.loadtxt(pathV+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		rmsarrayV = np.loadtxt(pathV+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		errfluxarrayV = np.loadtxt(pathV+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)

		misfit_centralfreqarrayV = np.loadtxt(pathV+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfit_rmsarrayV = np.loadtxt(pathV+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	

		misfits_freqV = np.loadtxt(pathV+ForT+"_misfits_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		misfits_zeroV = np.loadtxt(pathV+ForT+"_misfits_zero_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
		"""
		




		#, uplims=True, lolims=True  ####str(len(fluxarray))+
		plt.figure()
		plt.plot(freqx0,freqy0,c='green')
		plt.errorbar(freqarray,fluxarray*1000.0, xerr = final_slice_radius*(0.31/0.39), yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections stokes I', fmt='.') 
		plt.errorbar(freqarray,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius*(0.31/0.39), yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS stokes I',fmt='.') 
		plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius*(0.31/0.39), yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection stokes I, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_freq, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_freq))+' misfit channels')
		
		
		#plt.errorbar(freqarrayV,-fluxarrayV*1000.0, xerr = final_slice_radiusV, yerr = errfluxarrayV*1000, solid_capstyle='projecting', capsize=4, c = 'purple', label = str(len(fluxarrayV))+'x '+str(SNR_accept)+r'$\sigma$ detections stokes V', fmt='.') 
		#plt.errorbar(freqarrayV,-rmsarrayV*1000.0*SNR_accept, xerr = final_slice_radiusV, yerr=rmsarrayV*100, c = 'pink', label = 'detection, '+str(SNR_accept)+'*local RMS stokes V',fmt='.') 
		#plt.errorbar(misfit_centralfreqarrayV,misfit_rmsarrayV*1000.0*SNR_accept, xerr = final_misfit_radiusV, yerr=misfit_rmsarrayV*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'yellow', label = 'non-detection stokes V, '+str(SNR_accept)+'*RMS',fmt='.') 

		plt.title('Frequency VS Fluxdensity in Stokes '+stokesI+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		#plt.title('Frequency vs Flux, CR_Dra0, stokes I, SC=5, AT=4, n=10000)
		plt.legend(loc='best')
		plt.xlabel('Frequency (MHz)')
		plt.ylabel('Fluxdensity (mJy/Beam)')

		plt.show()		
	

	if ForT == "TIME": 
		timearray = np.loadtxt(path+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		fluxarray = np.loadtxt(path+ForT+"_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		rmsarray = np.loadtxt(path+ForT+"_rms_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		errfluxarray = np.loadtxt(path+ForT+"_err_flux_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)

		misfit_centraltimearray = np.loadtxt(path+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		misfit_rmsarray = np.loadtxt(path+ForT+"_final_misfit_rms_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

		#misfits_time = np.loadtxt(path+ForT+"_misfits_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		#misfits_zero = np.loadtxt(path+ForT+"_misfits_zero_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
		
		# This part will fix the SNR acceptance of 5sigma. It also makes the huge errorbars smaller if above a certain boundary.
		print(len(fluxarray))
		print(len(timearray))
		print(len(rmsarray))
		print(len(final_slice_radius))
		print(len(errfluxarray))

		
		#fluxlist = []
		#rmslist = [] 
		#timelist = []
		#radiuslist = [] 
		errfluxlist = [] 
		for i in range(len(fluxarray)):
			#if abs(fluxarray[i])/abs(rmsarray[i]) >= 5:
				#fluxlist.append(fluxarray[i])
				#print("fluxlist", fluxlist)
				#rmslist.append(rmsarray[i])
				#timelist.append(timearray[i])
				#print("timelist", timelist)
				#radiuslist.append(final_slice_radius[i])
			if errfluxarray[i] > 10*10**(-3):
				errfluxlist.append(0.1*fluxarray[i])
			else:	
				errfluxlist.append(errfluxarray[i])
		#fluxarray = np.asarray(fluxlist)
		#rmsarray = np.asarray(rmslist)
		##timearray = np.asarray(timelist)
		#final_slice_radius = np.asarray(radiuslist)
		errfluxarray = np.asarray(errfluxlist)

		print(len(fluxarray))
		print(len(timearray))
		print(len(rmsarray))
		print(len(final_slice_radius))
		print(len(errfluxarray))
			
		
		"""
		plt.figure()
		plt.plot(timex0,timey0,c='green')
		plt.errorbar(timearray,fluxarray*1000.0, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(timearray,rmsarray*1000.0*SNR_accept, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray,misfit_rmsarray*1000.0*SNR_accept, yerr=misfit_rmsarray*800, lolims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_time, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_time))+' misfit channels')
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlim(left=misfit_centraltimearray[0], right=misfit_centraltimearray[0]+8)
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()	
		"""	
	
		if timearray[0] >= misfit_centraltimearray[0]:
			lowest_time = misfit_centraltimearray[0]-final_misfit_radius[0]
		else: 
			lowest_time = timearray[0]-final_slice_radius[0] 
		
		if timearray[-1] >= misfit_centraltimearray[-1]:
			highest_time = timearray[-1]+final_slice_radius[-1]
		else: 
			highest_time = misfit_centraltimearray[-1]+final_misfit_radius[-1]

		print(lowest_time)
		print(highest_time)
		timex0 = [lowest_time-lowest_time, highest_time-lowest_time] 
		timey0=[0,0]

	
		plt.figure()
		plt.plot(timex0,timey0,c='green')
		plt.errorbar(timearray-lowest_time,fluxarray*1000.0, xerr = final_slice_radius, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+'x '+str(SNR_accept)+r'$\sigma$ detections', fmt='.') 
		plt.errorbar(timearray-lowest_time,rmsarray*1000.0*SNR_accept, xerr = final_slice_radius, yerr=rmsarray*100, c = 'aqua', label = 'detection, '+str(SNR_accept)+'*local RMS',fmt='.') 
		plt.errorbar(misfit_centraltimearray-lowest_time,misfit_rmsarray*-1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
		#plt.scatter(misfits_time, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_time))+' misfit channels')
		plt.title('Time VS Fluxdensity in Stokes '+stokes+' for '+starname+', field: '+starfield+' n='+str(n_iter)+' AT='+str(autothreshold)+' SC='+str(seedclip))
		plt.legend(loc='best')
		plt.xlabel('Time (hrs)')
		plt.ylabel('Fluxdensity (mJy/Beam)')
		plt.show()		




#hallo_plot()




	   





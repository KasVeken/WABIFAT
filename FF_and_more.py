import os,sys
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
import glob
import matplotlib.pyplot as plt 

#['CR_Dra0', 'CR_Dra0', 'CR_Dra0', 'CR_Dra0', 'CR_Dra4', 'CR_Dra4', 'CR_Dra4', 'CR_Dra4', 'CR_Dra11', 'CR_Dra11', 'CR_Dra11', 'CR_Dra11'], ["FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_9april', 'freq_stokesV_14april', 'time_stokesI_4june', 'time_stokesV_16june' ], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "9", "9", "9", "9", "9", "9", "9", "9", "3", "4" , "9", "9" ], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]

#Forced fitting mechanism for WABIFAT results.


#['CR_Dra0', 'CR_Dra0', 'CR_Dra4', 'CR_Dra4', 'CR_Dra11', 'CR_Dra11'], ["TIME", "TIME", "TIME", "TIME", "TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_7june', 'time_stokesI_31may', 'time_stokesV_7june', 'time_stokesI_4june', 'time_stokesV_16june' ], ["I", "V", "I", "V", "I", "V"], [ "9", "9", "9", "9", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000]


#['CR_Dra1', 'CR_Dra1', 'CR_Dra2', 'CR_Dra2', 'CR_Dra2', 'CR_Dra2', 'CR_Dra3', 'CR_Dra3', 'CR_Dra3', 'CR_Dra3', 'CR_Dra5', 'CR_Dra5', 'CR_Dra5', 'CR_Dra5'], ["TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_10april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april' , 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june'], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "9", "9", "2", "4", "9", "9", "2", "4", "9", "9", "2", "4", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000, 10000, 10000, 10000]

#['CR_Dra6', 'CR_Dra6', 'CR_Dra6', 'CR_Dra6', 'CR_Dra7', 'CR_Dra7', 'CR_Dra7', 'CR_Dra7', 'CR_Dra8', 'CR_Dra8', 'CR_Dra8', 'CR_Dra8'], ["FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_10april', 'freq_stokesV_14april' , 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june'], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "2", "4", "9", "9", "2", "4", "9", "9", "2", "4", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000, 10000, 10000, 10000]

#['LP212-62', 'LP212-62'], ["TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_31may'], ["I", "V"], [ "9", "9"], [1000, 1000]
#['2MASS_J14333139+3417472', '2MASS_J14333139+3417472', '2MASS_J14333139+3417472', '2MASS_J14333139+3417472'], ["FREQ", "FREQ", "TIME", "TIME"], ['freq_stokesI_20june', 'freq_stokesV_20june', 'time_stokesI_20june', 'time_stokesV_20june '], ["I", "V", "I", "V"], ["9", "9", "9", "9"], [10000, 10000, 10000, 10000]
def forced_fitter():
	forced_fitter_path = '/home/kas/Documents/MRP1/FF_txtfiles/'
	scripts_path = '/home/kas/Documents/MRP1/scripts/'
	
	negative = '--negative '

	hdulist = fits.open(scripts_path+'WABIFAT_input_data_16july.fits') #reads in propermotion corrected coordinates.
	tbdata = hdulist[1].data
	hdulist.close()
	names = tbdata['star']
	ra_lotss = tbdata['RA_gaia_to_lotss_epoch'] #pm-corrected coordinates
	dec_lotss = tbdata['DEC_gaia_to_lotss_epoch'] #pm-corrected coordinates
	lotss_fields = tbdata['field']
	for name, ForT, thing, stokes, lofar_number, n_iter in zip(['GJ3729', 'GJ3729'], ["TIME", "TIME"], ['time_stokesI_15july', 'time_stokesV_15july'], ["I", "V"], ["9", "9"], [10000, 10000]):
		print(name)
		ind = np.where(name == names)
		print(ind)
		ra = ra_lotss[ind]
		print(ra)
		dec = dec_lotss[ind] 
		print(dec)
		field = lotss_fields[ind][0]
		
		wabifatname = name+'_'+thing
		os.system('mkdir '+forced_fitter_path+wabifatname)
		os.system('mkdir '+forced_fitter_path+wabifatname+'/tables')
		os.system('mkdir '+forced_fitter_path+wabifatname+'/results')
		os.system('scp veken@pczaal2.strw.leidenuniv.nl:/net/lofar'+lofar_number+'/data2/veken/wabifat_folder/'+wabifatname+'/goodfits/*image.fits '+forced_fitter_path+wabifatname)
		os.system('scp veken@pczaal2.strw.leidenuniv.nl:/net/lofar'+lofar_number+'/data2/veken/wabifat_folder/'+wabifatname+'/misfits/*image.fits '+forced_fitter_path+wabifatname)
		filenames = glob.glob(forced_fitter_path+wabifatname+'/*-image.fits')
		
		for filename in filenames:
			print(filename)
			os.system('BANE '+filename)
			# read header for psf
			hdu = fits.open(filename)
			bmaj = np.array([hdu[0].header['bmaj']*3600])
			bmin = np.array([hdu[0].header['bmin']*3600]) 
			bpa = np.array([hdu[0].header['bpa']])
			print('bmaj', bmaj)
			print('bmin', bmin)
			print('bpa', bpa)
			# writing input file for Aegean
			cols = fits.ColDefs([
			fits.Column(name='ra', format='D',
			array=ra),
			fits.Column(name='dec', format='D',
			array=dec),
			fits.Column(name='peak_flux', format='E',
			array=np.array([1.])),
			fits.Column(name='a', format='E',
			array=bmaj),
			fits.Column(name='b', format='E',
			array=bmin),
			fits.Column(name='pa', format='E',
			array=bpa),
			fits.Column(name='psf_a', format='E',
			array=bmaj),
			fits.Column(name='psf_b', format='E',
			array=bmin),
			fits.Column(name='psf_pa', format='E',
			array=bpa),
			])
			print('cols', cols)
			hdu_save = fits.BinTableHDU.from_columns(cols)
			hdu_save.writeto(forced_fitter_path+wabifatname+'/'+wabifatname+'_input.fits',overwrite=True)
			print("Here Aegean:")
			os.system('aegean --version')
			os.system('aegean --autoload --priorized 1 '+negative+' --input '+forced_fitter_path+wabifatname+'/'+wabifatname+'_input.fits --floodclip -1 --table '+filename.split('image.fits')[0]+'_output.fits '+filename)
		
		#move the comp.fits tables to a separate folder.
		os.system('mv '+forced_fitter_path+wabifatname+'/*comp* '+forced_fitter_path+wabifatname+'/tables/')
		
		#sorts the filenames in a list.
		list_random = []
		list_sorted = [] 
		for random_name in os.listdir(forced_fitter_path+wabifatname+'/tables/'):
			list_random.append(str(random_name))
		for sorted_name in sorted(list_random):
			list_sorted.append(sorted_name)

		#saves the peakflux, err_peakflux, local_rms from the comp.fits tables.
		peak_flux_list = []
		err_peak_flux_list = []
		local_rms_list = []
		for table in list_sorted:  
			fitstable = fits.open(forced_fitter_path+wabifatname+'/tables/'+table)  
			peak_flux = fitstable[1].data['peak_flux']
			err_peak_flux = fitstable[1].data['err_peak_flux']
			local_rms = fitstable[1].data['local_rms']
			if peak_flux < 0:
				local_rms = (-1)*local_rms
			peak_flux_list.append(peak_flux)
			err_peak_flux_list.append(err_peak_flux)
			local_rms_list.append(local_rms)
		np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_peak_flux.txt', peak_flux_list, delimiter=',')	
		np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_err_peak_flux.txt', err_peak_flux_list, delimiter=',')	
		np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_local_rms.txt', local_rms_list, delimiter=',')	
		
		os.system('rm -r '+forced_fitter_path+wabifatname+'/tables')
		os.system('rm -r '+forced_fitter_path+wabifatname+'/*rms.fits')
		os.system('rm -r '+forced_fitter_path+wabifatname+'/*bkg.fits')
		os.system('rm -r '+forced_fitter_path+wabifatname+'/*image.fits')

		ForT = ForT
		starname = name #Don't use spaces
		starname_forced = wabifatname
		starfield = field
		n_iter = n_iter
		seedclip = 4.0
		autothreshold = 2
		SNR_accept = 5
		stokes = stokes
		path = "/home/kas/Documents/MRP1/results/"+name+"/"+thing+"/results/"
		print(starfield)
		#Frequency Run:
		if ForT == "FREQ":
			#reads in the three files made by the forced_fitter definition.
			fluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_peak_flux.txt", unpack = True)
			rmsarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_local_rms.txt" ,unpack=True)
			errfluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_err_peak_flux.txt",unpack=True)
			print(len(fluxarray))
			print(len(rmsarray))
			print(len(errfluxarray))
	
			freqarray = np.loadtxt(path+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
			final_slice_radius = np.loadtxt(path+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
			misfit_centralfreqarray = np.loadtxt(path+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)
			final_misfit_radius = np.loadtxt(path+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokes+".txt",unpack=True)	

			

			det_freq= freqarray.tolist()
			if isinstance(det_freq, float):
				det_freq = [det_freq]
			nondet_freq = misfit_centralfreqarray.tolist()
			if isinstance(nondet_freq, float):
				nondet_freq = [nondet_freq]
			det_slice = final_slice_radius.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radius.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

		
			totalfreq_not_sorted = det_freq+nondet_freq
			totalslice_not_sorted = det_slice+nondet_slice
			print(len(totalfreq_not_sorted))
			print(len(totalslice_not_sorted))

			for z in range(len(totalfreq_not_sorted)): 
				swap = z + np.argmin(totalfreq_not_sorted[z:])
				(totalfreq_not_sorted[z], totalfreq_not_sorted[swap]) = (totalfreq_not_sorted[swap], totalfreq_not_sorted[z])
				(totalslice_not_sorted[z], totalslice_not_sorted[swap]) = (totalslice_not_sorted[swap], totalslice_not_sorted[z])

			totalfreq_sorted = np.asarray(totalfreq_not_sorted)
			totalslice_sorted = np.asarray(totalslice_not_sorted)

			np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_freq.txt', totalfreq_sorted, delimiter=',')	
			np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_slice.txt', totalslice_sorted, delimiter=',')

			lowest_freq = totalfreq_sorted[0]-totalslice_sorted[0]
			highest_freq = totalfreq_sorted[-1]+totalslice_sorted[-1]
			freqx0 = [lowest_freq, highest_freq] 
			freqy0=[0,0]			

			print(len(totalfreq_sorted))
			print(len(totalslice_sorted))
			#plot
			plt.figure()
			plt.plot(freqx0,freqy0,c='green')
			plt.errorbar(totalfreq_sorted, fluxarray*1000.0, xerr = totalslice_sorted, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+' detections', fmt='.') 
			plt.errorbar(totalfreq_sorted, rmsarray*1000.0, xerr = totalslice_sorted, c = 'aqua', label = 'detection, local RMS',fmt='.') 
			#plt.errorbar(misfit_centralfreqarray, misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
			plt.title('Forced Fit in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
			plt.legend(loc='best')
			plt.xlabel('Frequency (MHz)')
			plt.ylabel('Fluxdensity (mJy/Beam)')
			plt.savefig(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_FF.eps',format='eps')
			#plt.show()		
			
		#Time run:
		if ForT == "TIME": 
			#reads in the three files made by the forced_fitter definition.
			fluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_peak_flux.txt", unpack = True)
			rmsarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_local_rms.txt" ,unpack=True)
			errfluxarray = np.loadtxt(forced_fitter_path+starname_forced+"/results/"+starname_forced+"_forced_err_peak_flux.txt",unpack=True)
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
			if isinstance(det_time, float):
				det_time = [det_time]
			nondet_time = misfit_centraltimearray.tolist()
			if isinstance(nondet_time, float):
				nondet_time = [nondet_time]
			det_slice = final_slice_radius.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radius.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

			totaltime_not_sorted = det_time+nondet_time
			totalslice_not_sorted = det_slice+nondet_slice

			for z in range(len(totaltime_not_sorted)): 
				swap = z + np.argmin(totaltime_not_sorted[z:])
				(totaltime_not_sorted[z], totaltime_not_sorted[swap]) = (totaltime_not_sorted[swap], totaltime_not_sorted[z])
				(totalslice_not_sorted[z], totalslice_not_sorted[swap]) = (totalslice_not_sorted[swap], totalslice_not_sorted[z])

			totaltime_sorted = np.asarray(totaltime_not_sorted)
			totalslice_sorted = np.asarray(totalslice_not_sorted)
			
			np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_time.txt', totaltime_sorted, delimiter=',')	
			np.savetxt(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_forced_slice.txt', totalslice_sorted, delimiter=',')

			#determining the lowest time value and highest time value such that we can plot the green y=0 line I always plot. 
			lowest_time = totaltime_sorted[0]-totalslice_sorted[0]
			highest_time = totaltime_sorted[-1]+totalslice_sorted[-1]
			timex0 = [lowest_time, highest_time] 
			timey0=[0,0]
			
			#plot
			plt.figure()
			plt.plot(timex0,timey0,c='green')
			plt.errorbar(totaltime_sorted-lowest_time, fluxarray*1000.0, xerr = totalslice_sorted, yerr = errfluxarray*1000, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(fluxarray))+' detections', fmt='.') 
			plt.errorbar(totaltime_sorted-lowest_time, rmsarray*1000.0, xerr = totalslice_sorted, c = 'aqua', label = 'detection, local RMS',fmt='.') 
			#plt.errorbar(misfit_centraltimearray-lowest_time, misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection, '+str(SNR_accept)+'*RMS',fmt='.') 
			plt.title('Forced Fit in Stokes '+stokes+' for '+starname+', field: '+starfield+' (n='+str(n_iter)+', AT='+str(autothreshold)+', SC='+str(seedclip)+')')
			plt.legend(loc='best')
			plt.xlabel('Time (hrs)')
			plt.ylabel('Fluxdensity (mJy/Beam)')
			plt.savefig(forced_fitter_path+wabifatname+'/results/'+wabifatname+'_FF.eps',format='eps')
			#plt.show()		
			
forced_fitter()











#['CR_Dra0', 'CR_Dra0', 'CR_Dra0', 'CR_Dra0', 'CR_Dra4', 'CR_Dra4', 'CR_Dra4', 'CR_Dra4', 'CR_Dra11', 'CR_Dra11', 'CR_Dra11', 'CR_Dra11'], ["FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_9april', 'freq_stokesV_14april', 'time_stokesI_4june', 'time_stokesV_16june' ], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "9", "9", "9", "9", "9", "9", "9", "9", "3", "4" , "9", "9" ], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]

#Forced fitting mechanism for WABIFAT results.


#['CR_Dra0', 'CR_Dra0', 'CR_Dra4', 'CR_Dra4', 'CR_Dra11', 'CR_Dra11'], ["TIME", "TIME", "TIME", "TIME", "TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_7june', 'time_stokesI_31may', 'time_stokesV_7june', 'time_stokesI_4june', 'time_stokesV_16june' ], ["I", "V", "I", "V", "I", "V"], [ "9", "9", "9", "9", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000]


#['CR_Dra1', 'CR_Dra1', 'CR_Dra2', 'CR_Dra2', 'CR_Dra2', 'CR_Dra2', 'CR_Dra3', 'CR_Dra3', 'CR_Dra3', 'CR_Dra3', 'CR_Dra5', 'CR_Dra5', 'CR_Dra5', 'CR_Dra5'], ["TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_10april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april' , 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june'], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "9", "9", "2", "4", "9", "9", "2", "4", "9", "9", "2", "4", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000, 10000, 10000, 10000]

#['CR_Dra6', 'CR_Dra6', 'CR_Dra6', 'CR_Dra6', 'CR_Dra7', 'CR_Dra7', 'CR_Dra7', 'CR_Dra7', 'CR_Dra8', 'CR_Dra8', 'CR_Dra8', 'CR_Dra8'], ["FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME", "FREQ", "FREQ", "TIME", "TIME"], ['freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_10april', 'freq_stokesV_14april' , 'time_stokesI_31may', 'time_stokesV_7june', 'freq_stokesI_8april', 'freq_stokesV_14april', 'time_stokesI_31may', 'time_stokesV_7june'], ["I", "V", "I", "V", "I", "V", "I", "V", "I", "V", "I", "V"], [ "2", "4", "9", "9", "2", "4", "9", "9", "2", "4", "9", "9"], [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000, 10000, 10000, 10000]

#['LP212-62', 'LP212-62'], ["TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_31may'], ["I", "V"], [ "9", "9"], [1000, 1000]


#zip(['CR_Dra0', 'CR_Dra1', 'CR_Dra2', 'CR_Dra3', 'CR_Dra4'], ["ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1", "ELAIS-N1"], ["FREQ", "FREQ", "FREQ", "FREQ", "FREQ"], ['freq_stokesI_24june', 'freq_stokesI_8april', 'freq_stokesI_10april', 'freq_stokesI_8april', 'freq_stokesI_8april'], ['freq_stokesV_14april', 'freq_stokesV_14april', 'freq_stokesV_14april', 'freq_stokesV_14april', 'freq_stokesV_14april'] )

#CR_Dra0: 'freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june'
#CR_Dra4: 'freq_stokesI_2july', 'freq_stokesV_2july', 'time_stokesI_31may', 'time_stokesV_7june'
#CR_Dra11: 'freq_stokesI_9april', 'freq_stokesV_14april', 'time_stokesI_4june', 'time_stokesV_16june'

#['LP212-62', 'LP212-62'], ["TIME", "TIME"], ['time_stokesI_31may', 'time_stokesV_31may'], ["I", "V"], [ "9", "9"], [1000, 1000]
#['LP212-62', 'LP212-62', 'CR_Dra0', 'CR_Dra0', 'CR_Dra1', 'CR_Dra1', 'CR_Dra2', 'CR_Dra2', 'CR_Dra3', 'CR_Dra3', 'CR_Dra4', 'CR_Dra4', 'CR_Dra5', 'CR_Dra5', 'CR_Dra6', 'CR_Dra6', 'CR_Dra7', 'CR_Dra7', 'CR_Dra8', 'CR_Dra8', 'CR_Dra9', 'CR_Dra10', 'CR_Dra11', 'CR_Dra11', 'CR_Dra12', 'CR_Dra14', 'CR_Dra15', 'CR_Dra15', 'CR_Dra16', 'CR_Dra18', 'CR_Dra19', 'CR_Dra19', 'CR_Dra20', 'CR_Dra20', '2MASS_J14333139+3417472', '2MASS_J14333139+3417472', 'FK_Com2', 'FK_Com2', 'AD_Leo', 'AD_Leo', 'BD+42_2437', 'BD+42_2437', 'CW_Uma2', 'CW_Uma2', 'GJ450', 'GJ450', 'GJ1151', 'GJ1151', 'OU_And2', 'OU_And2', 'HAT_182-00605', '2MASS_J09481615+511451', '2MASS_J09481615+511451', 'G_240-45', 'GJ412', 'YY_Gem', 'GJ3729'], ["FREQ", "TIME", "FREQ", "TIME","FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "FREQ", "TIME", "FREQ", "FREQ", "FREQ", "FREQ"], ['freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'time_stokesI_31may', 'time_stokesI_31may', 'freq_stokesI_9april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_10april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_20june', 'time_stokesI_20june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_18june', 'time_stokesI_20june', 'freq_stokesI_22june', 'time_stokesI_22june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_22june', 'time_stokesI_22june', 'freq_stokesI_24june', 'time_stokesI_24june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_15july', 'freq_stokesI_15july', 'time_stokesI_15july', 'freq_stokesI_15july', 'freq_stokesI_15july', 'freq_stokesI_17july', 'freq_stokesI_15july'], ['freq_stokesV_3july', 'time_stokesV_31may', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'time_stokesV_7june', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_20june', 'time_stokesV_20june', 'freq_stokesV_19june', 'time_stokesV_5july', 'freq_stokesV_19june', 'time_stokesV_10july', 'freq_stokesV_22june', 'time_stokesV_6july', 'freq_stokesV_19june', 'time_stokesV_9july', 'freq_stokesV_22june', 'time_stokesV_6july', 'freq_stokesV_24june', 'time_stokesV_24june', 'freq_stokesV_19june', 'time_stokesV_5july', 'freq_stokesV_15july', 'freq_stokesV_15july', 'time_stokesV_15july', 'freq_stokesV_15july', 'freq_stokesV_15july', 'freq_stokesV_17july', 'freq_stokesV_15july'], [100, 1000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]
#Polarization Fraction Calculator. 
def polfrac_forced_plot():
	forced_fitter_path = '/home/kas/Documents/MRP1/FF_txtfiles/'
	scripts_path = '/home/kas/Documents/MRP1/scripts/'
	polfractions_path = '/home/kas/Documents/MRP1/polfrac_txtfiles2/'

	hdulist = fits.open(scripts_path+'WABIFAT_input_data_16july.fits') #reads in propermotion corrected coordinates.
	tbdata = hdulist[1].data
	hdulist.close()
	names = tbdata['star']
	channelwidth_freq = tbdata['MHz_per_chan']
	channelwidth_time = tbdata['sec_per_chan']
	lotss_fields = tbdata['field']
	
	for name, ForT, thingI, thingV, n_iter in zip(['GJ3729'], ["TIME"], ['time_stokesI_15july'], ['time_stokesV_15july'], [10000]):
		print(name, ForT)
		ind = np.where(name == names)
		print(ind)
		channelwidth_freq_now = 0.9*channelwidth_freq[ind]
		channelwidth_time_now = (0.9*channelwidth_time[ind])/3600
		field = lotss_fields[ind][0]

		wabifatnameI = name+'_'+thingI
		wabifatnameV = name+'_'+thingV		

		ForT = ForT
		starname = name #Don't use spaces
		starfield = field 
		
		starname_forcedI = wabifatnameI
		starname_forcedV = wabifatnameV
		
		stokesI = "I"
		stokesV = "V"
		
		n_iter = n_iter
		seedclip = 4.0
		autothreshold = 2
		SNR_accept = 5
		
		pathI = "/home/kas/Documents/MRP1/results/"+name+"/"+thingI+"/results/"
		pathV = "/home/kas/Documents/MRP1/results/"+name+"/"+thingV+"/results/"

		os.system('mkdir '+polfractions_path+name+'_'+ForT) 

		if ForT == "FREQ":
			print("if ForT = FREQ is satisfied")
			#STOKES I:
			#loads in the three files made by the forced_fitter definition.
			fluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_peak_flux.txt", unpack = True)
			rmsarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_local_rms.txt" ,unpack=True)
			errfluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_err_peak_flux.txt",unpack=True)
			#print(len(fluxarrayI))
			#print(len(rmsarrayI))
			#print(len(errfluxarrayI))
			#loads in the stuff from WABIFAT results which are not changed by the force fitter
			freqarrayI = np.loadtxt(pathI+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			final_slice_radiusI = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			misfit_centralfreqarrayI = np.loadtxt(pathI+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			final_misfit_radiusI = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	
			#print(len(freqarrayI))
			#print(len(final_slice_radiusI))
			#print(len(misfit_centralfreqarrayI))
			#print(len(final_misfit_radiusI))

			########################################################## Combines det+non_det ###########################################################
			det_freq = freqarrayI.tolist()
			if isinstance(det_freq, float):
				det_freq = [det_freq]
			nondet_freq = misfit_centralfreqarrayI.tolist()
			if isinstance(nondet_freq, float):
				nondet_freq = [nondet_freq]
			det_slice = final_slice_radiusI.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radiusI.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

		
			totalfreq_not_sortedI = det_freq+nondet_freq
			totalslice_not_sortedI = det_slice+nondet_slice

			for z in range(len(totalfreq_not_sortedI)): 
				swap = z + np.argmin(totalfreq_not_sortedI[z:])
				(totalfreq_not_sortedI[z], totalfreq_not_sortedI[swap]) = (totalfreq_not_sortedI[swap], totalfreq_not_sortedI[z])
				(totalslice_not_sortedI[z], totalslice_not_sortedI[swap]) = (totalslice_not_sortedI[swap], totalslice_not_sortedI[z])

			totalfreq_sortedI = np.asarray(totalfreq_not_sortedI)
			totalslice_sortedI = np.asarray(totalslice_not_sortedI)

			print(len(totalfreq_sortedI))
			print(len(totalslice_sortedI))
			###########################################################################################################################################
			

			#STOKES V:
			#loads in the three files made by the forced_fitter definition.
			fluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_peak_flux.txt", unpack = True)
			rmsarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_local_rms.txt" ,unpack=True)
			errfluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_err_peak_flux.txt",unpack=True)
		#	print(len(fluxarrayV))
			#print(len(rmsarrayV))
			#print(len(errfluxarrayV))
			
			#loads in the stuff from WABIFAT results which are not changed by the force fitter
			freqarrayV = np.loadtxt(pathV+ForT+"_freq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			misfit_centralfreqarrayV = np.loadtxt(pathV+ForT+"_final_misfits_centralfreq_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	
			#print(len(freqarrayV))
			#print(len(final_slice_radiusV))
			#print(len(misfit_centralfreqarrayV))
			#print(len(final_misfit_radiusV))

			########################################################## Combines det+non_det ###########################################################
			det_freq = freqarrayV.tolist()
			if isinstance(det_freq, float):
				det_freq = [det_freq]
			nondet_freq = misfit_centralfreqarrayV.tolist()
			if isinstance(nondet_freq, float):
				nondet_freq = [nondet_freq]
			det_slice = final_slice_radiusV.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radiusV.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

		
			totalfreq_not_sortedV = det_freq+nondet_freq
			totalslice_not_sortedV = det_slice+nondet_slice

			for z in range(len(totalfreq_not_sortedV)): 
				swap = z + np.argmin(totalfreq_not_sortedV[z:])
				(totalfreq_not_sortedV[z], totalfreq_not_sortedV[swap]) = (totalfreq_not_sortedV[swap], totalfreq_not_sortedV[z])
				(totalslice_not_sortedV[z], totalslice_not_sortedV[swap]) = (totalslice_not_sortedV[swap], totalslice_not_sortedV[z])

			totalfreq_sortedV = np.asarray(totalfreq_not_sortedV)
			totalslice_sortedV = np.asarray(totalslice_not_sortedV)
			print(len(totalfreq_sortedV))
			print(len(totalslice_sortedV))
			###########################################################################################################################################


			jump_position_list = []
			leftI0 = totalfreq_sortedI[0]-totalslice_sortedI[0]
			leftV0 = totalfreq_sortedV[0]-totalslice_sortedV[0]
			if  leftV0 >= leftI0 :
				jump_position_list.append(leftI0)
			else: 
				jump_position_list.append(leftV0)
			for i in range(len(totalfreq_sortedI)):
				rightI = totalfreq_sortedI[i]+totalslice_sortedI[i]
				jump_position_list.append(rightI)
			for j in range(len(totalfreq_sortedV)):
				rightV = totalfreq_sortedV[j]+totalslice_sortedV[j]
				jump_position_list.append(rightV)
			print("jump_position_list:", jump_position_list)
			print("jump_position_list, length:", len(jump_position_list))
				
			jump_position_array = np.asarray(jump_position_list)
			for z in range(len(jump_position_array)): 
				swap = z + np.argmin(jump_position_array[z:])
				(jump_position_array[z], jump_position_array[swap]) = (jump_position_array[swap], jump_position_array[z])
			sorted_jump_list = jump_position_array.tolist()
			print("sorted_jump_list:", sorted_jump_list)
			print("sorted_jump_list, length:", len(sorted_jump_list))
				
			#final_jump_list = [sorted_jump_list[0]]
			#for k in range(len(sorted_jump_list)-1):
			#	diff = sorted_jump_list[k+1] - sorted_jump_list[k]
			#	if sorted_jump_list[k] != sorted_jump_list[k+1] and diff > channelwidth_freq_now:
			#		final_jump_list.append(sorted_jump_list[k+1])
			#final_jump_list = sorted_jump_list
			final_jump_list = [sorted_jump_list[0]]
			for k in range(len(sorted_jump_list)-1):
				diff = sorted_jump_list[k+1] - sorted_jump_list[k]
				if diff < channelwidth_freq_now:
					a, b = sorted_jump_list.index(sorted_jump_list[k+1]), sorted_jump_list.index(sorted_jump_list[k])
					sorted_jump_list[b], sorted_jump_list[a] = sorted_jump_list[a], sorted_jump_list[b]
			final_jump_list = sorted_jump_list
			print("final_jump_list:", final_jump_list)
			print("final_jump_list, length", len(final_jump_list))
			
			frac_fluxlist = []
			frac_freqlist = []
			frac_radiuslist = []
			frac_errlist = []
			jump_diff_list = []
			SNR_list = []
			for l in range(len(final_jump_list)-1):
				jump_diff = final_jump_list[l+1] - final_jump_list[l]
				if jump_diff < channelwidth_freq_now:
					print('skip l for loop because jump_diff =', jump_diff)
					continue
				jump_diff_list.append(jump_diff)
				print("final_jump_list[l+1]:", final_jump_list[l+1])
				print("final_jump_list[l]:", final_jump_list[l])
				print("jump_diff:", jump_diff)
				#if sorted_jump_list[l+2] != sorted_jump_list[-1] and sorted_jump_list[l+1] != sorted_jump_list[-1]:
				#	next_diff = sorted_jump_list[l+2] - sorted_jump_list[l]
				#else: 
				#	next_diff = 0
				
				for m in range(len(totalfreq_sortedI)):
					for n in range(len(totalfreq_sortedV)):
						if rmsarrayI[m] == 0 or rmsarrayV[n] == 0:
							continue
						SNR_I = fluxarrayI[m]/rmsarrayI[m]
						SNR_V =	fluxarrayV[n]/rmsarrayV[n]
						#print('SNR_I = ',SNR_I, 'SNR_V = ', SNR_V)
						if SNR_I <= 3 or SNR_V <= 3:
							SNR_list.append([totalfreq_sortedI[m], totalfreq_sortedV[n], SNR_I, SNR_V])
							continue
						new_radius = (final_jump_list[l+1]-final_jump_list[l])/2.0
						new_freq = final_jump_list[l]+new_radius
						polfrac = abs(fluxarrayV[n]/fluxarrayI[m])

						#if next_diff >= channelwidth_freq_now: 
						#	leftI = totalfreq_sortedI[m]-totalslice_sortedI[m]
						#	rightI = totalfreq_sortedI[m]-totalslice_sortedI[m]
						#	leftV = totalfreq_sortedV[n]-totalslice_sortedV[n]
						#	rightV = totalfreq_sortedV[n]-totalslice_sortedV[n]	
						#else: 
						leftI = totalfreq_sortedI[m]-totalslice_sortedI[m]
						rightI = totalfreq_sortedI[m]+totalslice_sortedI[m]
						leftV = totalfreq_sortedV[n]-totalslice_sortedV[n]
						rightV = totalfreq_sortedV[n]+totalslice_sortedV[n]	
						
				
						diff_leftI = abs(leftI-final_jump_list[l])
						diff_rightI = abs(rightI-final_jump_list[l+1])
						diff_leftV = abs(leftV-final_jump_list[l])
						diff_rightV = abs(rightV-final_jump_list[l+1])
						

						#intervalI and intervalV are the same width.
						if diff_leftI < channelwidth_freq_now and diff_rightI < channelwidth_freq_now and diff_leftV < channelwidth_freq_now and diff_rightV < channelwidth_freq_now: 
							print("1: if diff_leftI < channelwidth_freq_now and diff_rightI < channelwidth_freq_now and diff_leftV < channelwidth_freq_now and diff_rightV < channelwidth_freq_now:")
							polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_freqlist.append(new_freq)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)
								print("1: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						#IntervalI corresponds with jumps, intervalV will be larger.
						if diff_leftI < channelwidth_freq_now and diff_rightI < channelwidth_freq_now:
							print("2: elif diff_leftI < channelwidth_freq_now and diff_rightI < channelwidth_freq_now:")
							#intervalV is longer than intervalI to the right OR left. 
							if (diff_leftV < channelwidth_freq_now and diff_rightV >= channelwidth_freq_now) or (diff_leftV >= channelwidth_freq_now and diff_rightV < channelwidth_freq_now): 
								print("2.1")
								fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fitfraction**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_freqlist.append(new_freq)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("2.1: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
							#intervalV is bigger than intervalI to the left AND right.
							if diff_leftV >= channelwidth_freq_now and leftV < final_jump_list[l] and diff_rightV >= channelwidth_freq_now and rightV > final_jump_list[l+1]: 
								print("2.2")
								fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fitfraction**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_freqlist.append(new_freq)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("2.2: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
						#IntervalV corresponds with jumps, intervalI will be larger. 
						if diff_leftV < channelwidth_freq_now and diff_rightV < channelwidth_freq_now:
							print("3: elif diff_leftV < channelwidth_freq_now and diff_rightV < channelwidth_freq_now:")
							#intervalI is longer than intervalV to the right OR left. 
							if (diff_leftI < channelwidth_freq_now and diff_rightI >= channelwidth_freq_now) or (diff_leftI >= channelwidth_freq_now and diff_rightI < channelwidth_freq_now): 
								print("3.1")
								#fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								fitfraction = (totalslice_sortedV[n]/totalslice_sortedI[m])**(-1.0)
								polfrac_err = polfrac*(   ((errfluxarrayI[m]*fitfraction**(0.5))/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_freqlist.append(new_freq)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("3.1: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
							#intervalV is bigger than intervalI to the left AND right.
							if diff_leftI >= channelwidth_freq_now and leftI < final_jump_list[l] and diff_rightI >= channelwidth_freq_now and rightI > final_jump_list[l+1]:
								print("3.2") 
								#fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								fitfraction = (totalslice_sortedV[n]/totalslice_sortedI[m])**(-1.0)
								polfrac_err = polfrac*(   ((errfluxarrayI[m]*fitfraction**(0.5))/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_freqlist.append(new_freq)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("3.2: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
						#Leftboundary of intervalI corresponds with jump, rightboundary of intervalV corresponds with jump.
						if diff_leftI < channelwidth_freq_now and diff_rightI >= channelwidth_freq_now and diff_leftV >= channelwidth_freq_now and diff_rightV < channelwidth_freq_now:
							print("4: elif diff_leftI < channelwidth_freq_now and diff_rightI >= channelwidth_freq_now and diff_leftV >= channelwidth_freq_now and diff_rightV < channelwidth_freq_now:")
							fraction_ofI = (2*totalslice_sortedI[m])/jump_diff
							fraction_ofV = (2*totalslice_sortedV[n])/jump_diff
							polfrac_err = polfrac*(   ((errfluxarrayI[m]*fraction_ofI**(0.5))/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fraction_ofV**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_freqlist.append(new_freq)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)	
								print("4: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						#Leftboundary of intervalV corresponds with jump, rightboundary of intervalI corresponds with jump.
						if diff_leftV < channelwidth_freq_now and diff_rightV >= channelwidth_freq_now and diff_leftI >= channelwidth_freq_now and diff_rightI < channelwidth_freq_now:
							print("5: elif diff_leftV < channelwidth_freq_now and diff_rightV >= channelwidth_freq_now and diff_leftI >= channelwidth_freq_now and diff_rightI < channelwidth_freq_now:")
							fraction_ofI = (2*totalslice_sortedI[m])/jump_diff
							fraction_ofV = (2*totalslice_sortedV[n])/jump_diff
							polfrac_err = polfrac*(   ((errfluxarrayI[m]*fraction_ofI**(0.5))/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fraction_ofV**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_freqlist.append(new_freq)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)	
								print("5: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_flux.txt', frac_fluxlist, delimiter=',')	
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_freq.txt', frac_freqlist, delimiter=',')	
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_radius.txt', frac_radiuslist, delimiter=',')
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_err.txt', frac_errlist, delimiter=',')				
				
			print('final_jump_list', final_jump_list)
			print('jump_diff_list:',jump_diff_list)
			frac_fluxarray = np.asarray(frac_fluxlist)
			frac_freqarray = np.asarray(frac_freqlist)
			frac_radiusarray = np.asarray(frac_radiuslist)
			frac_errarray = np.asarray(frac_errlist)
			print("fraction:", frac_fluxarray)
			print("freq:", frac_freqarray)
			print("radius:", frac_radiusarray)
			print("error:", frac_errarray)
			print("SNR_list (freqI, freqV, SNR_I, SNR_V):", SNR_list)

			freqx0 = [final_jump_list[0], final_jump_list[-1]]
			freqy0 = [0,0]

			
			#, uplims=True, lolims=True
			plt.figure()
			plt.plot(freqx0,freqy0,c='green')
			plt.errorbar(frac_freqarray,frac_fluxarray, xerr = frac_radiusarray, yerr=frac_errarray, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(frac_fluxarray))+"x stokes(V/I)", fmt='.') 
			plt.title("Polarization Fraction (V/I) for star: "+starname+", in the frequency domain")
			plt.legend(loc='best')
			plt.xlabel('Frequency (MHz)')
			plt.ylabel('Fraction V/I')
			plt.savefig(polfractions_path+name+'_'+ForT+'/'+name+'_polfrac.eps',format='eps')
			plt.close()
			#plt.show()
			#plt.errorbar(misfit_centralfreqarray,misfit_rmsarray*1000.0*SNR_accept, xerr = final_misfit_radius, yerr=misfit_rmsarray*800, uplims=True, solid_capstyle='projecting', capsize=4, c = 'red', label = 'non-detection stokes I, '+str(SNR_accept)+'*RMS',fmt='.') 
			#plt.scatter(misfits_freq, misfits_zero, marker = 'x',c = 'red', label = str(len(misfits_freq))+' misfit channels')			
		






		
		if ForT == "TIME":
			print("if ForT = TIME is satisfied")
			#STOKES I:
			#loads in the three files made by the forced_fitter definition.
			fluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_peak_flux.txt", unpack = True)
			rmsarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_local_rms.txt" ,unpack=True)
			errfluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_err_peak_flux.txt",unpack=True)
			#print(len(fluxarrayI))
			#print(len(rmsarrayI))
			#print(len(errfluxarrayI))
			#loads in the stuff from WABIFAT results which are not changed by the force fitter
			timearrayI = np.loadtxt(pathI+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			final_slice_radiusI = np.loadtxt(pathI+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			misfit_centraltimearrayI = np.loadtxt(pathI+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)
			final_misfit_radiusI = np.loadtxt(pathI+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesI+".txt",unpack=True)	
			#print(len(timearrayI))
			#print(len(final_slice_radiusI))
			#print(len(misfit_centraltimearrayI))
			#print(len(final_misfit_radiusI))

			########################################################## Combines det+non_det ###########################################################
			det_time = timearrayI.tolist()
			if isinstance(det_time, float):
				det_time = [det_time]
			nondet_time = misfit_centraltimearrayI.tolist()
			if isinstance(nondet_time, float):
				nondet_time = [nondet_time]
			det_slice = final_slice_radiusI.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radiusI.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

		
			totaltime_not_sortedI = det_time+nondet_time
			totalslice_not_sortedI = det_slice+nondet_slice

			for z in range(len(totaltime_not_sortedI)): 
				swap = z + np.argmin(totaltime_not_sortedI[z:])
				(totaltime_not_sortedI[z], totaltime_not_sortedI[swap]) = (totaltime_not_sortedI[swap], totaltime_not_sortedI[z])
				(totalslice_not_sortedI[z], totalslice_not_sortedI[swap]) = (totalslice_not_sortedI[swap], totalslice_not_sortedI[z])

			totaltime_sortedI = np.asarray(totaltime_not_sortedI)
			totalslice_sortedI = np.asarray(totalslice_not_sortedI)

			print(len(totaltime_sortedI))
			print(len(totalslice_sortedI))
			###########################################################################################################################################
			

			#STOKES V:
			#loads in the three files made by the forced_fitter definition.
			fluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_peak_flux.txt", unpack = True)
			rmsarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_local_rms.txt" ,unpack=True)
			errfluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_err_peak_flux.txt",unpack=True)
		#	print(len(fluxarrayV))
			#print(len(rmsarrayV))
			#print(len(errfluxarrayV))
			
			#loads in the stuff from WABIFAT results which are not changed by the force fitter
			timearrayV = np.loadtxt(pathV+ForT+"_time_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			final_slice_radiusV = np.loadtxt(pathV+ForT+"_slice_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			misfit_centraltimearrayV = np.loadtxt(pathV+ForT+"_final_misfits_centraltime_list_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)
			final_misfit_radiusV = np.loadtxt(pathV+ForT+"_final_misfit_radius_"+starname+"_field"+starfield+"_n"+str(n_iter)+"AT"+str(autothreshold)+"SC"+str(int(round(seedclip)))+"pol"+stokesV+".txt",unpack=True)	
			#print(len(timearrayV))
			#print(len(final_slice_radiusV))
			#print(len(misfit_centraltimearrayV))
			#print(len(final_misfit_radiusV))

			########################################################## Combines det+non_det ###########################################################
			det_time = timearrayV.tolist()
			if isinstance(det_time, float):
				det_time = [det_time]
			nondet_time = misfit_centraltimearrayV.tolist()
			if isinstance(nondet_time, float):
				nondet_time = [nondet_time]
			det_slice = final_slice_radiusV.tolist()
			if isinstance(det_slice, float):
				det_slice = [det_slice]
			nondet_slice = final_misfit_radiusV.tolist()
			if isinstance(nondet_slice, float):
				nondet_slice = [nondet_slice]

		
			totaltime_not_sortedV = det_time+nondet_time
			totalslice_not_sortedV = det_slice+nondet_slice

			for z in range(len(totaltime_not_sortedV)): 
				swap = z + np.argmin(totaltime_not_sortedV[z:])
				(totaltime_not_sortedV[z], totaltime_not_sortedV[swap]) = (totaltime_not_sortedV[swap], totaltime_not_sortedV[z])
				(totalslice_not_sortedV[z], totalslice_not_sortedV[swap]) = (totalslice_not_sortedV[swap], totalslice_not_sortedV[z])

			totaltime_sortedV = np.asarray(totaltime_not_sortedV)
			totalslice_sortedV = np.asarray(totalslice_not_sortedV)
			print(len(totaltime_sortedV))
			print(len(totalslice_sortedV))
			###########################################################################################################################################


			jump_position_list = []
			leftI0 = totaltime_sortedI[0]-totalslice_sortedI[0]
			leftV0 = totaltime_sortedV[0]-totalslice_sortedV[0]
			if  leftV0 >= leftI0 :
				jump_position_list.append(leftI0)
			else: 
				jump_position_list.append(leftV0)
			for i in range(len(totaltime_sortedI)):
				rightI = totaltime_sortedI[i]+totalslice_sortedI[i]
				jump_position_list.append(rightI)
			for j in range(len(totaltime_sortedV)):
				rightV = totaltime_sortedV[j]+totalslice_sortedV[j]
				jump_position_list.append(rightV)
			print("jump_position_list:", jump_position_list)
			print("jump_position_list, length:", len(jump_position_list))
				
			jump_position_array = np.asarray(jump_position_list)
			for z in range(len(jump_position_array)): 
				swap = z + np.argmin(jump_position_array[z:])
				(jump_position_array[z], jump_position_array[swap]) = (jump_position_array[swap], jump_position_array[z])
			sorted_jump_list = jump_position_array.tolist()
			print("sorted_jump_list:", sorted_jump_list)
			print("sorted_jump_list, length:", len(sorted_jump_list))
				
			#final_jump_list = [sorted_jump_list[0]]
			#for k in range(len(sorted_jump_list)-1):
			#	diff = sorted_jump_list[k+1] - sorted_jump_list[k]
			#	if sorted_jump_list[k] != sorted_jump_list[k+1] and diff > channelwidth_time_now:
			#		final_jump_list.append(sorted_jump_list[k+1])
			#final_jump_list = sorted_jump_list
			final_jump_list = [sorted_jump_list[0]]
			for k in range(len(sorted_jump_list)-1):
				diff = sorted_jump_list[k+1] - sorted_jump_list[k]
				if diff < channelwidth_time_now:
					a, b = sorted_jump_list.index(sorted_jump_list[k+1]), sorted_jump_list.index(sorted_jump_list[k])
					sorted_jump_list[b], sorted_jump_list[a] = sorted_jump_list[a], sorted_jump_list[b]
			final_jump_list = sorted_jump_list
			print("final_jump_list:", final_jump_list)
			print("final_jump_list, length", len(final_jump_list))
			
			frac_fluxlist = []
			frac_timelist = []
			frac_radiuslist = []
			frac_errlist = []
			jump_diff_list = []
			SNR_list = []
			for l in range(len(final_jump_list)-1):
				jump_diff = final_jump_list[l+1] - final_jump_list[l]
				if jump_diff < channelwidth_time_now:
					print('skip l for loop because jump_diff =', jump_diff)
					continue
				jump_diff_list.append(jump_diff)
				print("final_jump_list[l+1]:", final_jump_list[l+1])
				print("final_jump_list[l]:", final_jump_list[l])
				print("jump_diff:", jump_diff)
				#if sorted_jump_list[l+2] != sorted_jump_list[-1] and sorted_jump_list[l+1] != sorted_jump_list[-1]:
				#	next_diff = sorted_jump_list[l+2] - sorted_jump_list[l]
				#else: 
				#	next_diff = 0
				for m in range(len(totaltime_sortedI)):
					for n in range(len(totaltime_sortedV)):
						if rmsarrayI[m] == 0 or rmsarrayV[n] == 0:
							continue
						SNR_I = fluxarrayI[m]/rmsarrayI[m]
						SNR_V =	fluxarrayV[n]/rmsarrayV[n]
						#print('SNR_I = ',SNR_I, 'SNR_V = ', SNR_V)
						if SNR_I <= 3 or SNR_V <= 3:
							SNR_list.append([totaltime_sortedI[m], totaltime_sortedV[n], SNR_I, SNR_V])
							continue
						new_radius = (final_jump_list[l+1]-final_jump_list[l])/2.0
						new_time = final_jump_list[l]+new_radius
						polfrac = abs(fluxarrayV[n]/fluxarrayI[m])

						#if next_diff >= channelwidth_time_now: 
						#	leftI = totaltime_sortedI[m]-totalslice_sortedI[m]
						#	rightI = totaltime_sortedI[m]-totalslice_sortedI[m]
						#	leftV = totaltime_sortedV[n]-totalslice_sortedV[n]
						#	rightV = totaltime_sortedV[n]-totalslice_sortedV[n]	
						#else: 
						leftI = totaltime_sortedI[m]-totalslice_sortedI[m]
						rightI = totaltime_sortedI[m]+totalslice_sortedI[m]
						leftV = totaltime_sortedV[n]-totalslice_sortedV[n]
						rightV = totaltime_sortedV[n]+totalslice_sortedV[n]	
						
				
						diff_leftI = abs(leftI-final_jump_list[l])
						diff_rightI = abs(rightI-final_jump_list[l+1])
						diff_leftV = abs(leftV-final_jump_list[l])
						diff_rightV = abs(rightV-final_jump_list[l+1])
						

						#intervalI and intervalV are the same width.
						if diff_leftI < channelwidth_time_now and diff_rightI < channelwidth_time_now and diff_leftV < channelwidth_time_now and diff_rightV < channelwidth_time_now: 
							print("1: if diff_leftI < channelwidth_time_now and diff_rightI < channelwidth_time_now and diff_leftV < channelwidth_time_now and diff_rightV < channelwidth_time_now:")
							polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_timelist.append(new_time)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)
								print("1: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						#IntervalI corresponds with jumps, intervalV will be larger.
						if diff_leftI < channelwidth_time_now and diff_rightI < channelwidth_time_now:
							print("2: elif diff_leftI < channelwidth_time_now and diff_rightI < channelwidth_time_now:")
							#intervalV is longer than intervalI to the right OR left. 
							if (diff_leftV < channelwidth_time_now and diff_rightV >= channelwidth_time_now) or (diff_leftV >= channelwidth_time_now and diff_rightV < channelwidth_time_now): 
								print("2.1")
								fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fitfraction**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_timelist.append(new_time)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("2.1: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
							#intervalV is bigger than intervalI to the left AND right.
							if diff_leftV >= channelwidth_time_now and leftV < final_jump_list[l] and diff_rightV >= channelwidth_time_now and rightV > final_jump_list[l+1]: 
								print("2.2")
								fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								polfrac_err = polfrac*(   (errfluxarrayI[m]/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fitfraction**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_timelist.append(new_time)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("2.2: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
						#IntervalV corresponds with jumps, intervalI will be larger. 
						if diff_leftV < channelwidth_time_now and diff_rightV < channelwidth_time_now:
							print("3: elif diff_leftV < channelwidth_time_now and diff_rightV < channelwidth_time_now:")
							#intervalI is longer than intervalV to the right OR left. 
							if (diff_leftI < channelwidth_time_now and diff_rightI >= channelwidth_time_now) or (diff_leftI >= channelwidth_time_now and diff_rightI < channelwidth_time_now): 
								print("3.1")
								#fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								fitfraction = (totalslice_sortedV[n]/totalslice_sortedI[m])**(-1.0)
								polfrac_err = polfrac*(   ((errfluxarrayI[m]*fitfraction**(0.5))/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_timelist.append(new_time)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("3.1: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
							#intervalV is bigger than intervalI to the left AND right.
							if diff_leftI >= channelwidth_time_now and leftI < final_jump_list[l] and diff_rightI >= channelwidth_time_now and rightI > final_jump_list[l+1]:
								print("3.2") 
								#fitfraction = totalslice_sortedV[n]/totalslice_sortedI[m]
								fitfraction = (totalslice_sortedV[n]/totalslice_sortedI[m])**(-1.0)
								polfrac_err = polfrac*(   ((errfluxarrayI[m]*fitfraction**(0.5))/fluxarrayI[m])**2.0  + (errfluxarrayV[n]/fluxarrayV[n])**2.0  )**(0.5)
								if abs(polfrac)-abs(polfrac_err) <= 1: 
									frac_fluxlist.append(polfrac)
									frac_timelist.append(new_time)
									frac_radiuslist.append(new_radius)
									frac_errlist.append(polfrac_err)
									print("3.2: Appended for:", final_jump_list[l], final_jump_list[l+1])
									print("leftI:", leftI)
									print("rightI:", rightI)
									print("leftV:", leftV)
									print("rightV:", rightV)
									print("dlI:", diff_leftI)
									print("drI:", diff_rightI)
									print("dlV:", diff_leftV)
									print("drV:", diff_rightV)
						#Leftboundary of intervalI corresponds with jump, rightboundary of intervalV corresponds with jump.
						if diff_leftI < channelwidth_time_now and diff_rightI >= channelwidth_time_now and diff_leftV >= channelwidth_time_now and diff_rightV < channelwidth_time_now:
							print("4: elif diff_leftI < channelwidth_time_now and diff_rightI >= channelwidth_time_now and diff_leftV >= channelwidth_time_now and diff_rightV < channelwidth_time_now:")
							fraction_ofI = (2*totalslice_sortedI[m])/jump_diff
							fraction_ofV = (2*totalslice_sortedV[n])/jump_diff
							polfrac_err = polfrac*(   ((errfluxarrayI[m]*fraction_ofI**(0.5))/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fraction_ofV**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_timelist.append(new_time)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)	
								print("4: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						#Leftboundary of intervalV corresponds with jump, rightboundary of intervalI corresponds with jump.
						if diff_leftV < channelwidth_time_now and diff_rightV >= channelwidth_time_now and diff_leftI >= channelwidth_time_now and diff_rightI < channelwidth_time_now:
							print("5: elif diff_leftV < channelwidth_time_now and diff_rightV >= channelwidth_time_now and diff_leftI >= channelwidth_time_now and diff_rightI < channelwidth_time_now:")
							fraction_ofI = (2*totalslice_sortedI[m])/jump_diff
							fraction_ofV = (2*totalslice_sortedV[n])/jump_diff
							polfrac_err = polfrac*(   ((errfluxarrayI[m]*fraction_ofI**(0.5))/fluxarrayI[m])**2.0  + ((errfluxarrayV[n]*fraction_ofV**(0.5))/fluxarrayV[n])**2.0  )**(0.5)
							if abs(polfrac)-abs(polfrac_err) <= 1: 
								frac_fluxlist.append(polfrac)
								frac_timelist.append(new_time)
								frac_radiuslist.append(new_radius)
								frac_errlist.append(polfrac_err)	
								print("5: Appended for:", final_jump_list[l], final_jump_list[l+1])
								print("leftI:", leftI)
								print("rightI:", rightI)
								print("leftV:", leftV)
								print("rightV:", rightV)
								print("dlI:", diff_leftI)
								print("drI:", diff_rightI)
								print("dlV:", diff_leftV)
								print("drV:", diff_rightV)
						
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_flux.txt', frac_fluxlist, delimiter=',')	
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_time.txt', frac_timelist, delimiter=',')	
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_radius.txt', frac_radiuslist, delimiter=',')
			np.savetxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_err.txt', frac_errlist, delimiter=',')				
				
			print('final_jump_list', final_jump_list)
			print('jump_diff_list:',jump_diff_list)
			frac_fluxarray = np.asarray(frac_fluxlist)
			frac_timearray = np.asarray(frac_timelist)
			frac_radiusarray = np.asarray(frac_radiuslist)
			frac_errarray = np.asarray(frac_errlist)
			print("fraction:", frac_fluxarray)
			print("time:", frac_timearray)
			print("radius:", frac_radiusarray)
			print("error:", frac_errarray)
			print("SNR_list (timeI, timeV, SNR_I, SNR_V):", SNR_list)
			timex0 = [final_jump_list[0], final_jump_list[-1]]
			timey0 = [0,0]

			#, uplims=True, lolims=True
			plt.figure()
			plt.plot(timex0,timey0,c='green')
			plt.errorbar(frac_timearray-final_jump_list[0],frac_fluxarray, xerr = frac_radiusarray, yerr=frac_errarray, solid_capstyle='projecting', capsize=4, c = 'blue', label = str(len(frac_fluxarray))+"x stokes(V/I)", fmt='.') 
			plt.title("Polarization Fraction (V/I) for star: "+starname+", in the time domain")
			plt.legend(loc='best')
			plt.xlabel('time (hrs)')
			plt.ylabel('Fraction V/I')
			plt.savefig(polfractions_path+name+'_'+ForT+'/'+name+'_polfrac.eps',format='eps')
			plt.close()
			#plt.show()
			
polfrac_forced_plot()





#['OU_And2'], ["TIME"], ['time_stokesI_24june'], ['time_stokesV_5july'], [10000]
#['LP212-62', 'LP212-62', 'CR_Dra0', 'CR_Dra0', 'CR_Dra1', 'CR_Dra1', 'CR_Dra2', 'CR_Dra2', 'CR_Dra3', 'CR_Dra3', 'CR_Dra4', 'CR_Dra4', 'CR_Dra5', 'CR_Dra5', 'CR_Dra6', 'CR_Dra6', 'CR_Dra7', 'CR_Dra7', 'CR_Dra8', 'CR_Dra8', 'CR_Dra9', 'CR_Dra10', 'CR_Dra11', 'CR_Dra11', 'CR_Dra12', 'CR_Dra14', 'CR_Dra15', 'CR_Dra15', 'CR_Dra16', 'CR_Dra18', 'CR_Dra19', 'CR_Dra19', 'CR_Dra20', 'CR_Dra20', '2MASS_J14333139+3417472', '2MASS_J14333139+3417472', 'FK_Com2', 'FK_Com2', 'AD_Leo', 'BD+42_2437', 'CW_Uma2', 'GJ450', 'GJ1151', 'OU_And2', 'OU_And2'], ["FREQ", "TIME", "FREQ", "TIME","FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "FREQ", "FREQ", "FREQ", "FREQ", "FREQ", "TIME"], ['freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'time_stokesI_31may', 'time_stokesI_31may', 'freq_stokesI_9april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_10april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_20june', 'time_stokesI_20june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_18june', 'freq_stokesI_22june', 'freq_stokesI_18june', 'freq_stokesI_22june', 'freq_stokesI_24june', 'freq_stokesI_18june', 'time_stokesI_24june'], ['freq_stokesV_3july', 'time_stokesV_31may', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'time_stokesV_7june', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_20june', 'time_stokesV_20june', 'freq_stokesV_19june', 'time_stokesV_5july', 'freq_stokesV_19june', 'freq_stokesV_22june', 'freq_stokesV_19june', 'freq_stokesV_22june', 'freq_stokesV_24june', 'freq_stokesV_19june', 'time_stokesV_5july'], [100, 1000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]




#['LP212-62', 'LP212-62', 'CR_Dra0', 'CR_Dra0', 'CR_Dra1', 'CR_Dra1', 'CR_Dra2', 'CR_Dra2', 'CR_Dra3', 'CR_Dra3', 'CR_Dra4', 'CR_Dra4', 'CR_Dra5', 'CR_Dra5', 'CR_Dra6', 'CR_Dra6', 'CR_Dra7', 'CR_Dra7', 'CR_Dra8', 'CR_Dra8', 'CR_Dra9', 'CR_Dra10', 'CR_Dra11', 'CR_Dra11', 'CR_Dra12', 'CR_Dra14', 'CR_Dra15', 'CR_Dra15', 'CR_Dra16', 'CR_Dra18', 'CR_Dra19', 'CR_Dra19', 'CR_Dra20', 'CR_Dra20', '2MASS_J14333139+3417472', '2MASS_J14333139+3417472', 'FK_Com2', 'FK_Com2', 'AD_Leo', 'AD_Leo', 'BD+42_2437', 'BD+42_2437', 'CW_Uma2', 'CW_Uma2', 'GJ450', 'GJ450', 'GJ1151', 'GJ1151', 'OU_And2', 'OU_And2', 'HAT_182-00605', 'HAT_182-00605', '2MASS_J09481615+511451', '2MASS_J09481615+511451', 'G_240-45', 'GJ412', 'YY_Gem', 'YY_Gem', 'GJ3729', 'LP259-39'], ["FREQ", "TIME", "FREQ", "TIME","FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "TIME", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "TIME", "FREQ", "FREQ", "FREQ", "TIME", "FREQ", "TIME"], ['freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_2july', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'freq_stokesI_10april', 'time_stokesI_31may', 'freq_stokesI_8april', 'time_stokesI_31may', 'time_stokesI_31may', 'time_stokesI_31may', 'freq_stokesI_9april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_10april', 'time_stokesI_4june', 'time_stokesI_4june', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_9april', 'time_stokesI_4june', 'freq_stokesI_20june', 'time_stokesI_20june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_18june', 'time_stokesI_20june', 'freq_stokesI_22june', 'time_stokesI_22june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_22june', 'time_stokesI_22june', 'freq_stokesI_24june', 'time_stokesI_24june', 'freq_stokesI_18june', 'time_stokesI_24june', 'freq_stokesI_15july', 'time_stokesI_15july', 'freq_stokesI_15july', 'time_stokesI_15july', 'freq_stokesI_15july', 'freq_stokesI_15july', 'freq_stokesI_17july', 'time_stokesI_17july', 'freq_stokesI_15july', 'time_stokesI_16july'], ['freq_stokesV_3july', 'time_stokesV_31may', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_2july', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_7june', 'time_stokesV_7june', 'time_stokesV_7june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'time_stokesV_16june', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_14april', 'time_stokesV_16june', 'freq_stokesV_20june', 'time_stokesV_20june', 'freq_stokesV_19june', 'time_stokesV_5july', 'freq_stokesV_19june', 'time_stokesV_10july', 'freq_stokesV_22june', 'time_stokesV_6july', 'freq_stokesV_19june', 'time_stokesV_9july', 'freq_stokesV_22june', 'time_stokesV_6july', 'freq_stokesV_24june', 'time_stokesV_24june', 'freq_stokesV_19june', 'time_stokesV_5july', 'freq_stokesV_15july', 'time_stokesV_15july', 'freq_stokesV_15july', 'time_stokesV_15july', 'freq_stokesV_15july', 'freq_stokesV_15july', 'freq_stokesV_17july', 'time_stokesV_17july', 'freq_stokesV_15july', 'time_stokesV_16july'], [100, 1000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]


def FF_yerr_bar_fixer(fluxarray, errfluxarray):
	errfluxlist = []
	for i in range(len(fluxarray)):
		if np.isnan(errfluxarray[i]) == True or abs(errfluxarray[i]) == 0:
			errfluxlist.append(0.1*(fluxarray[i]))
		else: 
			errfluxlist.append(errfluxarray[i])
	errfluxarray = np.asarray(errfluxlist)
	return errfluxarray

def frac_yerr_bar_fixer(errarray, fluxarray, x_array, radiusarray):
	errlist = errarray.tolist() 
	fluxlist = fluxarray.tolist()
	x_list = x_array.tolist()
	radiuslist = radiusarray.tolist()
	new_errlist = []
	new_fluxlist = []
	new_x_list = []
	new_radiuslist = []
	print("errlist:", errlist)
	if isinstance(errlist, float):
		print("isinstance ==TRue")
		errlist = [errlist]
		fluxlist = [fluxlist]
		x_list = [x_list]
		radiuslist = [radiuslist] 
	for i in range(len(errlist)):
		if 2*errlist[i] <= fluxlist[i]:
			new_errlist.append(errlist[i])
			new_fluxlist.append(fluxlist[i])
			new_x_list.append(x_list[i])
			new_radiuslist.append(radiuslist[i])
	return np.asarray(new_errlist), np.asarray(new_fluxlist), np.asarray(new_x_list), np.asarray(new_radiuslist)
	
def negative_stokes_I_deleter(fluxarray, x_array, errfluxarray, final_slice_radius, rmsarray):
	x_list = []
	flux_list = []
	final_slice_radius_list = []
	errflux_list = []
	rms_list = []
	for i in range(len(fluxarray)):
		uppervalue = fluxarray[i]+errfluxarray[i]
		if uppervalue >= 0:
			x_list.append(x_array[i])
			flux_list.append(fluxarray[i])
			errflux_list.append(errfluxarray[i])
			final_slice_radius_list.append(final_slice_radius[i])
			rms_list.append(rmsarray[i])
	return np.asarray(x_list), np.asarray(flux_list), np.asarray(errflux_list), np.asarray(final_slice_radius_list), np.asarray(rms_list)


#Plotting mechanism for the combination plot of Stokes I, Stokes V and Stokes V/I.
def FFI_FFV_polfrac_plotter():
	forced_fitter_path = '/home/kas/Documents/MRP1/FF_txtfiles/'
	scripts_path = '/home/kas/Documents/MRP1/scripts/'
	polfractions_path = '/home/kas/Documents/MRP1/polfrac_txtfiles2/'
	triple_path = '/home/kas/Documents/MRP1/triple_plots/'

	hdulist = fits.open(scripts_path+'WABIFAT_input_data_16july.fits') #reads in propermotion corrected coordinates.
	tbdata = hdulist[1].data
	hdulist.close()
	names = tbdata['star']
	channelwidth_freq = tbdata['MHz_per_chan']
	channelwidth_time = tbdata['sec_per_chan']
	lotss_fields = tbdata['field']
	freq_chans = tbdata['nr_freq_chans_for_msfile']
	time_chans = tbdata['nr_time_chans_for_msfile']
	freq_bandwidth_low = tbdata['bandwidth_low']
	freq_bandwidth_up = tbdata['bandwidth_up']
	for name, ForT, thingI, thingV, n_iter in zip(['GJ3729', 'GJ3729'], ["TIME", "TIME"], ['time_stokesI_15july', 'time_stokesI_15july'], ['time_stokesV_15july', 'time_stokesV_15july'], [10000, 10000]):
		print(name, ForT)
		ind = np.where(name == names)
		print(ind)
		field = lotss_fields[ind][0]

		wabifatnameI = name+'_'+thingI
		wabifatnameV = name+'_'+thingV		

		ForT = ForT
		starname = name #Don't use spaces
		starfield = field 

		starname_forcedI = wabifatnameI
		starname_forcedV = wabifatnameV

		stokesI = "I"
		stokesV = "V"

		n_iter = n_iter
		seedclip = 4.0
		autothreshold = 2
		SNR_accept = 5


		if ForT == "FREQ":
			FF_fluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_peak_flux.txt", unpack = True)
			FF_rmsarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_local_rms.txt" ,unpack=True)
			FF_errfluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_err_peak_flux.txt",unpack=True)
			FF_freqarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_freq.txt", unpack=True)
			FF_slicearrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_slice.txt", unpack=True)

			FF_fluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_peak_flux.txt", unpack = True)
			FF_rmsarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_local_rms.txt" ,unpack=True)
			FF_errfluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_err_peak_flux.txt",unpack=True)
			FF_freqarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_freq.txt", unpack=True)
			FF_slicearrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_slice.txt", unpack=True)

			frac_fluxarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_flux.txt', unpack = True)
			frac_freqarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_freq.txt', unpack = True)
			frac_radiusarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_radius.txt', unpack = True)
			frac_errarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_err.txt', unpack = True)
	

			#FF_freqarrayI, FF_fluxarrayI, FF_errfluxarrayI, FF_slicearrayI, FF_rmsarrayI = negative_stokes_I_deleter(FF_fluxarrayI, FF_freqarrayI, FF_errfluxarrayI, FF_slicearrayI, FF_rmsarrayI)
			FF_errfluxarrayI = FF_yerr_bar_fixer(FF_fluxarrayI, FF_errfluxarrayI)
			FF_errfluxarrayV = FF_yerr_bar_fixer(FF_fluxarrayV, FF_errfluxarrayV)
			frac_errarray, frac_fluxarray, frac_freqarray, frac_radiusarray = frac_yerr_bar_fixer(frac_errarray, frac_fluxarray, frac_freqarray, frac_radiusarray)

			#number_of_channels = freq_chans[ind]
			#left_freq = bandwidth_low[ind]
			#right_freq = bandwidth_up[ind]
			#freqx0 = [left_freq, right_freq]
			print("mean polfrac "+name+" "+ForT+":", np.mean(frac_fluxarray))

			freqx0 = [FF_freqarrayI[0]-FF_slicearrayI[0], FF_freqarrayI[-1]+FF_slicearrayI[-1]]
			freqy0 = [0,0]
			freqy100= [100, 100]
			tickrange = np.arange(round(min(freqx0)), round(max(freqx0))+1, 5.0)
			print(tickrange)

			f,axs = plt.subplots(3, 1, sharey='row', sharex='col', facecolor='w',figsize=(18.0, 20.0),gridspec_kw={'hspace': 0})
			axI, axV, axFrac = axs

			#Plot FF_Stokes_I
			axI.plot(freqx0,freqy0,c='grey', linestyle='dashed', lw=2)
			axI.errorbar(FF_freqarrayI, FF_fluxarrayI*1000.0, xerr = FF_slicearrayI, yerr = FF_errfluxarrayI*1000, solid_capstyle='projecting', capsize=4, c = 'black', label = 'forced fit '+str(SNR_accept)+r'$\sigma$ detections', fmt='.', ms=10, lw=2)
			axI.errorbar(FF_freqarrayI, FF_rmsarrayI*1000.0, xerr = FF_slicearrayI, c = 'dodgerblue', label = 'local RMS',fmt='.', ms=10, lw=2) 

			#Plot FF_Stokes_V
			axV.plot(freqx0,freqy0,c='grey', linestyle='dashed', lw=2)
			axV.errorbar(FF_freqarrayV, FF_fluxarrayV*1000.0, xerr = FF_slicearrayV, yerr = FF_errfluxarrayV*1000, solid_capstyle='projecting', capsize=4, c = 'black', label = 'forced fit '+str(SNR_accept)+r'$\sigma$ detections', fmt='.', ms=10, lw=2)
			axV.errorbar(FF_freqarrayV, FF_rmsarrayV*1000.0, xerr = FF_slicearrayV, c = 'dodgerblue', label = 'local RMS',fmt='.', ms=10, lw=2) 

			#Plot Polfraction V/I
			axFrac.plot(freqx0,freqy0,c='grey', linestyle='dashed', ms=10, lw=2)
			axFrac.plot(freqx0,freqy100, c = 'grey', linestyle='dashdot', lw=2)
			axFrac.errorbar(frac_freqarray,frac_fluxarray*100, xerr = frac_radiusarray, yerr=frac_errarray*100, solid_capstyle='projecting', capsize=4, c = 'black', label = "Polarization Fraction (V/I)", fmt='.', ms=10, lw=2)

			for axis in ['top','bottom','left','right']:
				axI.spines[axis].set_linewidth(2)
				axV.spines[axis].set_linewidth(2)
				axFrac.spines[axis].set_linewidth(2)

			axI.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axI.tick_params(axis='both',which='major',length=8,width=1.5)
			axI.tick_params(axis='both',which='minor',length=5,width=1.5)
			axV.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axV.tick_params(axis='both',which='major',length=8,width=1.5)
			axV.tick_params(axis='both',which='minor',length=5,width=1.5)
			axFrac.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axFrac.tick_params(axis='both',which='major',length=8,width=1.5)
			axFrac.tick_params(axis='both',which='minor',length=5,width=1.5)
	
			axI.xaxis.set_ticks(tickrange)
			axV.xaxis.set_ticks(tickrange)
			axFrac.xaxis.set_ticks(tickrange)

			#axFrac.yaxis.set_ticks([0,10,20,30,40,50,60,70,80,90,100])
			#axFrac.set_yticklabels([0,'',20,'',40,'',60,'',80,'',100])

			axI.set_ylabel(r'$S_{I}$ (mJy)', fontsize = 20)
			axV.set_ylabel(r'$S_{V}$ (mJy)', fontsize = 20)
			axFrac.set_ylabel(r'$|S_{V}/S_{I}|\:(\%)$', fontsize = 20)
			axFrac.set_xlabel('Frequency (MHz)', fontsize = 20)
			#ax4.xaxis.set_label_coords(1.05, -0.235)

			#f.suptitle('Star: '+starname+', Field: '+field, fontsize=20)
			plt.minorticks_on()
			plt.savefig(triple_path+name+'_'+ForT+'_triple_30july.eps',format='eps')
			plt.close()


		if ForT == "TIME":
			FF_fluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_peak_flux.txt", unpack = True)
			FF_rmsarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_local_rms.txt" ,unpack=True)
			FF_errfluxarrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_err_peak_flux.txt",unpack=True)
			FF_timearrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_time.txt", unpack=True)
			FF_slicearrayI = np.loadtxt(forced_fitter_path+starname_forcedI+"/results/"+starname_forcedI+"_forced_slice.txt", unpack=True)

			FF_fluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_peak_flux.txt", unpack = True)
			FF_rmsarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_local_rms.txt" ,unpack=True)
			FF_errfluxarrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_err_peak_flux.txt",unpack=True)
			FF_timearrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_time.txt", unpack=True)
			FF_slicearrayV = np.loadtxt(forced_fitter_path+starname_forcedV+"/results/"+starname_forcedV+"_forced_slice.txt", unpack=True)

			frac_fluxarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_flux.txt', unpack = True)
			frac_timearray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_time.txt', unpack = True)
			frac_radiusarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_radius.txt', unpack = True)
			frac_errarray = np.loadtxt(polfractions_path+name+'_'+ForT+'/'+name+'_frac_err.txt', unpack = True)

			#FF_timearrayI, FF_fluxarrayI, FF_errfluxarrayI, FF_slicearrayI, FF_rmsarrayI = negative_stokes_I_deleter(FF_fluxarrayI, FF_timearrayI, FF_errfluxarrayI, FF_slicearrayI, FF_rmsarrayI)
			FF_errfluxarrayI = FF_yerr_bar_fixer(FF_fluxarrayI, FF_errfluxarrayI)
			FF_errfluxarrayV = FF_yerr_bar_fixer(FF_fluxarrayV, FF_errfluxarrayV)
			frac_errarray, frac_fluxarray, frac_timearray, frac_radiusarray = frac_yerr_bar_fixer(frac_errarray, frac_fluxarray, frac_timearray, frac_radiusarray)

			print("mean polfrac "+name+" "+ForT+":", np.mean(frac_fluxarray))
			timex0 = [0, (FF_timearrayI[-1]+FF_slicearrayI[-1])-(FF_timearrayI[0]-FF_slicearrayI[0])]
			timey0 = [0,0]
			timey100= [100, 100]
			tickrange = np.arange(min(timex0), round(max(timex0))+0.1 , 1.0)
			print(tickrange)
			print(len(FF_timearrayI))
			print(len(FF_fluxarrayI)) 
			print(len(FF_slicearrayI))
			print(len(FF_errfluxarrayI))
			print('timex0', timex0)
			print('lowest time:', (FF_timearrayI[0]-FF_slicearrayI[0]))
			f,axs = plt.subplots(3, 1, sharey='row', sharex='col', facecolor='w',figsize=(18.0,20.0),gridspec_kw={'hspace': 0})
			axI, axV, axFrac = axs

			#Plot FF_Stokes_I
			axI.plot(timex0,timey0,c='grey', linestyle='dashed', lw=2)
			axI.errorbar(FF_timearrayI-(FF_timearrayI[0]-FF_slicearrayI[0]), FF_fluxarrayI*1000.0, xerr = FF_slicearrayI, yerr = FF_errfluxarrayI*1000, solid_capstyle='projecting', capsize=4, c = 'black', label = 'forced fit '+str(SNR_accept)+r'$\sigma$ detections', fmt='.', ms=10, lw=2)
			axI.errorbar(FF_timearrayI-(FF_timearrayI[0]-FF_slicearrayI[0]), FF_rmsarrayI*1000.0, xerr = FF_slicearrayI, c = 'dodgerblue', label = 'local RMS',fmt='.', ms=10, lw=2) 

			#Plot FF_Stokes_V
			axV.plot(timex0,timey0,c='grey', linestyle='dashed', lw=2)
			axV.errorbar(FF_timearrayV-(FF_timearrayI[0]-FF_slicearrayI[0]), FF_fluxarrayV*1000.0, xerr = FF_slicearrayV, yerr = FF_errfluxarrayV*1000, solid_capstyle='projecting', capsize=4, c = 'black', label = 'forced fit '+str(SNR_accept)+r'$\sigma$ detections', fmt='.', ms=10, lw=2)
			axV.errorbar(FF_timearrayV-(FF_timearrayI[0]-FF_slicearrayI[0]), FF_rmsarrayV*1000.0, xerr = FF_slicearrayV, c = 'dodgerblue', label = 'local RMS',fmt='.', ms=10, lw=2) 

			#Plot Polfraction V/I
			axFrac.plot(timex0,timey0,c='grey', linestyle='dashed', ms=10, lw=2)
			axFrac.plot(timex0,timey100, c = 'grey', linestyle='dashdot', lw=2)
			axFrac.errorbar(frac_timearray-(FF_timearrayI[0]-FF_slicearrayI[0]),frac_fluxarray*100, xerr = frac_radiusarray, yerr=frac_errarray*100, solid_capstyle='projecting', capsize=4, c = 'black', label = "Polarization Fraction (V/I)", fmt='.', ms=10, lw=2)

			for axis in ['top','bottom','left','right']:
				axI.spines[axis].set_linewidth(2)
				axV.spines[axis].set_linewidth(2)
				axFrac.spines[axis].set_linewidth(2)

			axI.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axI.tick_params(axis='both',which='major',length=8,width=1.5)
			axI.tick_params(axis='both',which='minor',length=5,width=1.5)
			axV.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axV.tick_params(axis='both',which='major',length=8,width=1.5)
			axV.tick_params(axis='both',which='minor',length=5,width=1.5)
			axFrac.tick_params(axis='both',which='both',labelsize=15,direction='in')
			axFrac.tick_params(axis='both',which='major',length=8,width=1.5)
			axFrac.tick_params(axis='both',which='minor',length=5,width=1.5)
			
			axI.xaxis.set_ticks(tickrange)
			axV.xaxis.set_ticks(tickrange)
			axFrac.xaxis.set_ticks(tickrange)

			#axFrac.yaxis.set_ticks([0,10,20,30,40,50,60,70,80,90,100])
			#axFrac.set_yticklabels([0,'',20,'',40,'',60,'',80,'',100])

			axI.set_ylabel(r'$S_{I}$ (mJy)', fontsize = 20)
			axV.set_ylabel(r'$S_{V}$ (mJy)', fontsize = 20)
			axFrac.set_ylabel(r'$|S_{V}/S_{I}|\:(\%)$', fontsize = 20)
			axFrac.set_xlabel('Time (hrs)', fontsize = 20)
			#ax4.xaxis.set_label_coords(1.05, -0.235)

			#f.suptitle('Star: '+starname+', Field: '+field, fontsize=20)
			plt.minorticks_on()
			plt.savefig(triple_path+name+'_'+ForT+'_triple_30july.eps',format='eps')
			plt.close()



FFI_FFV_polfrac_plotter()			





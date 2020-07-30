# WABIFAT
WABIFAT (WSClean, Aegean, Bane Interaction to create fluxdensity versus Frequency And Time plots) is a software capable of creating spectra and light curves of radio observations done telescopes such as LOFAR.

This software creates spectra and light curves from input .MS files resulting from radio observations of e.g LOFAR. 

It is not user friendly at the moment. I will be working on that. 
- The WABIFAT_FINAL file contains the code which adaptively bins channels to find detections or non-detections. 
- FF_and_more is the forced fitting mechanism which applies the priorized fitting option of Aegean to the detections and non-detections resulting from WABIFAT. 
  Also it calculates the circular polarization fraction of the resulting spectra and light curves. Followed by plotting the Stokes I, Stokes V, fraction plots.
- plot.py creates the input fits file necessary for WABIFAT and FF_and_more. And contains some functions on plotting the results from WABIFAT. 

Kas Veken

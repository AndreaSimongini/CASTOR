# CASTOR

[THIS GITHUB PAGE IS STILL WORK IN PROGRESS]

Building templates and reconstructing parameters for Core Collapse SNe. 

You will need: 
- Python3 and all the basics packages such as numpy, scipy, ...
- The package **george** (https://george.readthedocs.io/en/latest/) 

## Data

Every photometric and spectrometric data for each supernova in the training set is taken from Open Supernova Catalog https://github.com/astrocatalogs/supernovae and WISeREP https://www.wiserep.org/ respectively. You will find the data for every SN out of the training set in the /Data/ directory. 

### data_lightcurves

Every light curve point is contained in a file.dat that **must** follow this rules in order to be read by CASTOR with no problems: 
- The file.dat is named after the SN (i.e. SN2014G.dat)
- The first four columns must be Time (in MJD), magnitude, error of the magnitude and filter.
- The file must be in the right directory which is Localpath/Data/data_lightcurves/

### data_spectra 

Every spectral point is contained in a file.dat that **must** follow this rules in order to be read by CASTOR with no problems: 
- The file.dat is named after the SN + the epoch (in MJD) at which was observed (i.e. ASASSN14jb_56957.39907407408.dat).
- In the case of two spectra observed in the same day, you must differentiate the names of the two files adding some minutes (i.e. iPTF13bvn_56460.0.dat and iPTF13bvn_56460.914675925924.dat)
- The first four columns must be Wavalength (in AA), error of the wavelength (in AA), flux (in erg/s/cm^2/A) and error of the flux (in erg/s/cm^2/A)
- The file must be in the right directory which is Localpath/Data/data_spectra/SN_name/






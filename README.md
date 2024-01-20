# CASTOR

[THIS GITHUB PAGE IS STILL WORK IN PROGRESS]

Building templates and reconstructing parameters for Core Collapse SNe. 

You will need: 
- Python3 and all the basics packages such as numpy, scipy, ...
- The package **george** (https://george.readthedocs.io/en/latest/)

## Reference

In this directoy you will find two .xlcs files:
- SN_names.xlcs containes a column with the names of the SN contained in the training set.
- bandpass.xlcs containes the needed information for basics filters, such as effective wavelength, fwhm and transformation from Jhonson to AB system.
  

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

## How to use CASTOR step by step: 

**STEP 0**: clone this github repository. The enviromental variable LOCALPATH should automatically set the right directory. The *Data* and *Reference* folders must be downloaded and placed in the same directory as the CASTOR code. 

**Step 1**: prepare your own Supernova. 

You should create in the same directory a *Application/SN_name* (where SN_name is the name of your new SN i.e. Application/SN2015ap/) folder where to put your new light curve data respecting this nomenclature: 
- The file must be in .dat format
- The name of the file must be the SN_name (i.e. SN2015ap.dat)
- The first four columns of the file must be time (in MJD), magnitude, error and filter
- The filter must follow the same nomenclature as in the bandpass.xlcs file. - You have to set the new_name variable to the name of your new supernova.

  **Step 2**: Building templates  






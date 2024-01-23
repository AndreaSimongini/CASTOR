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
- The filter must follow the same nomenclature as in the bandpass.xlcs file.
- You have to set the new_name variable to the name of your new supernova.

**Step 2**: get into CASTOR  

If you followed every step correctly so far, the software will now operate autonomously. What to know about every section and what is your role. 
- Section 0: only functions for reading and analyze data. 
- Section 1: preliminary studies. If the training set has not changed, leave the "execute" binary variable to False setting.
- Section 2: comparison between YOUR SN and the training set. You must have followed **Step 1** info.
- Section 3: creation of spectral templates. This is auto, but you can choose either to execute it or not.
- Section 4: parameter estimation. This is partially automatic, but you can choos either to execute it or not.

**Step 3**: velocity fitting 

The entire software is *almost* completely authomatic, except for the velocity fitting routine, which requires the hands and the eyese of an active user. 
Here the instructions to follow in order to correctly estimate velocity, redshift and so distance and spectral class. First the list of the questions you'll be asked, then the list of the answers and **how** you are supposed to answer in order to continue. 

1. First, a spectrum is shown to you. You have to decide whether you see or not see both a P-Cygni profile and an Absorption line. If you don't see any, the second spectrum will be shown, and so on. 
2. Then, the selected spectrum is shown with lines above it. You have to decide which lines (one for P-Cygni and one for absorption) you want to study.
3. Then, you have to define the interval containing first the P-Cygni and then the absorption.
4. Finally, the code will automatically analyze every spectra on the basis of your choices. 

Now the answers: note that you must respect also the spaces in between words if you have to put more than one. 

1. YES / NO (it works fine also with yes / no)
2. Line 1, Line 2 (where Line 1 is the P-Cygni line and Line 2 is the absorption line. You **must** follow the same nomenclature)
3. a, b (it will add +a to the left and +b to the right. Either way you can reduce the space adding a - before the number)

An example: 
1. Can you see at least one P-Cygni and one emission line? [Yes]
2. Select the P-Cygni and the emission lines [He Ia, Fe IIa]
3. Please insert values to enlarge or reduce the interval [100, -15]
4. Please insert values to enlarge or reduce the interval [0, 0]


Once every question is answered the software will work automatically and give you the results. 





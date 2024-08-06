# CASTOR (**Core collApse Supernovae parameTers estimatOR**) 

CASTOR is a full open access software that takes as input the *light curves* of a core collapse supernova and estimates the parameters belonging to the ejecta, to the progenitor star and to the event itself. We define the **case-of-study** supernova as the object we want to analyse, **training set** as the list of published supernovae used as reference and **reference** supernova as the most resembling supernova to our case-of-study out of the training set. To fully recover the parametric map of the case-of-study supernova, we apply the following procedures: 
1. We compare the light curves of the case-of-study supernova with the light curves of every SN from the training set.
2. By applying a chi-squared test, we select the reference supernova.
3. Using the light curves of the case-of-study supernova and the spectra of the reference supernova, we create synthetic spectra.
4. Using the synthetic light curves and spectra of the case-of-study supernova we estimate parameters.

There are some general physical assumptions we made to estimate the parameters:
- Spherical symmetry
- Mass conservation
- Canonical partition of energy between photons and neutrinos
- Modified-black-body emission

The list of filter allowed in CASTOR is: 
- Bessell system (U, B, V, R, I)
- SDSS and modified SDSS systems (u, g, r, i, z and u', g', r', i', z')
- Johnson-Glass system (Y, J, H, K, Ks)
- SWIFT/UVOT system (w1, w2, m2)

## 1 - Data

All data needed for the analysis are available in the `Training_Set`. Data are divided into:
- `data_spectra`: containing spectral data of each supernova
- `data_lightcurves`: containing photometric data of each supernova

Every light curve file is in the form `name.dat` (i.e. SN2015ap.dat) and contains data points divided in four columns: 
- time: in MJD units
- magnitude: in the AB system
- error: in the AB system
- filter: in the common nomenclature (i.e. B, u, K, w1)

Every spectrum file is in the form `name_MJD.dat` (i.e. SN2015ap_57273.8433.dat) and contains data points divided in two columns: 
- wavelength: in $ \AA $ units
- flux: in erg/s/cm $^2$ /\AA units 








## Ingredients 

In order to correctly use CASTOR you will need a few ingredients. Here, in this repository, you will find the entire training set dataset. This contains: 
- A Training_Set.xlsx file with a list of names, types, redshift and other ancillary information.
- A data_lightcurves folder containing every lightcurve of each supernova in the training set
- A data_spectra folder containing every spectra of each supernova in the training set

The light curve files are all of the form **name.dat** containing four columns: time (in MJD), apparent magnitude, error and filter. 
The spectra files are all of the form **epochMJD.dat** containing three columns: wavelength (in A), flux (in erg/s/cm2/A) and error. 

What **you** will need to add to these data is a **name.dat** file (i.e. *SN2015ap.dat*) containing all your photometric points with correct units and following the column scheme of time, magnitude, error and filter. Note that:
- time: has to be in MJD units
- magnitude: has to be in AB system
- error: has to be in AB system
- filter: has to follow the previous nomenclature (i.e. B, u, K, w1)



`castor_source.py`
`calibration.py`
`execute_castor.ipynb`


## Instructions 

Once every directory is ready for usage, you will simply need to execute the **castor_source.py** file available in this repository. We prepared also a **castor_usage.ipynb** notebook to help. The information that need to be given as input by the user are: 
- name of the supernova (be sure that the name is exactly the same as defined in the name.dat file)
- path of the training set data
- path of the name.dat file
- path to the output directory

Once every input has been defined, you can execute the notebook and CASTOR will automatically produce the spectral templates and estimate the parameters. If you are using the **castor_source.py** you will be asked to set some important information regarding the P-Cygni profile fitting process that will need a user selection choice. 

CASTOR will automatically generate the following outputs: 
- a comparison.png image containing the comparison of the case-of-study and the reference supernovae light curves
- a templates.png image which shows the synthetic spectra reconstruction
- a velocity.png image which shows the evolution of expansion velocity
- a parameters.txt text file with all the parameters with their relative uncertainty. 

## Debugging 

There are a few ways to see if your analysis went fine or there are some problems. The first thing to check is the templates.png file: if the synthetic spectra are looking bad, there is probably a problem in the name.dat file that can be of the form: 
- incorrect normalization of time
- incorrect definition of magnitudes
- incorrect definition of filters
- bad sampling

If the spectral templates look fine, you should check the velocity.png image: the velocity should have a rather normal behaviour: if you see something weird it's probably due to a bad selection of P-Cygni lines to study. An important check is also the distance: since many parameters directly depend on the estimation of distance, if their results look odd, it's probably due to an incorrect estimation of the distance parameter. 



## Aknowledgments








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
- wavelength: in Å units
- flux: in erg/s/cm2/Å units

Every information regarding the training set supernovae is contained in the `Training_Set.xlsx` file, which needs to be in the same directory as the data files. 

## 2 - Scripts 

The scripts needed for this analysis are all deposited in the `Scripts` directory you find in this repository. This contains: 
- the `castor_source.py` file containing the entire source code of CASTOR for the user-dependent analysis
- the `execute_castor.ipynb` notebook containing the commands and the input needed for the analysis
- the `calibration.py` file showing how the calibration of the training set is performed
- the `auto_castor_source.ipynb` file containing the entire source code of CASTOR for the user-independent analysis 

## 3 - User-dependent analysis  

### 2.0 - Your data 

CASTOR can be employed in two ways: 
- Analysis of a supernova for which only photometry is available (Simongini et al. 2024)
- Analysis of a supernova for which photometry and spectroscopy are available (LST Collaboration 2024, in prep)

In the first case, you will only need to prepare your `name.dat` file containing the light curves of the case-of-study supernova following the same scheme as the training set data (four columns, time, magnitude, error and filter). 

In the second case, other than the light curve file, you will need to prepare also the `name_MJD.dat` files containing your spectra following the same scheme as the training set data (two columns, wavelegnth and flux). Then, you have to put your data together with those of the training set in a folder called as your supernova. Finally, a little hack is needed in the software: you replace the `final_name` variable with the name of your source, which will be then used as both case-of-study and reference supernovae. 

### 2.1 - User input

The information that need to be given as input by the user are: 
- name of the supernova (be sure that the name is exactly the same as defined in the name.dat file)
- path of the training set data
- path of the name.dat file
- path to the output directory

### 2.2 - Output 

CASTOR will automatically generate the following outputs: 
- a comparison.png image containing the comparison of the case-of-study and the reference supernovae light curves
- a templates.png image which shows the synthetic spectra reconstruction
- a spectrum_MJD.dat file for each of the synthetic spectra 
- a velocity.png image which shows the evolution of expansion velocity
- a results.txt text file with all the parameters with their relative uncertainty. 

### 2.3 Profile fitting

For the user-dependent analysis, you will be asked to provide some easy inputs, helping the software with your own eye for the P-Cygni fitting analysis. In order, you will be asked to select the epoch in which to start the analysis. Then, every available line will be plotted in the selected spectrum: you will need to select one P-Cygni and one absorption line (the list of lines can be found in the head of the code). For both lines you will be asked to select the wavelength range in which to fit the profile, adding respectively to the left and to the right. Finally, you will be asked to select the class: red lines are Hydrogen lines and blue lines are Helium lines. 

These are the Q/A you will need to answer: 
1. Can you see at least one P-Cygni and one emission line? `yes / no`
2. Select the P-Cygni and the emission lines. The answer should be comma-space separated i.e. `Helium Ia, Hydrogen a`
3. Please insert values to enlarge or reduce the interval i.e. `10, -35`
4. If the interval is fine, please type `0, 0`
5. Which class is this? Choose between II, IIP, IIn, IIb, Ib, Ib/c, Ic.



## Debugging 

There are a few ways to see if your analysis went fine or there are some problems. The first thing to check is the templates.png file: if the synthetic spectra are looking bad, there is probably a problem in the name.dat file that can be of the form: 
- incorrect normalization of time
- incorrect definition of magnitudes
- incorrect definition of filters
- bad sampling

If the spectral templates look fine, you should check the velocity.png image: the velocity should have a rather normal behaviour: if you see something weird it's probably due to a bad selection of P-Cygni lines to study. An important check is also the distance: since many parameters directly depend on the estimation of distance, if their results look odd, it's probably due to an incorrect estimation of the distance parameter. 



## Aknowledgments

Please if you used CASTOR for your analysis or your scientific results cite the following paper: 
--Inserire paper quando sarà pubblicato-- 
And link this github page in the footnotes 







# CASTOR (**Core collApse Supernovae parameTers estimatOR**) 

CASTOR is a open access software that takes as input the *light curves* of a core collapse supernova and estimates the parameters belonging to the ejecta, to the progenitor star and to the event itself. We define the **case-of-study** supernova as the object we want to analyse, **training set** as the list of published supernovae used as reference and **reference** supernova as the most resembling supernova to our case-of-study out of the training set. To fully recover the parametric map of the case-of-study supernova, we apply the following procedures, extensively described in [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911) : 
1. We compare the light curves of the case-of-study supernova with the light curves of every SN from the training set.
2. By applying a chi-squared test, we select the reference supernova.
3. Using the light curves of the case-of-study supernova and the spectra of the reference supernova, we create synthetic spectra.
4. Using the synthetic light curves and spectra of the case-of-study supernova we estimate parameters.

There are some general physical assumptions we made to estimate the parameters:
- Spherical symmetry
- Mass conservation
- Canonical partition of energy between photons and neutrinos
- Modified-black-body emission
- Perfect adiabaticity at the time of maximum luminosity

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
- (for the future) the `auto_castor_source.ipynb` file containing the entire source code of CASTOR for the user-independent analysis 

## 3 - User-dependent analysis  

### 3.1 - Input data 

CASTOR can be employed in two ways: 
- Analysis of a supernova for which only photometry is available (described in [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911)). 
- Analysis of a supernova for which photometry and spectroscopy are available (described in LST Collaboration 2024, in prep)

In the first case, you will only need to prepare your `name.dat` file containing the light curves of the case-of-study supernova following the same scheme as the training set data (four columns, time, magnitude, error and filter). 

In the second case, other than the light curve file, you will need to prepare also the `name_MJD.dat` files containing your spectra following the same scheme as the training set data (two columns, wavelegnth and flux). Then, you have to put your data together with those of the training set in a folder called as your supernova. Finally, a little hack is needed in the software: you replace the `final_name` variable with the name of your source, which will be then used as both case-of-study and reference supernovae. 

### 3.2 - Other input

The information that need to be given as input by the user are: 
- name of the supernova (be sure that the name is exactly the same as defined in the name.dat file)
- path of the training set data
- path of the name.dat file
- path to the output directory

### 3.3 - Output 

CASTOR will automatically generate the following outputs: 
- a comparison.png image containing the comparison of the case-of-study and the reference supernovae light curves
- a templates.png image which shows the synthetic spectra reconstruction
- a spectrum_MJD.dat file for each of the synthetic spectra 
- a velocity.png image which shows the evolution of expansion velocity
- a results.txt text file with all the parameters with their relative uncertainty. 

### 3.4 Profile fitting

For the user-dependent analysis, you will be asked to provide some easy inputs, helping the software with your own eye for the P-Cygni fitting analysis. In order, you will be asked to select the epoch in which to start the analysis. Then, every available line will be plotted in the selected spectrum: you will need to select one P-Cygni and one absorption line (the list of lines can be found in the head of the code). For both lines you will be asked to select the wavelength range in which to fit the profile, adding respectively to the left and to the right. Finally, you will be asked to select the class: red lines are Hydrogen lines and blue lines are Helium lines. 

These are the questions you will need to answer: 
1. Can you see at least one P-Cygni and one emission line? `yes / no`
2. Select the P-Cygni and the emission lines. The answer should be comma-space separated i.e. `Helium Ia, Hydrogen a`
3. Please insert values to enlarge or reduce the interval i.e. `10, -35`
4. If the interval is fine, please type `0, 0`
5. Which class is this? Choose between II, IIP, IIn, IIb, Ib, Ib/c, Ic.


## 4. Debugging 

As always, things can go pretty wrong and you will need to re-do your analysis. During the development and testing of CASTOR I faced several problems and learnt how to solve them. Here a few suggestions regarding the most common problems I met. 

1. The first thing you need to check is the `templates.png` file: are the spectra looking reasonable? If not, there can be some problems in your input data:
     - incorrect normalization of time (is the time column in MJD format?)
     - incorrect definition of magnitudes (is the magnitude column in AB system? Are magnitudes apparent and not absolute?)
     - filters (are the filters in the allowed `filterlist` above and with the correct nomenclature?)
2. Another important check is the estimate of the `t0` and the `tmax`. If the sampling of your data is bad, or they are collected far from the explosion, the estimate of `t0` may be wrong, thus you need to set it manually. The same problem can occur with the `tmax`: if your data are collected far from the maximum, then it is impossible to determine it with CASTOR. Another problem may occur if the light curves show two peaks with the second brighter than the first: CASTOR will select the second peak for the estimate of `tmax'. Note that many problems you may encounter are simply related to a bad estimation of these two parameters which depend strongly on how the light curves are taken.
3. Check the `velocity.png` plot and the `distance` value: are they reasonable and close to what you were expecting? If not, you may need to execute again the `parameter_estimation` routine and select different lines or different values for the range or different epochs. You may need to repeat this step several times for both features (P-Cygni and absorption) until you find the correct one that gives you the correct results, fine tuning the specifics. 4.
6. Check the `gp.png` plot. If the GPs failed in fitting the light curves, there can be several problems in templates building and parameter estimation. This problem may occur if:
   - the `t0` value is estimated badly.
   - data points are too scattered in time.
   - data points are too scattered in magnitude following a non-physical behaviour. 

## 5. Future updates for version 1.1
* safer selection of real spectral points from the reference supernova in the building_templates class
* update gaussian process fitting of light curves
* definition of a n_cores variable
* accurate scaling of real spectral points for templates reconstruction
* implementation of a user-independent pipeline of analysis
* dynamic definition of R$^2$ limit
* faster estimation of mass of nichel 

## References (main)

[1] [Aryan et al. 2021, MNRAS, 505, 2530](10.1093/mnras/stab1379). 
[2] [Branch & WHeeler 2017, Astronomy and Astrophysics Library](10.1007/978-3-662-55054-0). 
[3] [Filippenko 1997, ARA&A, 35, 309](10.1086/309659). 
[4] [Meza et al. 2019, A&A, 629, A57](10.1051/0004-6361/201834972). 
[5] [Nakamura et al. 2016, MNRAS, 461, 3296](10.1093/mnras/stw1453). 
[6] [Vincenzi et al. 2019, MNRAS, 489, 5802](10.1093/mnras/stz2448). 
[7] [Yaron & Gal-Yam, Publ. Astron. Soc. Pac., 124, 668](10.1086/666656).

## Aknowledgments

To acknowledge CASTOR in a publication, please cite [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911). 

Thank you!










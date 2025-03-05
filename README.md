# CASTOR 2.0
(**Core collApse Supernovae parameTers estimatOR**) 

CASTOR is a open access software that takes as input the *light curves* of a core collapse supernova and estimates the parameters belonging to the ejecta, to the progenitor star and to the event itself. We define the **case-of-study** supernova as the object we want to analyse, **training set** as the list of published supernovae used as reference and **reference** supernova as the most resembling supernova to our case-of-study out of the training set. To fully recover the parametric map of the case-of-study supernova, we apply the following procedures, extensively described in [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911) : 
1. We compare the light curves of the case-of-study supernova with the light curves of every SN from the training set.
2. By applying a chi-square test, we select the reference supernova.
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
- SDSS and modified SDSS systems (u, g, r, i, z)
- Johnson-Glass system (Y, J, H, K, Ks)
- SWIFT/UVOT system (w1, w2, m2)

## 1 - Data

All data needed for the analysis are available in the `Training_Set`. Data are divided into:
- `data_spectra`: containing spectral data of each supernova
- `data_lightcurves`: containing photometric data of each supernova

Note that for storage limit only we divided the spectral data into `data_spectra_1`, `data_spectra_2`, `data_spectra_3`. Make sure to put all files in a single folder called `data_spectra`. 

Every light curve file is in the form `name.dat` (i.e. SN1987A.dat) and contains data points divided in four columns: 
- time: in MJD units
- magnitude: in the AB system
- error: in the AB system
- filter: in the common nomenclature (i.e. B, u, K, w1)

Every spectrum file is in the form `name_MJD.dat` (i.e. SN1987A_0.0.dat) and contains data points divided in two columns: 
- wavelength: in Å units
- flux: in erg/s/cm2/Å units

Every light curve and spectrum is already scaled with respect to the time of explosion of the object. 

Every information regarding the training set supernovae is contained in the `Training_Set.xlsx` file, which needs to be in the same directory as the data files. 

## 2 - Scripts 

CASTOR v2.0 uses only 1 script for the entire analysis, which can be found in the `Scripts` directory with the name: 
- `castor2.0.py`
This file can be easily executed from terminal previa being in the same directory of all the other files. 

## 3 - Analaysis  

### 3.1 - Input data 

CASTOR can be employed in two ways: 
- Analysis of a supernova for which only photometry is available (described in [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911)). 
- Analysis of a supernova for which photometry and spectroscopy are available (described in LST Collaboration 2024, in prep)

In the first case, you will only need to prepare your `name.dat` file containing the light curves of the case-of-study supernova following the same scheme as the training set data (four columns, time, magnitude, error and filter). 

In the second case, other than the light curve file, you will need to prepare also the `name_MJD.dat` files containing your spectra following the same scheme as the training set data (two columns, wavelegnth and flux). Then, you have to put your data together with those of the training set in a folder called as your supernova. 

An example of the standard directory scheme can be found in the folder `Analysis`

### 3.2 - Other input

The user will be directly asked to insert some inputs to facilitate the analysis and have better results. The user can decide not to insert any of the parameters:
- Name of the case-of-study supernova (must be the same nomenclature used in defining the data files)
- Time of explosion in MJD.
- Spectral availability (yes or no)
- Redshift
- Distance in Mpc
- Extinction 

### 3.3 - Output 

CASTOR will automatically create a `output_sn_name` folder containing all plots and results: 
- comparison.png: a figure containing the comparison with the reference supernova.
- light_curves.png : a figure containing the GP interpolated light curves.
- luminosity.png: a figure containing bolometric and pseudo luminosity.
- photosphere.png: a figure containing temperature and radius of photosphere.
- spectra.png: a figure containing synthetic spectra
- velocity.png: a figure containign the time evolution of line velocity.
- results.txt: a file text containing all estimated parameters. 


## 4. Debugging 

As always, things can go pretty wrong and you will need to re-do your analysis. During the development and testing of CASTOR I faced several problems and learnt how to solve them. Here a few suggestions regarding the most common problems I met. 

1. The first thing you need to check is the `spectra.png` file: are the spectra looking reasonable? If not, there can be some problems in your input data:
     - incorrect normalization of time (is the time column in MJD format?)
     - incorrect definition of magnitudes (is the magnitude column in AB system? Are magnitudes apparent and not absolute?)
     - filters (are the filters in the allowed `filterlist` above and with the correct nomenclature?)
2. Another important check is the estimate of the `t0` and the `tmax`. If the sampling of your data is bad, or they are collected far from the explosion, the estimate of `t0` may be wrong, thus you need to set it manually. The same problem can occur with the `tmax`: if your data are collected far from the maximum, then it is impossible to determine it with CASTOR. Another problem may occur if the light curves show two peaks with the second brighter than the first: CASTOR will select the second peak for the estimate of `tmax'. Note that many problems you may encounter are simply related to a bad estimation of these two parameters which depend strongly on how the light curves are taken.
3. Check the `velocity.png` plot and the `distance` value: are they reasonable and close to what you were expecting? If not, you may need to execute again the `parameter_estimation` routine and select different lines or different values for the range or different epochs. You may need to repeat this step several times for both features (P-Cygni and absorption) until you find the correct one that gives you the correct results, fine tuning the specifics. 4.
6. Check the `light_curves.png` plot. If the GPs failed in fitting the light curves, there can be several problems in templates building and parameter estimation. This problem may occur if:
   - the `t0` value is estimated badly.
   - data points are too scattered in time.
   - data points are too scattered in magnitude following a non-physical behaviour.
7. It's always a good idea to insert redshift, distance and extinction instead of letting the software estimating them automatically. It should be fine but better not to risk. 

## 5. What's new
We made some changes in the usage and in the software itself with respect to CASTOR v1.0. We highlight here the most important. A better explanation is given in `Updates2.0.pdf`
* New supernovae added to the training set (a total of 150 supernovae, references in the Updates file). 
* light curves and spectra of the training set are now scaled wrt their time of explosion. 
* A new calibration has been performed.
* Velocity is now estimated automatically.
* The user can now insert automatically redshift, distance, time of explosion and extinction.
* The user can now choose to use only the observed spectra of the case-of-study source.
* 



## Future application 
* LST Collaboration, in prep: the analysis of SN2024bch.
* Simongini, Ragosta et al., in prep: the analysis of Vera C. Rubin simulations
* Analysis of new CCSNe events 



## References (main)

1. [Aryan et al. 2021, MNRAS, 505, 2530](10.1093/mnras/stab1379).
2. [Branch & WHeeler 2017, Astronomy and Astrophysics Library](10.1007/978-3-662-55054-0).
3. [Filippenko 1997, ARA&A, 35, 309](10.1086/309659).
4. [Meza et al. 2019, A&A, 629, A57](10.1051/0004-6361/201834972).
5. [Nakamura et al. 2016, MNRAS, 461, 3296](10.1093/mnras/stw1453).
6. [Vincenzi et al. 2019, MNRAS, 489, 5802](10.1093/mnras/stz2448).
7. [Yaron & Gal-Yam, Publ. Astron. Soc. Pac., 124, 668](10.1086/666656).

## Aknowledgments

To acknowledge CASTOR in a publication, please cite [Simongini et al. 2024](https://doi.org/10.1093/mnras/stae1911). 
The current version of the software is published in Zenodo. Please cite also [Simongini et al. 2025](https://doi.org/10.5281/zenodo.14761308)

Thank you!










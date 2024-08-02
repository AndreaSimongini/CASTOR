# CASTOR (**Core collApse Supernovae parameTers estimatOR**) 

Welcome to the CASTOR repository!

## Introduction 

CASTOR is an open access software that takes as input light curves (and spectra) of a core collapse supernova and estimates the parameters of the event, of the photosphere and of the progenitor star. We define *case-of-study supernova* the object we want to analyze, *training set* the list of published supernovae we use as a reference and *reference supernova* the most resembling supernova to our case-of-study. To fully recover the parametric map of the case-of-study supernova, we apply the following procedure: 
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

## Ingredients 

In order to correctly use CASTOR you will need a few ingredients. Here, in this repository, you will find the entire training set dataset. This contains: 
- A Training_Set.xlsx file with a list of names, types, redshift and other ancillary information.
- A data_lightcurves folder containing every lightcurve of each supernova in the training set
- A data_spectra folder containing every spectra of each supernova in the training set

The light curve files are all of the form **name.dat** containing four columns: time (in MJD), apparent magnitude, error and filter. 
The spectra files are all of the form **epochMJD.dat** containing three columns: wavelength (in A), flux (in erg/s/cm2/A) and error. 

What **you** will need to add to these date is a **name.dat** file (i.e. *SN2015ap.dat*) containing all your data with correct units and following the column scheme of time, magnitude, error and filter. 

--WORK IN PROGRESS-- 
## Instructions 

## Debugging 

## Aknowledgments








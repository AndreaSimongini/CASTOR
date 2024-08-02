# CASTOR (**Core collApse Supernovae parameTers estimatOR**) 

Welcome to the CASTOR repository!

## Introduction 

CASTOR is an open access software that takes as input light curves (and spectra) of a core collapse supernova and estimates the parameters of the event, of the photosphere and of the progenitor star. We define \textit{case-of-study supernova} the object we want to analyze, \textit{training set} the list of published supernovae we use as a reference and \textit{reference supernova} the most resembling supernova to our case-of-study. To fully recover the parametric map of the case-of-study supernova, we apply the following procedure: 
1. We compare the light curves of the case-of-study supernova with the light curves of every SN from the training set.
2. By applying a chi-squared test, we select the reference supernova.
3. Using the light curves of the case-of-study supernova and the spectra of the reference supernova, we create synthetic spectra.
4. Using the synthetic light curves and spectra of the case-of-study supernova we estimate parameters.
There are some general physical assumptions we made to estimate the parameters:
- Spherical symmetry
- Mass conservation
- Canonical partition of energy between photons and neutrinos
- Modified-black-body emission

--WORK IN PROGRESS-- 
## Ingredients 

In order to correctly use CASTOR you will need a few ingredients. Here, in this repository, you will find the entire training set dataset. This contains: 
- A Training_Set.xlsx file with a list of names, types, redshift and other ancillary information.
- A data_lightcurves folder containing every lightcurve of each supernova in the training set
- A data_spectra folder containing every spectra of each supernova in the training set

The light curve files are all of the form \textbf{name.dat} containing three columns: time (in MJD), apparent magnitude and error. 
The spectra files are all of the form \textbf{epochMJD.dat} containing three columns: wavelength (in \AA), flux (in erg/s/cm$^2$/\AA) and error. 

## Instructions 

## Debugging 

## Aknowledgments








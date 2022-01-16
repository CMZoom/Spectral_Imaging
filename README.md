# Spectral Imaging Scripts
This repository contains the various codes used in the process of imaging and analysing all current CMZoom spectral line data.
I am by no means an experienced programmer and these codes were written and added to over several years, so many of them
are likely inefficient and could be improved upon greatly.

I will list the general process of running these scripts to produce the figures in the paper here.

## Step 1: Image Production

### Step 1a: Imaging ([CASA_tclean_test.py](https://github.com/CMZoom/Spectral_Imaging/blob/ff1b1caecc1666326ca32a3213d5aab70a68644f/CASA_tclean_test.py))
- The primary code for imaging the data from measurement sets. To remedy issues with gaps in coverage
between SMA spectral windows, this script begins by concatenating the individual uvfits files into one .ms file per sideband.
- Performs continuum subtraction by masking out major CO transitions in cube.
- For each sideband, creates a dirty image of size [100,100,100] surrounding the phasecenter. This dirty image is used to
calculate the noise threshold for the full clean.
- If region is <1000 in (x,y), cleans a full sideband using this 5x this noise threshold. Otherwise, cleans a subcube
surrounding the key transitions.

### Step 1b: Galactic Projection ([beam_fix.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/beam_fix.py))
- A number of cubes have issues with varying beamsizes throughout. While not a major concern in most cases, some cubes had
severa variation. This script reprojects all cubes into galactic coordinates, masks out channels that are at least 30% larger
than the average, and convolves rest of cube to a uniform beam size.

### Step 1c: Cube Trimming ([fits_cubes.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/fits_cubes.py))
- For smaller regions that had entire sidebands imaged, this script produces subcubes of the same transitions as the previous
step makes for the larger cubes.

## Step 2: Spectra

### Step 2a: Spectrum Extraction ([all_spectra.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/all_spectra.py))
- For every transition, this script extracts a spectrum for each leaf in the dendrogram catalog and adds it to a fits file
in the form [leaf_ID, [intensity], [velocity]].

### Step 2b: Spectrum Fitting ([scouse_fitting.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/scouse_fitting.py))
CURRENT CODE SEEMS TO BE OLD VERSION, NO PANDAS INVOLVED. WILL UPDATE
- Has two modes; single and full.
  - Single takes a leaf ID and a transition and updates the scouse fitting.
  - Full allows user to fit all spectra for a given transition.
- Uses scousepy to fit all spectra in survey and compile the fit values in a csv file using pandas.

## Step 3: Moment Maps

### Step 3a: Moment Map Production ([moment_maps_new.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/moment_maps_new.py))
-Uses binary dilation to produce moment maps

### Step 3b: Moment Map Figures ([moment_figures.py](https://github.com/CMZoom/Spectral_Imaging/blob/b11087996ebfe51cbc9fa7ebd79db43aef31b646/moment_figures.py))
- Takes moment maps generated from previous steps and generates complete moment map figures for each transition and region.

## Step 4: Analysis

### Step 4a (DC_CMZoom_fits_analysis.py)
This is currently an unorganised mess as this developed slowly over the entire analysis of the paper and I didn't start breaking things up until near the end. Currently, it:
- fixes small issues with data
- produces latex tables
- quality control
- corrects mislabeled peaks (especially the methanol line) NEEDS IMPROVED, CURRENTLY ONLY WORKS WITH PRIOR KNOWLEDGE OF PANDAS INDICES.
- uses splatalogue to identify peaks that aren't key transition
- makes all histogram figures in paper
- makes correlation matrix
- identifies single velocity for each core
#### TO-DO:
  - Split above steps up into individual, uniquely responsible scripts.
  - Improve data fixing so it requires no prior knowledge and remains consistent (or at least identifies cores that need fixing in comments
  - Remove extra, unneeded analysis.

### Step 4b (Larsons.py)
- For each unique core, with mass and radius identified in catalog, calculates alpha parameter.

### Step 4c (pressure.py)
- For each core, using SF tracers in Perry's paper, redoes alpha parameter analysis with pressure factored in.

 ## TO-DO:
- clean and comment codes for ease of use
- update older codes
- add moment ratio code

# Spectral Imaging Scripts
This repository contains the various codes used in the process of imaging and analysing all current CMZoom spectral line data.
I am by no means an experienced programmer and these codes were written and added to over several years, so many of them
are likely inefficient and could be improved upon greatly.

I will list the general process of running these scripts to produce the figures in the paper here.

## Step 1: Imaging (CASA_tclean_test.py)
- The primary code for imaging the data from measurement sets. To remedy issues with gaps in coverage
between SMA spectral windows, this script begins by concatenating the individual uvfits files into one .ms file per sideband.
- Performs continuum subtraction by masking out major CO transitions in cube.
- For each sideband, creates a dirty image of size [100,100,100] surrounding the phasecenter. This dirty image is used to
calculate the noise threshold for the full clean.
- If region is <1000 in (x,y), cleans a full sideband using this 5x this noise threshold. Otherwise, cleans a subcube
surrounding the key transitions.

### Step 1b: Galactic Projection (beam_fix.py)
- A number of cubes have issues with varying beamsizes throughout. While not a major concern in most cases, some cubes had
severa variation. This script reprojects all cubes into galactic coordinates, masks out channels that are at least 30% larger
than the average, and convolves rest of cube to a uniform beam size.

### Step 1c: Cube Trimming (fits_cubes.py)
- For smaller regions that had entire sidebands imaged, this script produces subcubes of the same transitions as the previous
step makes for the larger cubes.

## Step 2: Spectra

### Step 2a: Spectrum Extraction (all_spectra.py)
- For every transition, this script extracts a spectrum for each leaf in the dendrogram catalog and adds it to a fits file
in the form [leaf_ID, [intensity], [velocity]].

### Step 2b: Spectrum Fitting (scouse_fitting.py)
- Has two modes; single and full.
  - Single takes a leaf ID and a transition and updates the scouse fitting.
  - Full allows user to fit all spectra for a given transition.
-

TO-DO:
- clean and comment codes for ease of use

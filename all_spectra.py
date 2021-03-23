from astropy.io import fits
from astropy.table import Table
from astropy.stats import mad_std
import astropy.units as u
import astropy.constants as constants
import numpy as np
import matplotlib.pyplot as plt
from pyspeckit import Spectrum, spectrum
import os
from astroquery.splatalogue import Splatalogue
from astropy import constants
import aplpy
from astropy.coordinates import FK5
from astropy.coordinates import SkyCoord
import sys
import matplotlib.gridspec as gridspec
import glob
import math as m

def line_spectrum(data, original_header, line, leaf, mask):
	spectral_mask = np.where(mask == leaf)
	
	line_vel = np.linspace(original_header['crval3'], original_header['crval3']+(original_header['naxis3']*original_header['cdelt3']), original_header['naxis3'])

	if len(str(int(line_vel[0]))) > 3:
		line_vel = line_vel/1e3

	line_vel = line_vel.tolist()
	spectrum = [np.nanmean(data[i][spectral_mask]) for i in range(len(line_vel))]

	if len(spectrum) < 369:
		n = 369 - len(spectrum)
		spectrum = n * [0] + spectrum
		line_vel = np.pad(line_vel, (n,0), 'linear_ramp', end_values=(line_vel[0]+n*(line_vel[0]-line_vel[1])))
	if len(spectrum) > 369:
		n = len(spectrum)-369
		spectrum = spectrum[:(-1*n)]
		line_vel = line_vel[:(-1*n)]

	spectrum = np.array(spectrum)
	NaN = np.isnan(spectrum)
	spectrum[NaN] = 0

	if np.any(spectrum):
		spectrum = spectrum.tolist()
		l = len(spectrum)
		spectrum.insert(0, leaf)
		for x in line_vel:
			spectrum.append(x)
		savefile = '../'+line+'_spectrum.fits'
		if os.path.exists(savefile):
			print 'file found, adding to it.'
			hdu = fits.open(savefile)[0].data
			if hdu.ndim == 1:
				hdu = np.vstack((hdu,spectrum))
			elif spectrum[0] in hdu[:,0]:
				hdu[np.where(hdu[:,0] == spectrum[0])] = spectrum
			else:			
				hdu = np.vstack((hdu,spectrum))
			print np.shape(hdu)
			hdu1 = fits.PrimaryHDU(hdu)
			hdu1.writeto(savefile, overwrite=True)
		else:
			print 'No file found, making it.'
			hdu = fits.PrimaryHDU(spectrum)
			hdu.writeto(savefile)
	else:
		print 'Empty spectrum'

dendro_mask = fits.open('dendrogram_mask_pruned_rms3e6_k14_dv3_dd1_dp17_pp6_pm2_gal_10-23-19.fits')
mask = dendro_mask[0].data
mask_reference = aplpy.FITSFigure('dendrogram_mask_pruned_rms3e6_k14_dv3_dd1_dp17_pp6_pm2_gal_10-23-19.fits', convention='wells')

for sources in glob.glob('G0.014+0.021_CASA/'):
	sourcename = sources[:-6]
	print sourcename
	if sourcename == 'G359.863-0.069' or sourcename == 'G0.699-0.028' or sourcename == 'G0.619+0.012' or sourcename == 'G359.948-0.052' or sourcename == 'G359.863-0.069':
		continue
	os.chdir(sources)

	if '+' in sourcename:
		lon = float(sourcename.split('+')[0][1:])
		lat = float(sourcename.split('+')[1])
	else:
		lon = float(sourcename.split('-')[0][1:])
		lat = -1*float(sourcename.split('-')[1])

	for filename in glob.glob('*GHz.fits'):
		line = filename[len(sourcename)+1:-5]
		print sourcename,line
		cube = fits.open(filename)
		data = cube[0].data
		print np.shape(data)
		header = cube[0].header
		
		ImageSize = np.shape(header['NAXIS1'])
		lon_range = [float(header['CRVAL1']-(header['CRPIX1']*header['CDELT1'])), float(header['CRVAL1']+((header['NAXIS1']-header['CRPIX1'])*header['CDELT1']))]
		lat_range = [float(header['CRVAL2']-(header['CRPIX2']*header['CDELT2'])), float(header['CRVAL2']+((header['NAXIS2']-header['CRPIX2'])*header['CDELT2']))]

		a = mask_reference.world2pixel(lon_range[0]*u.degree, lat_range[0]*u.degree)
		b = mask_reference.world2pixel(lon_range[1]*u.degree, lat_range[1]*u.degree)

		Mask = mask[int(a[1]):int(b[1]),int(a[0]):int(b[0])]

		leaves = np.unique(Mask)
		leaves = leaves[leaves != -1]
		leaves = leaves.tolist()

		for leaf in leaves:
			line_spectrum(data, header, line, leaf, Mask)

	os.chdir('../')	
		

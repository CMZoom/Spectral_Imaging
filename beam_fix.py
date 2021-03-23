from spectral_cube import VaryingResolutionSpectralCube
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import astropy.units as u
import os
from reproject import reproject_interp
import shutil
from radio_beam import Beams
from spectral_cube.cube_utils import beams_to_bintable
from spectral_cube import LazyMask
from astropy import wcs
import glob

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
rc('text',usetex=True)
#rc('text.latex',preamble=r'\usepackage{bbding}')
rc('patch',antialiased=False)
from matplotlib import rcParams
rcParams['font.size'] = 15.5
rcParams['xtick.major.size'] = 6
rcParams['ytick.major.size'] = 6
rcParams['xtick.minor.size'] = 3
rcParams['ytick.minor.size'] = 3

sourcename = sys.argv[1]

from guppy import hpy
h = hpy()

os.chdir(sourcename+'_CASA/')

beam = []
freq = []
conv = []
for filename in glob.glob('*.deep.test.pbcor.fits'):
	lsb_cube = VaryingResolutionSpectralCube.read(filename)
	lsb_freq = lsb_cube.spectral_axis
	freq.append(lsb_freq)
	print beams_to_bintable(lsb_cube.beams)
	try:
		beamtable = beams_to_bintable(lsb_cube.beams)

		lsb_beam = lsb_cube.beams.sr
		beam.append(lsb_beam)
		lsb_cube.allow_huge_operations=True
	
		if os.path.isfile(sourcename+'.continuum.fits'):
			outhdr = fits.open(sourcename+'.continuum.fits')[0].header
			hdu = fits.open(filename)
			header = hdu[0].header
	
			outhdr['NAXIS'] = hdu[0].header['NAXIS']
		
			#outhdr.insert('NAXIS2', 'NAXIS3', after=True)
			outhdr['NAXIS3'] = hdu[0].header['NAXIS3']
			outhdr['CTYPE3'] = hdu[0].header['CTYPE3']
			outhdr['CDELT3'] = hdu[0].header['CDELT3']
			outhdr['CRVAL3'] = hdu[0].header['CRVAL3']
			outhdr['CRPIX3'] = hdu[0].header['CRPIX3']
			outhdr['CASAMBM'] = hdu[0].header['CASAMBM']
#		del outhdr['BMAJ']
#		del outhdr['BMIN']
#		del outhdr['BPA']
			outhdr['RESTFRQ'] = hdu[0].header['RESTFRQ']
			outhdr['RADESYS'] = hdu[0].header['RADESYS']
			outhdr['EQUINOX'] = hdu[0].header['EQUINOX']
			outhdr['SPECSYS'] = hdu[0].header['SPECSYS']

			shape = (outhdr['NAXIS3'], outhdr['NAXIS2'], outhdr['NAXIS1'])
			print outhdr
			outarray = np.memmap(filename='output.np', mode='w+', shape=shape, dtype='float32')

			rslt = reproject_interp(hdu, outhdr, output_array=outarray, return_footprint=False, independent_celestial_slices=True)
			comp_mask = LazyMask(np.isfinite, data=outarray, wcs=wcs.WCS(outhdr))
			new_cube = VaryingResolutionSpectralCube(data=outarray, header=outhdr, wcs=wcs.WCS(outhdr), mask=comp_mask, beam_table=beamtable.data, beam_threshold=lsb_cube.beam_threshold)
			new_cube.allow_huge_operations=True
			masked_lsb_cube = new_cube.mask_out_bad_beams(0.3, criteria=['sr','major','minor'])
			save_file = 'convolved.galactic'
		else:
			masked_lsb_cube = lsb_cube.mask_out_bad_beams(0.3, criteria=['sr','major','minor'])
			save_file = '.convolved'
		masked_lsb_cube.allow_huge_operations = True
		lsb_common_beam = masked_lsb_cube.beams[masked_lsb_cube._goodbeams_mask].common_beam(tolerance=1e-5)	
		convolved_lsb_cube = masked_lsb_cube.convolve_to(lsb_common_beam)
		lsb_convolve = convolved_lsb_cube.beam.sr.value
		conv.append(lsb_convolve)
		#convolved_lsb_cube.write(filename[:-10]+save_file+'.pbcor.fits', format='fits', overwrite=True)
	except	AttributeError:
		if os.path.isfile(sourcename+'.continuum.fits'):
			outhdr = fits.open(sourcename+'.continuum.fits')[0].header
			hdu = fits.open(filename)
			header = hdu[0].header
	
			outhdr['NAXIS'] = hdu[0].header['NAXIS']
		
			#outhdr.insert('NAXIS2', 'NAXIS3', after=True)
			outhdr['NAXIS3'] = hdu[0].header['NAXIS3']
			outhdr['CTYPE3'] = hdu[0].header['CTYPE3']
			outhdr['CDELT3'] = hdu[0].header['CDELT3']
			outhdr['CRVAL3'] = hdu[0].header['CRVAL3']
			outhdr['CRPIX3'] = hdu[0].header['CRPIX3']
#		del outhdr['BMAJ']
#		del outhdr['BMIN']
#		del outhdr['BPA']
			outhdr['RESTFRQ'] = hdu[0].header['RESTFRQ']
			outhdr['RADESYS'] = hdu[0].header['RADESYS']
			outhdr['EQUINOX'] = hdu[0].header['EQUINOX']
			outhdr['SPECSYS'] = hdu[0].header['SPECSYS']

			shape = (outhdr['NAXIS3'], outhdr['NAXIS2'], outhdr['NAXIS1'])
			outarray = np.memmap(filename='output.np', mode='w+', shape=shape, dtype='float32')

			rslt = reproject_interp(hdu, outhdr, output_array=outarray, return_footprint=False, independent_celestial_slices=True)
			comp_mask = LazyMask(np.isfinite, data=outarray, wcs=wcs.WCS(outhdr))
			new_cube = SpectralCube(data=outarray, header=outhdr, wcs=wcs.WCS(outhdr), mask=comp_mask)
			new_cube.allow_huge_operations=True
			save_file = 'convolved.galactic'
		else:
			save_file = 'convolved'
		#new_cube.write(filename[:-10]+save_file+'.pbcor.fits', format='fits', overwrite=True)
		

else:
	print 'No Lower Sideband Image present'
'''
if os.path.isfile(sourcename+'_CASA/'+sourcename+'.asic.usb.multiscale.deep.test.pbcor.fits'):
	usb_cube = VaryingResolutionSpectralCube.read(sourcename+'_CASA/'+sourcename+'.asic.usb.multiscale.deep.test.pbcor.fits')
	usb_freq = usb_cube.spectral_axis

	beamtable = beams_to_bintable(usb_cube.beams)

	usb_beam = usb_cube.beams.sr
	usb_cube.allow_huge_operations=True
	if os.path.isfile(sourcename+'_CASA/'+sourcename+'.continuum.fits'):
		outhdr = fits.open(sourcename+'_CASA/'+sourcename+'.continuum.fits')[0].header
		hdu = fits.open(sourcename+'_CASA/'+sourcename+'.asic.usb.multiscale.deep.test.pbcor.fits')
		header = hdu[0].header
	
		outhdr['NAXIS'] = hdu[0].header['NAXIS']
#		outhdr.insert('NAXIS2', 'NAXIS3', after=True)
		outhdr['NAXIS3'] = hdu[0].header['NAXIS3']
		outhdr['CTYPE3'] = hdu[0].header['CTYPE3']
		outhdr['CDELT3'] = hdu[0].header['CDELT3']
		outhdr['CRVAL3'] = hdu[0].header['CRVAL3']
		outhdr['CRPIX3'] = hdu[0].header['CRPIX3']
		outhdr['CASAMBM'] = hdu[0].header['CASAMBM']
		del outhdr['BMAJ']
		del outhdr['BMIN']
		del outhdr['BPA']
		outhdr['RESTFRQ'] = hdu[0].header['RESTFRQ']
		outhdr['RADESYS'] = hdu[0].header['RADESYS']
		outhdr['EQUINOX'] = hdu[0].header['EQUINOX']
		outhdr['SPECSYS'] = hdu[0].header['SPECSYS']

		shape = (outhdr['NAXIS3'], outhdr['NAXIS2'], outhdr['NAXIS1'])
		
		outarray = np.memmap(filename='output.np', mode='w+', shape=shape, dtype='float32')

		rslt = reproject_interp(hdu, outhdr, output_array=outarray, return_footprint=False, independent_celestial_slices=True)
		comp_mask = LazyMask(np.isfinite, data=outarray, wcs=wcs.WCS(outhdr))
		new_cube = VaryingResolutionSpectralCube(data=outarray, header=outhdr, wcs=wcs.WCS(outhdr), mask=comp_mask, beam_table=beamtable.data, beam_threshold=usb_cube.beam_threshold)
		new_cube.allow_huge_operations=True
		masked_usb_cube = new_cube.mask_out_bad_beams(0.3, criteria=['sr','major','minor'])
		save_file = '.convolved.galactic'
	else:
		masked_usb_cube = usb_cube.mask_out_bad_beams(0.3, criteria=['sr','major','minor'])
		save_file = '.convolved'

	usb_common_beam = masked_usb_cube.beams[masked_usb_cube._goodbeams_mask].common_beam(tolerance=1e-5)	
	convolved_usb_cube = masked_usb_cube.convolve_to(usb_common_beam)
	usb_convolve = convolved_usb_cube.beam.sr.value
	convolved_usb_cube.write(sourcename+'_CASA/'+sourcename+'.asic.usb.multiscale.deep'+save_file+'.pbcor.fits', format='fits', overwrite=True)

else:
	print 'No Upper Sideband Image present'

usb_cube = None
gal_usb_cube = None
masked_usb_cube = None
convolved_usb_cube = None
'''

fig = plt.figure(figsize=(16,12))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
ax1.plot(freq[0]/1e9,beam[0]/1e-9)
ax1.axhline(conv[0]/1e-9,linestyle='--',color='black')
ax1.set_xlabel('Frequency (GHz)', fontsize=20)
ax1.set_ylabel(r'Beam Area (10$^{-9}$Sr)', fontsize=20)
#ax1.annotate(sourcename, (0.1, 0.9), xycoords='figure fraction', fontsize=20)
ax1.tick_params(axis='both')

ax2.plot(freq[1]/1e9,beam[1]/1e-9)
ax2.axhline(conv[1]/1e-9,linestyle='--',color='black')
ax2.set_xlabel('Frequency (GHz)', fontsize=20)
ax2.set_ylabel(r'Beam Area (10$^{-9}$Sr)', fontsize=20)
#ax2.annotate(sourcename, (0.1, 0.9), xycoords='figure fraction', fontsize=20)
ax2.tick_params(axis='both')

'''
if 'lsb_freq' in locals():
	ax1.plot(lsb_freq/1e9, lsb_beam/1e-9)
	ax1.axhline(lsb_convolve/1e-9, color='black')
	ax1.set_xlabel('Frequency (GHz)', fontsize=20)
	ax1.set_ylabel(r'Beam area (10$^{-9}$Sr)', fontsize=20)
	ax1.annotate(sourcename, (0.1, 0.9), xycoords='figure fraction', fontsize=20)
	ax1.tick_params(axis='both', labelsize=20)
if 'usb_freq' in locals():
	ax2.plot(usb_freq/1e9, usb_beam/1e-9)
	ax2.axhline(usb_convolve/1e-9, color='black')
	ax2.set_xlabel('Frequency (GHz)', fontsize=20)
	ax2.tick_params(axis='both', labelsize=20)
'''

plt.savefig(sourcename+'_Beam_Area.pdf')




import sys
import os
import moments
import exp_mask
import aplpy
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import astropy.units as au
import astropy.stats as stats
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
import glob
from matplotlib import cm
import additionalmods as admods
import rms_comp
from scipy.ndimage.morphology import binary_dilation, binary_closing
import pandas as pd

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
#rc('text',usetex=True)
#rc('text.latex',preamble=r'\usepackage{bbding}')
rc('patch',antialiased=False)
from matplotlib import rcParams
rcParams['font.size'] = 9
rcParams['xtick.major.size'] = 6
rcParams['ytick.major.size'] = 6
rcParams['xtick.minor.size'] = 3
rcParams['ytick.minor.size'] = 3

vels = pd.read_csv('Velocity_vals.csv')
vels['Leaf_ID'] = vels.apply(lambda row: row.id[:-1], axis = 1)

for filename in glob.iglob('G1.038*_CASA/'):
	sourcename = filename.split('_')[0]
	source_vel = np.array(vels[vels['Leaf_ID'] == sourcename]['v_peak'])
	if len(source_vel) > 0:
		velmin = np.min(source_vel)-20
		velmax = np.max(source_vel)+20
		cuberange = 'modded'
	else:
		cuberange = 'full'

	print(sourcename)
	
	os.chdir(filename)
	rms_check = 0

#	try:
#		lsb_rms = fits.open(sourcename+'_lsb_rms.fits')
#		usb_rms = fits.open(sourcename+'_usb_rms.fits')
#		rms_check = 1
#	except IOError:
#		rms_check = 0
		
	for input_cube in glob.iglob('*GHz.fits'):
		data = fits.open(input_cube)
		header = data[0].header
		rest_freq = header['RESTFRQ']
		if (rest_freq / 1000) < 1:
			header['RESTFRQ'] = header['RESTFRQ']*1e9
			fits.writeto(filename=input_cube, data=data[0].data, header=header, overwrite=True)
		else:
			print "Frequency is in Hz"

		velo = np.arange(header['CRVAL3'], header['CRVAL3']+(header['NAXIS3']*header['CDELT3']), header['CDELT3'])
		if header['CUNIT3'] == 'm/s':
			velo = velo / 1000
			header['CUNIT3'] = 'km/s'
		print header['NAXIS3']
		velo = velo.astype(int)

		if cuberange == 'modded':
			mom_velo = [velmin,velmax]
		else:
			mom_velo = [velo[-91],velo[90]]
		print mom_velo
		cube = moments.get_cube(input_cube, x_units=(au.km/au.s), y_units='Jy / beam')
		#cube = exp_mask.get_delchans(cube)

#		if rms_check == 1 and rest_freq < 230:
#			rms = lsb_rms
#			rms_data = rms[0].data
#			rms_header = rms[0].header
#		elif rms_check == 1 and rest_freq > 230:
#			rms = usb_rms
#			rms_data = rms[0].data
#			rms_header = rms[0].header
#		else:
#			rms_velo = [[velo[-1],velo[-51]],[velo[50],velo[0]]]
#			rms_cube = moments.get_rms(cube, rms_velo)
#			rms_data = rms_cube.data
#			rms_header = rms_cube.header

#		print np.shape(rms), np.shape(cube)

		if '12CO' in input_cube:
			rms_data = rms_comp.get_rmsrob(cube[-50:])
			rms = fits.PrimaryHDU(rms_data.data, rms_data.header)
		else:
			rms_1, rms_2 = rms_comp.get_rmsrob(cube[:50]), rms_comp.get_rmsrob(cube[-50:])
			rms_data = np.nanmean([rms_1.data, rms_2.data], axis=0)
			print np.shape(rms_data)
			rms = fits.PrimaryHDU(rms_data, rms_1.header)

		print np.shape(rms.data)
    		bmaj = cube.header['BMAJ']
    		bmin = cube.header['BMIN']
    		pix = np.absolute(cube.header['CDELT1'])
    		bmajp = bmaj/pix
   		bminp = bmin/pix
    		bareap = np.pi*bmajp*bminp
    		rad = np.floor(bmajp/2)
    		structure = admods.get_circmask(radius=rad)
		mask_l_, _ = moments.get_threshmask(cube, rms.data, thresh=3)
		mask_l = mask_l_.include()

		for i in range(6):
			j = 10 - i
    			mask_h_, _ = moments.get_threshmask(cube, rms.data, thresh=j)
    			mask_h = mask_h_.include()
			mask_d = binary_dilation(mask_h, iterations=-1, mask=mask_l)
    			mask_dc = binary_closing(mask_d, structure=structure, iterations=1)
    			mask_dcc = admods.get_prunemask(mask_dc, thresh=3*bareap)
			sigma=str(j)+'$\sigma$'
			if j == 5 and np.nansum(mask_dcc.flatten()) == 0:
				cube_masked = cube
				nchan = admods.get_nchanbool(cube_masked)
				sigma='unmasked'
			if np.nansum(mask_dcc.flatten()) != 0:
				nchan = admods.get_nchanbool(mask_dcc)
				cube_masked = cube.with_mask(mask_dcc)
				break

		
		fig = plt.figure(figsize=(12,12))
		ax = aplpy.FITSFigure(rms, figure=fig, subplot=[0.05,0.52,0.4,0.4])
		ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_axis_label_text('rms (K)')
		ax.add_label(0.23, 0.95, 'RMS', relative=True, weight='bold')
		ax.ticks.hide()
		ax.tick_labels.hide()
		ax.axis_labels.hide()

		levels = [3,10,15,20,30,40,50,60,70,80,90,120,150]

		mom = moments.get_momentmaps(cube_masked, mom_velocity=mom_velo)
		mom0_hdu = mom['mom0'].hdu
		mom1_hdu = mom['mom1'].hdu
		mom2_hdu = mom['sigma'].hdu
		
		if np.nansum(mom0_hdu.data) == 0:
			print 'check'
			cube_masked = cube
			nchan = admods.get_nchanbool(cube_masked)
			sigma='unmasked'
			mom = moments.get_momentmaps(cube_masked, mom_velocity=mom_velo)
			mom0_hdu = mom['mom0'].hdu

		mom0_hdr = mom0_hdu.header
		mom0_hdr = cube_masked[0].header
		mom0_hdr.append(('THRESH', sigma))
		cdelt = np.absolute(cube.header['CDELT3'])
   		mom0err_hdu = moments.get_mom0err(nchan, rms, cdelt)
		
    		s2n_data = np.array(mom['mom0'])/mom0err_hdu.data
    		ids = np.isnan(s2n_data)
    		s2n_data[ids] = 0
    		s2n_hdu = fits.PrimaryHDU(s2n_data, rms.header)

    		maxs2n_data = np.array(mom['max'])/rms.data
   		maxs2n_hdu = fits.PrimaryHDU(maxs2n_data, rms.header)
		ids = np.where(rms.data<0)
    		rms.data[ids] = np.nan
    		mom0_hdu.data[ids] = np.nan
    		mom0err_hdu.data[ids] = np.nan
    		s2n_hdu.data[ids] = np.nan
		
		maxs2n_hdu.data[np.isnan(maxs2n_hdu.data)] = 0
		s2n_hdu.data[np.isnan(s2n_hdu.data)] = 0
		print np.nanmean(mom2_hdu.data)
		print np.all(s2n_hdu.data == 0)
		if ~np.all(s2n_hdu.data == 0):
			#rms_out = input_cube.replace('.fits', '_rms.fits')
    			mom0_out = input_cube.replace('.fits', '_mom0.fits')
			mom1_out = input_cube.replace('.fits', '_mom1.fits')
			mom2_out = input_cube.replace('.fits', '_mom2.fits')
    			mom0err_out = input_cube.replace('.fits', '_mom0err.fits')
    			s2n_out = input_cube.replace('.fits', '_mom0s2n.fits')
    			#rms.writeto(rms_out, overwrite=True)
			fits.writeto(mom0_out, data=mom0_hdu.data, header=mom0_hdr, overwrite=True)
			fits.writeto(mom1_out, data=mom1_hdu.data, header=mom0_hdr, overwrite=True)
			fits.writeto(mom2_out, data=mom2_hdu.data, header=mom0_hdr, overwrite=True)
    			mom0err_hdu.writeto(mom0err_out, overwrite=True)
    			s2n_hdu.writeto(s2n_out, overwrite=True)

		if sigma == 'unmasked':
			fig.suptitle(input_cube+r' ('+str(sigma)+')')
			ax = aplpy.FITSFigure(s2n_hdu, figure=fig, subplot=[0.05,0.05,0.4,0.4])
			ax.show_colorscale(cmap=cm.gist_heat_r, interpolation='nearest')
			ax.add_colorbar()
			ax.add_label(0.23, 0.95, 'Signal/Noise', relative=True, weight='bold')
		else:
			fig.suptitle(input_cube+r' ('+str(sigma)+')')
			ax = aplpy.FITSFigure(s2n_hdu, figure=fig, subplot=[0.05,0.05,0.4,0.4])
			ax.show_colorscale(cmap=cm.gist_heat_r, interpolation='nearest', vmin=0, pmax=100)
			ax.add_colorbar()
			ax.add_label(0.23, 0.95, 'Signal/Noise', relative=True, weight='bold')

		ax.colorbar.set_axis_label_text(r'Signal/Noise')
		cont = fits.open(sourcename+'.continuum.fits')
		cont_data = cont[0].data
		cont_lvl = [np.nanpercentile(cont_data, 98),np.nanpercentile(cont_data, 99)]
		
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvl, colors='blue', linewidth=0.01)
		ax.set_nan_color('white')
		ax.ticks.set_xspacing(0.01)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.set_xformat('dd.dd')
		ax.tick_labels.set_yformat('dd.dd')

		pmin = np.nanpercentile(mom['max'].hdu.data, 5)
    		pmax = np.nanpercentile(mom['max'].hdu.data, 99.9)
		
		ax = aplpy.FITSFigure(mom['max'].hdu, figure=fig, subplot=[0.52,0.05,0.4,0.4])
		ax.add_label(0.23, 0.95, 'Maximum', relative=True, weight='bold')

		if sigma == 'unmasked':
			ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest')
		else:
			ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest', vmin=pmin, vmax=pmax)

		ax.add_colorbar()
		ax.colorbar.set_axis_label_text(r'Maximum (K)')
		ax.ticks.hide()
		ax.tick_labels.hide()
		ax.axis_labels.hide()

		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvl, colors='blue', linewidth=0.1)
	
		pmin = np.nanpercentile(mom0_hdu.data, 5)
  		pmax = np.nanpercentile(mom0_hdu.data, 99.9)
		
		ax = aplpy.FITSFigure(mom0_hdu, figure=fig, subplot=[0.52,0.52,0.4,0.4])
		ax.add_label(0.23, 0.95, 'Moment 0', relative=True, weight='bold')
		
		if sigma == 'unmasked':
			ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest')
		else:
			ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest', vmin=pmin, vmax=pmax)
		ax.add_colorbar()
		ax.colorbar.set_axis_label_text('mom0 (K kms$^{-1}$)')
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvl, colors='blue', linewidth=0.1)
		ax.ticks.hide()
		ax.tick_labels.hide()
		ax.axis_labels.hide()
		ax = fig.axes[0]
		for tick in ax.get_yticklabels():
			tick.set_rotation(90)
		if ~np.all(s2n_hdu.data == 0):
			fig.savefig(input_cube+'_amoment_maps.pdf', dpi=500, bbox_inches='tight')
		del rms
		del mom0_hdu
		del mom
		del s2n_hdu
	os.chdir('../')

plt.show()
#fig.savefig('rms.pdf', dpi=500, bbox_inches='tight')

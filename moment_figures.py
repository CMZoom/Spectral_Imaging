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
import glob
from matplotlib import cm
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
rc('text',usetex=True)
rc('text.latex',preamble=r'\usepackage{bbding}')
rc('patch',antialiased=False)
from matplotlib import rcParams
rcParams['font.size'] = 15.5
rcParams['xtick.major.size'] = 6
rcParams['ytick.major.size'] = 6
rcParams['xtick.minor.size'] = 3
rcParams['ytick.minor.size'] = 3
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

def tick_lvls(filename):
	data = fits.open(filename)[0].data
	ticks = [np.nanpercentile(data, 0), np.nanpercentile(data, 25), np.nanpercentile(data, 50),
		np.nanpercentile(data, 75), np.nanpercentile(data, 100)]
	return ticks
sourcename = 'G359.889-0.093'

if '359' in sourcename:
	nx = 0.04
else:
	nx = 0.02

cont = fits.open(sourcename+'.continuum.fits')[0].data
cont_lvls = [np.nanpercentile(cont, 99)]

fig = plt.figure(figsize=(20,15))
#fig.suptitle(sourcename)
ax = aplpy.FITSFigure(sourcename+'.continuum.fits', figure=fig, subplot=[0.06,0.67,0.2,0.3])
ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
ax.show_colorscale(cmap=cm.gist_heat_r,interpolation='nearest')
ax.add_colorbar()
ax.colorbar.set_axis_label_text('Jy/beam')
ax.colorbar.set_location('right')
ax.colorbar.set_pad(0)
ax.add_beam()
ax.add_scalebar(0.0084)
ax.scalebar.set_linewidth(2)
ax.scalebar.set_corner('bottom right')
ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
ax.scalebar.set_font_weight('bold')
ax.add_label(0.3, 0.95, 'CONTINUUM', relative=True, weight='bold')
ax.ticks.show()
ax.ticks.set_xspacing(0.03)
ax.ticks.set_yspacing(0.01)
ax.tick_labels.set_xformat('dd.dd')
ax.tick_labels.set_yformat('dd.dd')
ax.ticks.set_color('black')
ax.tick_labels.set_xposition('top')
ax.axis_labels.set_xposition('top')
ax.tick_labels.set_yposition('left')
ax.axis_labels.set_yposition('left')
ax.axis_labels.show()

moment = '2'

if moment == '0':
    cmap = cm.gist_heat_r
    label = r'Jy/beam'
elif moment == '1':
    cmap = cm.RdBu
    label = r'K km s$^{-1}$'
else:
    cmap = cm.plasma
    label = r'K km s$^{-1}$'

try:
    filename = fits.open(sourcename+'.12CO.230.5GHz_mom'+moment+'.fits')
    thresh = filename[0].header['thresh']
    if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
        print('No 12CO Emission')
        fig.text(0.32,0.82,'No 12CO Emission', fontsize=14, fontweight='bold')
    else:
        ax = aplpy.FITSFigure(sourcename+'.12CO.230.5GHz_mom'+moment+'.fits', figure=fig, subplot=[0.30,0.67,0.2,0.3])
        ax.show_contour(sourcename+'.continuum.fits', layer='cont',levels=cont_lvls, colors='black',linewidths=1, smooth=1)
        ax.show_contour(sourcename+'.12CO.230.5GHz_mom0s2n.fits', layer='s2n',levels=[50, 100, 150], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
        ax.show_colorscale(cmap=cmap,interpolation='nearest')
        ax.add_colorbar()
        ax.colorbar.set_location('right')
        ax.colorbar.set_pad(0)
        ax.add_beam()
        ax.add_scalebar(0.0084)
        ax.scalebar.set_linewidth(2)
        ax.scalebar.set_corner('bottom right')
        ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
        ax.scalebar.set_font_weight('bold')
        ax.add_label(0.3, 0.95, '12CO - '+thresh, color='black', relative=True, weight='bold')
        ax.ticks.show()
        ax.ticks.set_xspacing(nx)
        ax.ticks.set_yspacing(0.01)
        ax.ticks.set_color('black')
        ax.tick_labels.hide()
        ax.axis_labels.hide()

except IOError:
    print('non')
    print('No 12CO Data.')
    fig.text(0.32,0.82,'No 12CO Data',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.13CO.220.4GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No 13CO Emission')
		fig.text(0.55, 0.82,'No 13CO Emission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.13CO.220.4GHz_mom'+moment+'.fits', figure=fig, subplot=[0.53,0.67,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.13CO.220.4GHz_mom0s2n.fits', layer='s2n',levels=[20, 50, 75], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.3, 0.95, '13CO - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No 13CO Data.')
	fig.text(0.55,0.82,'No 13CO Data',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.C18O.219.6GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No C18O Emission')
		fig.text(0.78,0.82,'No C18O Emission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.C18O.219.6GHz_mom'+moment+'.fits', figure=fig, subplot=[0.76,0.67,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.C18O.219.6GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_axis_label_text(label)
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.4, 0.95, 'C18O - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No C18O Data.')
	fig.text(0.78,0.82,'No C18O Data',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.H2CO.218.2GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No H2CO (218.2 GHz) Emission')
		fig.text(0.08,0.49,'No H2CO (218.2 GHz)\nEmission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.H2CO.218.2GHz_mom'+moment+'.fits', figure=fig, subplot=[0.06,0.34,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.H2CO.218.2GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.5, 0.95, 'H2CO (218.2GHz) - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No H2CO (218.2 GHz) Data.')
	fig.text(0.08,0.49,'No H2CO (218.2 GHz)\nData',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.H2CO.218.5GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No H2CO (218.5 GHz) Emission')
		fig.text(0.32,0.49,'No H2CO (218.5 GHz)\nEmission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.H2CO.218.5GHz_mom'+moment+'.fits', figure=fig, subplot=[0.30,0.34,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.H2CO.218.5GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.5, 0.95, 'H2CO (218.5GHz) - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No H2CO (218.5 GHz) Data.')
	fig.text(0.32,0.49,'No H2CO (218.5 GHz)\nData',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.H2CO.218.8GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No H2CO (218.8 GHz) Emission')
		fig.text(0.55,0.49,'No H2CO (218.8 GHz)\nEmission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.H2CO.218.8GHz_mom'+moment+'.fits', figure=fig, subplot=[0.53,0.34,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.H2CO.218.8GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.5, 0.95, 'H2CO (218.8GHz) - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No H2CO (218.8 GHz) Data')
	fig.text(0.55,0.49,'No H2CO (218.8 GHz)\nData',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.OCS.218.9GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No OCS (218.9 GHz) Emission')
		fig.text(0.78,0.49,'No OCS (218.9 GHz)\nEmission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.OCS.218.9GHz_mom'+moment+'.fits', figure=fig, subplot=[0.76,0.34,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.OCS.218.9GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_axis_label_text(label)
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.5, 0.95, 'OCS (218.9GHz) - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No OCS (218.9 GHz) Data.')
	fig.text(0.78,0.49,'No OCS (218.9 GHz)\nData',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.OCS.231.1GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No OCS (231.1 GHz) Emission')
		fig.text(0.08,0.16,'No OCS (231.1 GHz)\nEmission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.OCS.231.1GHz_mom'+moment+'.fits', figure=fig, subplot=[0.06,0.01,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.OCS.231.1GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.5, 0.95, 'OCS (231.1GHz) - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
		
except IOError:
	print('No OCS (231.1 GHz) Data.')
	fig.text(0.08,0.16,'No OCS (231.1 GHz)\nData',fontsize=14, fontweight='bold')

try:	
	filename = fits.open(sourcename+'.SiO.217.1GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No SiO Emission')
		fig.text(0.32,0.16,'No SiO Emission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.SiO.217.1GHz_mom'+moment+'.fits', figure=fig, subplot=[0.30,0.01,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.SiO.217.1GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.28, 0.95, 'SiO - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No SiO Data.')
	fig.text(0.32,0.16,'No SiO Data',fontsize=14, fontweight='bold')

try:
	filename = fits.open(sourcename+'.SO.219.9GHz_mom'+moment+'.fits')
	thresh = filename[0].header['thresh']
	if thresh == 'unmasked' and moment == '1' or thresh == 'unmasked' and moment == '2':
		print('No SO Emission')
		fig.text(0.55,0.16,'No SO Emission', fontsize=14, fontweight='bold')
	else:
		ax = aplpy.FITSFigure(sourcename+'.SO.219.9GHz_mom'+moment+'.fits', figure=fig, subplot=[0.53,0.01,0.2,0.3])
		ax.show_contour(sourcename+'.continuum.fits', levels=cont_lvls, colors='black', linewidths=1, smooth=1)
		ax.show_contour(sourcename+'.SO.219.9GHz_mom0s2n.fits', layer='s2n',levels=[5, 7, 10], colors='grey', linestyles='dashed',linewidths=1, smooth=1)
		ax.show_colorscale(cmap=cmap,interpolation='nearest')
		ax.add_colorbar()
		ax.colorbar.set_axis_label_text(label)
		ax.colorbar.set_location('right')
		ax.colorbar.set_pad(0)
		ax.add_beam()
		ax.add_scalebar(0.0084)
		ax.scalebar.set_linewidth(2)
		ax.scalebar.set_corner('bottom right')
		ax.scalebar.set_label(str(int(0.0084*3600*0.1/3))+' pc')
		ax.scalebar.set_font_weight('bold')
		ax.add_label(0.28, 0.95, 'SO - '+thresh, color='black', relative=True, weight='bold')
		ax.ticks.show()
		ax.ticks.set_xspacing(nx)
		ax.ticks.set_yspacing(0.01)
		ax.ticks.set_color('black')
		ax.tick_labels.hide()
		ax.axis_labels.hide()
except IOError:
	print('No SO Data.')
	fig.text(0.55,0.16,'No SO Data',fontsize=14, fontweight='bold')

plt.subplots_adjust(wspace=0, hspace=0)
sourcename = sourcename.replace('.','')
plt.savefig(sourcename+'_moment_'+moment+'.pdf')
plt.show()


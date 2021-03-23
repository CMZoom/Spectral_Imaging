import os
import sys
import glob
import numpy as np
import time
import re

sourcename = sys.argv[3]
print sourcename
#Phase_Center = sys.argv[3] + " " + sys.argv[4] + " " + sys.argv[5]

if 'G359' in sourcename:
    L = sourcename[1:8]
    B = sourcename[8:14]
else:
    L = sourcename[1:6]
    B = sourcename[6:12]

Phase_Center = 'GALACTIC ' + L + ' ' +  B[1:]

imagesize=np.loadtxt('./ImageSize.txt',dtype='str')

#os.chdir('/reduction/czdata/line_cubes/')

for i in range(0,len(imagesize[:]),1):
        if sourcename in imagesize[i][0]:
                ImageSize = [int(imagesize[i][2]),int(imagesize[i][2])]

path = os.getcwd()

numbers = re.compile(r'(\d+)')
def numericalSort(value):
	parts= numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts
'''
for infile in sorted(glob.iglob('./'+sourcename+'_CASA/'+sourcename+'*.fits'), key=numericalSort):
	pathway = len('./'+sourcename+'_CASA/')
	t = os.path.splitext(infile)
	f = str(t[0])
	f = f[pathway:]
	importuvfits(fitsfile=infile,vis=path+'/'+sourcename+'_CASA/'+f+'.cube.ms')
#	cvel2(vis=path+'/'+sourcename+'_CASA/'+f+'.cube.ms', outputvis=path+'/'+sourcename+'_CASA/'+f+'.cube.lsrk.ms', mode='frequency', outframe='lsrk')
'''
print sourcename
os.chdir('./'+sourcename+'_CASA/')

##################
### Line Clean ###
##################

'''
line_usb_asic = []
line_lsb_asic = []
line_usb_swarm = []
line_lsb_swarm = []

print "Concatenating files..."

for infile in sorted(glob.iglob('*.cube.ms'), key=numericalSort):
        if 'usb_asic' in infile:
		line_usb_asic.append(infile)
        if 'lsb_asic' in infile:
                line_lsb_asic.append(infile)
        if 'usb_swarm' in infile:
                line_usb_swarm.append(infile)
        if 'lsb_swarm' in infile:
                line_lsb_swarm.append(infile)

print line_usb_asic, line_lsb_asic

concat(vis=line_usb_asic, concatvis=sourcename+'.line.usb_asic.concat.ms')
concat(vis=line_lsb_asic, concatvis=sourcename+'.line.lsb_asic.concat.ms')
concat(vis=line_usb_swarm, concatvis=sourcename+'.line.usb_swarm.concat.ms')
concat(vis=line_lsb_swarm, concatvis=sourcename+'.line.lsb_swarm.concat.ms')

#os.system('rm -rf *.cube.ms')

for filename in glob.iglob('*.concat.ms'):
#	os.system('rm -rf '+filename+'.ms.contsub')
	print filename
	if 'usb_asic' in filename:
		uvcontsub(vis=filename, fitspw='*:230.35~230.75GHz', excludechans=True, combine='spw')
	if 'lsb_asic' in filename:
		uvcontsub(vis=filename, fitspw='*:220.25~220.45GHz', excludechans=True, combine='spw')
	if 'lsb_swarm' in filename:
		uvcontsub(vis=filename, fitspw='*:219.05~219.09GHz', excludechans=True, combine='spw')
	if 'usb_swarm' in filename:
		uvcontsub(vis=filename, fitspw='*:234.69~234.71GHz', excludechans=True, combine='spw')
'''
print "Cleaning spectral cubes..."
#RestFrequency = '230.537970GHz'        # Rest frequency for data cubes.
CellSize = '0.5arcsec'                 # Specify cell size.
niter_line = 0		                    # No. iterations for line clean.
thresh_line = '120mJy'                 # Threshold for line clean.
phasecenter_txt = np.loadtxt('../phasecenter.txt', dtype=str)

for i in range(0,len(phasecenter_txt[:]),1):
        if sourcename in phasecenter_txt[i][0]:
                PhaseCenter = str(phasecenter_txt[i][1]+' '+phasecenter_txt[i][2]+' '+phasecenter_txt[i][3])

print PhaseCenter

if sourcename == 'G0.380+0.050':
	usb_vis = [sourcename+'.line.usb_asic.concat.ms.contsub', '../G0.376+0.040_CASA/G0.376+0.040.line.usb_asic.concat.ms.contsub']
	lsb_vis = [sourcename+'.line.lsb_asic.concat.ms.contsub', '../G0.376+0.040_CASA/G0.376+0.040.line.lsb_asic.concat.ms.contsub']
elif sourcename == 'G359.863-0.069' or sourcename == 'G359.889-0.093':
	usb_vis = ['../G359.863-0.069_CASA/G359.863-0.069.line.usb_asic.concat.ms.contsub', '../G359.889-0.093_CASA/G359.889-0.093.line.usb_asic.concat.ms.contsub']
	lsb_vis = ['../G359.863-0.069_CASA/G359.863-0.069.line.lsb_asic.concat.ms.contsub', '../G359.889-0.093_CASA/G359.889-0.093.line.lsb_asic.concat.ms.contsub']
else:
	usb_vis = sourcename+'.line.usb_asic.concat.ms.contsub'
	lsb_vis = sourcename+'.line.lsb_asic.concat.ms.contsub'

print lsb_vis
'''
print "Cleaning usb asic"

if not os.path.isdir(sourcename+'.asic.usb.multiscale.dirty.image'):
	tclean(vis=usb_vis,
	imagename = sourcename+'.asic.usb.multiscale.dirty',
	imsize = [100,100],
	cell = CellSize,
	phasecenter = PhaseCenter,
	specmode = 'cube',
	outframe = 'LSRK',
	gridder = 'mosaic',
	deconvolver = 'multiscale',
	scales=[0,3,9,27],
	weighting = 'briggs',
	start=50,
	nchan=100,
	width=1,
	robust = 0.5,
	niter = 0,
	gain = 0.1,
	threshold = '225mJy',
	interactive = False,
	mosweight=False,
	chanchunks=-1)

dirty_im_usb = imstat(sourcename+'.asic.usb.multiscale.dirty.image', axes=[0,1])

freq_usb = [230.538000, 231.060993]
line_usb = ['12CO', 'OCS']
#freq_usb = [231.060993]
#line_usb = ['OCS']

if ImageSize[0] > 100:
	for i in range(len(line_usb)):
		print freq_usb[i]-(205*0.0008125)
		tclean(vis=usb_vis,
		imagename = sourcename+'.asic.usb.multiscale.'+line_usb[i]+'.deep.test',
		imsize = ImageSize,
		cell = CellSize,
		phasecenter = PhaseCenter,
		specmode = 'cube',
		outframe = 'LSRK',
		gridder = 'mosaic',
		deconvolver = 'multiscale',
		scales = [0,3,9,27],
		weighting = 'briggs',
		start=str(freq_usb[i]-(184*0.0008125))+'GHz',
		width='0.8125MHz',
		nchan=369,
		restfreq=str(freq_usb[i])+'GHz',
		robust = 0.5,
		niter = 1000000,
		gain = 0.1,
		threshold = str(3*np.average(dirty_im_usb['rms']))+'Jy',
		interactive = False,
		mosweight=False,
		chanchunks=-1)
else:
	tclean(vis=usb_vis,
		imagename = sourcename+'.asic.usb.multiscale.deep.test',
		imsize = ImageSize,
		cell = CellSize,
		phasecenter = PhaseCenter,
		specmode = 'cube',
		outframe = 'LSRK',
		gridder = 'mosaic',
		deconvolver = 'multiscale',
		scales = [0,3,9,27],
		weighting = 'briggs',
		width='0.8125MHz',
		robust = 0.5,
		niter = 1000000,
		gain = 0.1,
		threshold = str(3*np.average(dirty_im_usb['rms']))+'Jy',
		interactive = False,
		mosweight=False,
		chanchunks=-1)
'''
print "Cleaning lsb asic"

if not os.path.isdir(sourcename+'.asic.lsb.multiscale.dirty.image'):
	tclean(vis=lsb_vis,
	imagename = sourcename+'.asic.lsb.multiscale.dirty',
	imsize = [100,100],
	cell = CellSize,
	phasecenter = PhaseCenter,
	specmode = 'cube',
	outframe = 'LSRK',
	gridder = 'mosaic',
	deconvolver = 'multiscale',
	scales=[0,3,9,27],
	weighting = 'briggs',
	start=50,
	nchan=100,
	width=1,
	robust = 0.5,
	niter = 0,
	gain = 0.1,
	threshold = '225mJy',
	interactive = False,
	mosweight=False,
	chanchunks=-1)

dirty_im_lsb = imstat(sourcename+'.asic.lsb.multiscale.dirty.image', axes=[0,1])

#line_lsb = ['H2CO_218.2', 'H2CO_218.5', 'H2CO_218.8', '13CO', 'C18O', 'SiO', 'OCS', 'SO']
#freq_lsb = [218.222192, 218.475632, 218.760066, 220.39868420, 219.56035410,  217.104980, 218.90335550, 219.94944200]
line_lsb = ['C18O']
freq_lsb = [219.56035410]

print ImageSize[0]

if ImageSize[0] > 100:
	for i in range(len(line_lsb)):
		print line_lsb[i], freq_lsb[i]
		tclean(vis=lsb_vis,
		imagename = sourcename+'.asic.lsb.multiscale.'+line_lsb[i]+'.deep.test.err',
		imsize = ImageSize,
		cell = CellSize,
		phasecenter = PhaseCenter,
		specmode = 'cube',
		outframe = 'LSRK',
		gridder = 'mosaic',
		deconvolver = 'multiscale',
		scales=[0,3,9,27],
		weighting = 'briggs',
		start=str(freq_lsb[i]-(184*0.0008125))+'GHz',
		width='0.8125MHz',
		nchan=369,
		restfreq=str(freq_lsb[i])+'GHz',
		robust = 0.5,
		niter = 1000000,
		gain = 0.1,
		threshold = str(3*np.average(dirty_im_lsb['rms']))+'Jy',
		interactive = False)
else:
	print 'Didnt work'
	tclean(vis=lsb_vis,
	imagename = sourcename+'.asic.lsb.multiscale.deep',
	imsize = ImageSize,
	cell = CellSize,
	phasecenter = PhaseCenter,
	specmode = 'cube',
	outframe = 'LSRK',
	gridder = 'mosaic',
	deconvolver = 'multiscale',
	scales=[0,3,9,27],
	weighting = 'briggs',
	start='216.9395GHz',
	width='0.8125MHz',
	robust = 0.5,
	niter = 1000000,
	gain = 0.1,
	threshold = str(3*np.average(dirty_im_lsb['rms']))+'Jy',
	interactive = False,
	mosweight=False,
	chanchunks=-1)
'''
tclean(vis=lsb_vis,
imagename = sourcename+'.asic.lsb.multiscale.deep.recheck',
imsize = ImageSize,
cell = CellSize,
phasecenter = PhaseCenter,
specmode = 'cube',
outframe = 'LSRK',
gridder = 'mosaic',
deconvolver = 'hogbom',
#scales=[0,3,9,27],
weighting = 'briggs',
start='216.9395GHz',
width='0.8125MHz',
robust = 0.5,
niter = 1000000,
gain = 0.1,
threshold = '225mJy',
interactive = False,
usemask='auto-multithresh',
sidelobethreshold=1.25,
noisethreshold=4.25,
minbeamfrac=0.1,
lownoisethreshold=1.50,
negativethreshold=0.0,
growiterations=75,
mosweight=False)
#chanchunks=-1)
'''
'''
print field
#partition(vis=sourcename+'.line.usb_asic.concat.ms.contsub', outputvis=sourcename+'.usb.input.ms', numsubms='auto', separationaxis='auto')
clean(vis=sourcename+'.line.usb_asic.concat.ms.contsub', 
	imagename = sourcename+'.asic.usb.pbcor.'+str(niter_line)+'.single_fielda.'+str(thresh_line),
	threshold = threshold,
	field = field,
        niter=150000,
        gain=0.1,
        psfmode='hogbom',
        imagermode=imagermode,
        interactive=False,
        imsize=ImageSize,
        cell="0.5arcsec",
#	resmooth = True,
        robust=0.5,
        weighting='briggs',
	mode = 'frequency',
#        specmode='cube',
	cyclespeedup=75,
        phasecenter = phasecenter)

clean(vis=sourcename+'.line.lsb_asic.concat.ms.contsub', 
	imagename = sourcename+'.asic.lsb.pbcor.'+str(niter_line)+'.single_fielda.'+str(thresh_line),
	threshold = threshold,
	field = field,
        niter=150000,
        gain=0.1,
        psfmode='hogbom',
        imagermode=imagermode,
        interactive=False,
        imsize=ImageSize,
        cell="0.5arcsec",
#	resmooth = True,
        robust=0.5,
        weighting='briggs',
	mode = 'frequency',
#        specmode='cube',
	cyclespeedup=75,
        phasecenter = phasecenter)

#	spw = '*:228.9~230.86GHz,230.9~232.87GHz',
	imsize = ImageSize,
	cell = CellSize,
	phasecenter = Phase_Center,
	specmode='cube',
#	start='228.9 GHz',
#	outframe = 'LSRK',
	restfreq = RestFrequency,
	pbcor = False,
	gridder = 'mosaic',
	deconvolver = 'hogbom',
	interactive = False,
	cyclefactor = 5,
	weighting = 'briggs',
	gain = 0.1,
	robust = 0.5,
	niter = niter_line,
	chanchunks = -1,
	threshold = thresh_line,
	restart=True,
	parallel=True)
'''
'''
print "Cleaning lsb asic"
partition(vis=sourcename+'.line.lsb_asic.concat.ms.contsub', outputvis=sourcename+'.lsb.input.ms', numsubms='auto', separationaxis='auto')
tclean(vis = sourcename+'.lsb.input.ms',
	imagename = sourcename+'.line.asic.lsb.pbcor.smoothed.'+str(niter_line)+'.'+str(thresh_line),
#	spw = '*:216.939~218.858GHz,218.895~220.86GHz',
	imsize = ImageSize,
	cell = CellSize,
	phasecenter = Phase_Center,
	specmode='cube',
#	start='216.9 GHz',
#	outframe = 'LSRK',
	restfreq = RestFrequency,
	pbcor = True,
	gridder = 'mosaic',
	deconvolver = 'hogbom',
	interactive = False,
	cyclefactor = 5,
	weighting = 'briggs',
	gain = 0.1,
	robust = 0.5,
	niter = niter_line,
	chanchunks = -1,
	threshold = thresh_line,
	restart=True,
	parallel=True)
print "Clean complete."
'''
'''
os.system('rm -rf ' + sourcename+'.line.swarm.usb.concat.dirty.0'+field +'*')
if len(line_usb_swarm) > 0:
#	concat(vis=line_usb_swarm, concatvis=sourcename+'.line.swarm.usb.concat.ms')
	print "Cleaning usb swarm"
	tclean(vis = sourcename+'.line.usb_swarm.concat.ms.contsub',
		imagename = sourcename+'.line.swarm.usb.concat.dirty.contsub.'+field,
	        field = field,
		imsize = ImageSize,
		cell = '0.5arcsec',
		phasecenter = Phase_Center,
		specmode='cube',
		start='230.75 GHz',
		width = '5 MHz',
		outframe = 'LSRK',
		gridder = 'standard',
		deconvolver = 'hogbom',
		interactive = False,
		cyclefactor = 5,
		weighting = 'briggs',
		gain = 0.1,
		robust = 0.5,
		niter = niter_line,
		threshold = thresh_line)
	print "Clean complete."

os.system('rm -rf ' + sourcename+'.line.swarm.lsb.concat.dirty.0'+field +'*')
if len(line_lsb_swarm) > 0:
#	concat(vis=line_lsb_swarm, concatvis=sourcename+'.line.swarm.lsb.concat.ms')
	print "Cleaning lsb swarm"
	tclean(vis = sourcename+'.line.lsb_swarm.concat.ms.contsub',
		imagename = sourcename+'.line.swarm.lsb.concat.dirty.contsub.'+field,
	        field = field,
		imsize = ImageSize,
		cell = '0.5arcsec',
		phasecenter = Phase_Center,
		specmode='cube',
		start='212.776 GHz',
		width = '5 MHz',
		outframe = 'LSRK',
		gridder = 'standard',
		deconvolver = 'hogbom',
		interactive = False,
		cyclefactor = 5,
		weighting = 'briggs',
		gain = 0.1,
		robust = 0.5,
		niter = niter_line,
		threshold = thresh_line)
	print "Clean complete."
'''

#os.system('rm -rf *.workdirectory')
#os.system('rm -rf *.input.ms.*')
'''
print "Converting images to Galactic Coordinates"
for filename in glob.iglob('*.image*'):
	print filename
	t = os.path.splitext(filename)
	f = str(t[0])
	imregrid(imagename=filename,template='GALACTIC',output=f+'.galactic')
#	exportfits(imagename=f+'.galactic',fitsimage=f+'.fits',dropdeg=True,velocity=True)
for filename in glob.iglob('*.residual*'):
	t = os.path.splitext(filename)
	f = str(t[0])
	imregrid(imagename=filename,template='GALACTIC',output=filename+'.galactic')
#	exportfits(imagename=filename+'.galactic',fitsimage=filename+'.fits',velocity=True)
'''

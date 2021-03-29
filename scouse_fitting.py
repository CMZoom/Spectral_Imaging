import sys
import numpy as np
from astropy.io import fits
import matplotlib
#matplotlib.use('Qt5Agg')
#matplotlib.rcParams['backend'] = 'Qt5Agg' 
matplotlib.rcParams['keymap.quit_all'] = 'Q'
matplotlib.rcParams['keymap.quit'] = 'q'
sys.path.append('scousepy-delta')
import matplotlib.pyplot as plt
import lmfit
import pyspeckit
from scousepy.scousefitter import ScouseFitter
import glob
import pickle
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.stats import mad_std

def writeDict(dict, filename, sep):
	with open(filename, "a") as f:
		for i in dict.keys():
			print(i, dict[i])            
			f.write(i + " " + sep.join([str(x) for x in dict[i]]) + "\n")

def openfile(filename):
	a = fits.open(filename)
	header = a[0].header
	data = a[0].data
	return header, data

mode = input('Select either Single or Full: ')

if mode.lower() == 'single':
	'''
	This allows the user to update the scouse fitting for a specific spectrum, as long as leaf number is known.

	TODO: add ability to select transition as well.
	'''
	leaf = int(input('Enter desired leaf number: '))
	trans = input('Enter transition: ')
	if trans == 'all':
		line = '*GHz_spectrum.fits'
	else:
		line = trans+'GHz_spectrum.fits'
	for filename in glob.glob(line):
		header, data = openfile(filename)

		with open(filename[:-4]+'test.dict.pkl', 'rb') as f:
			spec_out = pickle.load(f)

		leaves = data[:,0]
		try:
			l_id = np.where(leaves == leaf)[0][0]
			print(l_id)

			l = int(((np.shape(data)[1]-1)/2)+1)
			spectrum = data[l_id,1:l]
			spec_axis = data[l_id,l:]
			#print(type(spectrum))
			temp_data = []
			#kern = Gaussian1DKernel(10)
			#smoothed = convolve(spectrum,kern)
			#plt.figure(figsize=(8,8))
			#plt.plot(spec_axis, spectrum)
			#plt.plot(spec_axis,smoothed)
			#plt.show()
			#grad=np.gradient(smoothed)/(spec_axis[1]-spec_axis[0])
			#d2 = grad/(spec_axis[1]-spec_axis[0])
			#d3 = d2/(spec_axis[1]-spec_axis[0])
			#d4 = d3/(spec_axis[1]-spec_axis[0])
			#noise = mad_std(spectrum)
			#condition1 = np.array(spectrum>(5*noise), dtype='int')[1:]
			#condition2 = np.array(d2[1:]<0, dtype='int')
			#condition3 = np.array(d4[1:]>0, dtype='int')
			#condition4 = np.array(np.abs(np.diff(np.sign(d3)))!=0,dtype='int')
			#conditionmask = condition1*condition2*condition3*condition4
			#print(conditionmask)
			#print(np.size(np.where(conditionmask)))
			temp_data.append([spec_axis, spectrum])	
			#print(np.shape(temp_data))
			temp_dict = {}
			print(np.shape(temp_data))
			d = ScouseFitter(modelstore=temp_dict, spectra=l_id, method='individual', individual=np.array(temp_data),SNR=5, kernelsize=5, unit='jy beam-1')
			plt.show()
	
			#print(spec_out[l_id])
			temp_dict[0]['leaf_ID'] = leaf
			#print(temp_dict)
			spec_out[l_id] = temp_dict[0]
			print(spec_out[l_id])
	
			f = open(filename[:-4]+'test.dict.pkl', 'wb')
			pickle.dump(spec_out,f)
		except IndexError:
			print('No leaf found')
elif mode.lower() == 'full':
	'''
	This allows the user to do fit all spectra for a given transition.
	'''
	for filename in glob.glob('*_spectrum.fits'):
		spec_out = {}
		a = fits.open(filename)
		header = a[0].header
		data = a[0].data

		leaves = data[:,0]
		l = int(((np.shape(data)[1]-1)/2)+1)
		spectrum = data[:,1:l]
		spec_axis = data[:,l:]
		new_data = []

		for i in range(np.shape(spectrum)[0]):
			new_data.append([spec_axis[i],spectrum[i]])

		print(new_data)
		d = ScouseFitter(modelstore=spec_out, spectra=leaves, method='individual', individual=np.array(new_data), SNR=5, kernelsize=5)
		plt.show()

		for key in spec_out.keys():
			spec_out[key]['leaf_ID'] = int(leaves[key])

		print(spec_out)
		f = open(filename[:-4]+'test.dict.pkl','wb')
		pickle.dump(spec_out, f)
else:
	print('Select valid option (Single or Full)')
line_dict = {0: '12CO.230.5GHz', 1: '13CO.220.4GHz', 2: 'C18O.219.6GHz', 3: 'H2CO.218.2GHz', 4: 'H2CO.218.5GHz', 5: 'H2CO.218.8GHz', 6: 'OCS.218.9GHz', 7: 'OCS.231.1GHz', 8: 'SiO.217.1GHz', 9: 'SO.219.9GHz'}
'''
for filename in glob.glob('C18O*GHz_spectrum.fits'):
	spec_out = {}
	a = fits.open(filename)
	header = a[0].header
	data = a[0].data

	leaves = data[:,0]
	l = int(((np.shape(data)[1]-1)/2)+1)
	spectrum = data[:,1:l]
	spec_axis = data[:,l:]
	new_data = []

	for i in range(np.shape(spectrum)[0]):
		new_data.append([spec_axis[i],spectrum[i]])

	print(new_data)
	d = ScouseFitter(modelstore=spec_out,spectra=leaves,method='individual',individual=np.array(new_data),SNR=5,kernelsize=5)
	plt.show()
	
	for key in spec_out.keys():
		spec_out[key]['leaf_ID'] = int(leaves[key])
	print(spec_out)
	f = open(filename[:-4]+'test.dict.pkl','wb')
	pickle.dump(spec_out, f)
	#f = open(filename[:-4]+'test.dict.txt', 'w+')
	#f.write(json.dumps(spec_out, default= lambda x: None))
'''
'''
for row in data:
	leaf_dict = {}
	leaf = row[0]
	spectrum = row[1:149]
	spec_axis = row[149:]
	error = mad_std(spectrum)*np.ones_like(spectrum)

	d = ScouseFitter(modelstore=leaf_dict,method='individual',individual=np.array([[spec_axis,spectrum]]))
	print(leaf_dict)
	plt.show()
#d = Decomposer(spectral_axis=np.linspace(-50, 50, 100),spectrum=sp,rms=5)
#d.method = 'manual'
#d.fit_spectrum_manually()
#Decomposer.fit_spectrum_manual(spectrum)
'''




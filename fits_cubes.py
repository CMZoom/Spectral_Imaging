from astropy.io import fits
import numpy as np
import sys
import astropy.units as u
sourcename = sys.argv[1]

line_dict = {"H2CO": [["3(0,3)-2(0,2)", 218.222192*u.GHz],["3(2,2)-2(2,1)", 218.475632*u.GHz],["3(2,1)-2(2,0)", 218.760066*u.GHz]], 
             "13CO": [["J=2-1", 220.39868420*u.GHz]], "C18O": [["J=2-1", 219.56035410*u.GHz]], "SiO": [["J=5-4", 217.104980*u.GHz]], 
             "12CO": [["J=2-1", 230.538000*u.GHz]], "OCS": [["18-17", 218.90335550*u.GHz],["19-18", 231.06099340*u.GHz]],
             "SO": [["6-5", 219.94944200*u.GHz]]}

def fits_cube(data, original_header, line, freq):

	subcube_freq = np.arange(freq-0.15, freq+0.15+original_header['cdelt3']/1e9, original_header['cdelt3']/1e9)
	subcube_freq = subcube_freq[1:]
	print len(subcube_freq)
	subcube_vel = ((1. - subcube_freq[:-1]/freq)*3e5)

	line_freq = np.linspace(original_header['crval3'], original_header['crval3']+(original_header['naxis3']*original_header['cdelt3']), original_header['naxis3'])/1e9

	subcube = []

	for k in np.where((line_freq >= subcube_freq[0]) & (line_freq <= subcube_freq[-1]))[0]:
		subcube.append(data[k])

	subcube = np.array(subcube)

	new_header = original_header.copy()
	new_header['NAXIS3'] = len(np.where((line_freq >= subcube_freq[0]) & (line_freq <= subcube_freq[-1]))[0])
	new_header['CTYPE3'] = 'VRAD'
	new_header['CRVAL3'] = subcube_vel[0]
	new_header['CDELT3'] = subcube_vel[1] - subcube_vel[0]
	new_header['CUNIT3'] = 'km/s'
	new_header['RESTFRQ'] = freq*1e9

	fits.writeto(sourcename+'_CASA/'+sourcename+'.'+line+'.'+str(round(freq,1))+'GHz.fits', subcube, header=new_header, overwrite=True)

usb_line_cube = fits.open(sourcename+'_CASA/'+sourcename+'.asic.usb.multiscale.deep.convolved.galactic.pbcor.fits')
usb_line_cube_data = usb_line_cube[0].data
usb_line_cube_header = usb_line_cube[0].header
print usb_line_cube_header['CUNIT3']
lsb_line_cube = fits.open(sourcename+'_CASA/'+sourcename+'.asic.lsb.multiscale.deep.convolved.galactic.pbcor.fits')
lsb_line_cube_data = lsb_line_cube[0].data
lsb_line_cube_header = lsb_line_cube[0].header
print lsb_line_cube_header['CUNIT3']
for line in line_dict.keys():
	rest_freq = line_dict[line]
	for i in range(len(rest_freq)):
		print rest_freq[i][1].value
		if rest_freq[i][1].value > 216 and rest_freq[i][1].value < 221:
			print 'here'
			fits_cube(lsb_line_cube_data, lsb_line_cube_header, line, rest_freq[i][1].value)
		else:
			fits_cube(usb_line_cube_data, usb_line_cube_header, line, rest_freq[i][1].value)

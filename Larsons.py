#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 13:21:54 2020

@author: danielcallanan
"""

import pandas as pd
import os
import astropy.constants as c
import astropy.units as u
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

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

os.chdir('../Documents/scouse/')
vals = pd.read_csv('velocity_vals.csv')

catalog = fits.open('catalog_acc-2.fits')
cat_cols = catalog[1].columns
cat_data = catalog[1].data

new_catalog = fits.open('megacatalog_acc.fits')
new_cat_cols = new_catalog[1].columns
new_cat_data = new_catalog[1].data

def R(sig):
    return np.power(sig/1.1,1/0.5)

def M(sig):
    return np.power(sig/0.42,1/0.20)

def vel_disp(R_val, M_val):
    R_val = R_val*u.pc
    R_val = R_val.to(u.m)
    M_val = M_val*u.M_sun
    M_val = M_val.to(u.kg)
    G = c.G
    alp = 1.0
    return np.sqrt(alp*G*M_val/(5*R_val)).to(u.km/u.s).value

def alpha(sig, R_val, M_val):
    sig = sig*1000*(u.m/u.s)
    R_val = R_val*u.pc
    R_val = R_val.to(u.m)
    M_val = M_val*u.M_sun
    M_val = M_val.to(u.kg)
    G = c.G
    return (5*(sig**2)*R_val)/(G*M_val)

def alpha_unc(sig, sig_unc, R_val, M_val, M_unc):
    G = c.G
    sig = sig*1000*(u.m/u.s)
    sig_unc = sig_unc*100*(u.m/u.s)
    R_val = R_val*u.pc
    R_val = R_val.to(u.m)
    M_val = M_val*u.M_sun
    M_val = M_val.to(u.kg)
    M_unc = M_unc*u.M_sun
    M_unc = M_unc.to(u.kg)
    sa = ((10*sig*R_val)/(G*M_val))*sig_unc
    sb = (-(5*sig**2*R_val)/(G*M_val**2))*M_unc
    unc = np.sqrt(sa**2+sb**2)
    return unc

vel = []
fwhms = []
fwhm_ints = []
sigma = []

fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(16,8))
fig2,ax2 = plt.subplots(nrows=1,ncols=1,figsize=(8,8))
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim([0.01,100])
ax2.set_ylim([0.01,100])
for id in np.unique(vals['id']):
    fwhm = vals[vals['id'] == id]['h_fwhm'].to_numpy()
    fwhms.append(fwhm[0])
    fwhm_unc = vals[vals['id'] == id]['h_fwhm_unc'].to_numpy()
    fwhm_int = vals[vals['id'] == id]['h_fwhm_int'].to_numpy()
    fwhm_ints.append(fwhm_int[0])
    fwhm_int_unc = vals[vals['id'] == id]['h_fwhm_int_unc'].to_numpy()
    mass = new_cat_data[new_cat_data['leaf_ID'] == id]['mass']
    mass_unc = new_cat_data[new_cat_data['leaf_ID'] == id]['mass_unc']
    rad = new_cat_data[new_cat_data['leaf_ID'] == id]['r_eff_pc']
    vel.append(vel_disp(rad,mass)[0])
    sigma = new_cat_data[new_cat_data['leaf_ID'] == id]['Sigma']
    ax2.scatter(rad,fwhm,marker='x',color='black')
    ax2.scatter(rad,fwhm_int, marker='+',color='red')
    ax2.scatter(rad,vel_disp(rad,mass)[0],marker='o',color='blue')
    #if fwhm > 5:
        #print(id,fwhm)
    a = alpha(fwhm,rad,mass)
    a_unc = alpha_unc(fwhm,fwhm_unc,rad,mass,mass_unc)
    if np.isnan(fwhm_int) and ~np.isnan(fwhm):
        ax[0].errorbar(mass,a,a_unc,mass_unc,color='black',label='obs. velocity dispersion ($<$ vel. res.)')
    else:
        ax[0].errorbar(mass,a,a_unc,mass_unc,color='red',label='obs. velocity dispersion')
    ax[1].errorbar(fwhm,a,a_unc,fwhm_unc,color='red')
    a_int = alpha(fwhm_int,rad,mass)
    a_int_unc = alpha_unc(fwhm_int,fwhm_int_unc,rad,mass,mass_unc)
    ax[0].errorbar(mass,a_int,a_int_unc,mass_unc,color='blue',alpha=0.25,label='corrected velocity dispersion')
    ax[1].errorbar(fwhm_int,a_int,a_int_unc,fwhm_int_unc,color='blue',alpha=0.25)
    if a <= 1:
        print(id)
    #if np.isnan(fwhm_int) and ~np.isnan(fwhm):
            #ax[0].scatter(mass,a,marker='x',color='black')
 
handles, labels = ax[0].get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax[0].legend(by_label.values(), by_label.keys(),fontsize=12)
#ax[0].legend(['obs. velocity dispersion','int. velocity dispersion','obs. velocity dispersion (< vel. res.)'])
ax[1].legend(['obs. velocity dispersion','corrected velocity dispersion'],fontsize=12,loc='lower right')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
#ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[0].set_xlim(5)
ax[0].axhspan(0,1,color='gray',alpha=0.25)
ax[1].axhspan(0,1,color='gray',alpha=0.25)
ax[0].set_xlabel('Mass')
ax[1].set_xlabel('Velocity Dispersion')
ax[0].set_ylabel('Virial Parameter')
ax[1].axvline(1.057309424046252,color='black',linestyle='--')
ax[1].text(1.2,40,'Channel Width',rotation=90,weight='bold')
fig.savefig('quick_virial_unfilter.pdf')

plt.figure(figsize=(8,8))
plt.hist(vel,bins=10,label=r'if $\alpha$ = 1')
plt.hist(fwhm_ints,bins=10,alpha=0.5,label='Measured')
plt.axvline(1.057309424046252,color='black',linestyle='--')
plt.xlabel(r'Velocity Dispersion')
plt.text(1.2,15,'Channel Width',rotation=90)
plt.ylabel('Number')
plt.legend(loc='upper right')
#plt.savefig('Dispersion_hist.pdf')

plt.figure(figsize=(8,8))
plt.scatter(fwhms,np.array(fwhms)-np.array(fwhm_ints),marker='x')
plt.xlabel('Observed Velocity Dispersion')
plt.ylabel('Change in Velocity Dispersion')
#plt.savefig('Disp_comparison.pdf')

#plt.legend(['obs. velocity dispersion','int. velocity dispersion'])
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0.1)
#plt.axhline(1.0,color='gray',linestyle='--')
#plt.xlabel('Mass')
#plt.ylabel('Virial Parameter')
#plt.savefig('quick_virial_unfilter.pdf')

plt.figure(figsize=(8,8))
plt.scatter(cat_data['mass_bgsub'],cat_data['r_eff_pc'])
plt.xscale('log')
    

fwhm = np.array(vals['h_fwhm'])
R_vals = R(np.array(vals['h_fwhm']))
M_vals = M(np.array(vals['h_fwhm']))
    
a = alpha(fwhm,R_vals,M_vals)

#print(a)

    
    

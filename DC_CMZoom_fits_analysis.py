#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 17:00:01 2020

Script to analyse output of automated spectral line fitting of CMZooom 
datacubes

@author: stevenlongmore
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from astropy.io import ascii
import astropy.constants as c

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

line_dict = {"13CO.220.4": ["J=2-1", 220.39868420], "C18O.219.6": ["J=2-1", 219.56035410],
             "H2CO.218.2": ["3(0,3)-2(0,2)",218.222192], "H2CO.218.5": ["3(2,2)-2(2,1)",218.475632],
             "H2CO.218.8": ["3(2,1)-2(2,0)",218.760066], "SiO.217.1": ["J=5-4", 217.104980],
             "12CO.230.5": ["J=2-1", 230.538000], "OCS.218.9": ["18-17", 218.90335550],
             "OCS.231.1":["19-18", 231.06099340],"SO.219.9": ["6-5", 219.94944200]}
# -------------------------------
# read in fitting output

data = pd.read_csv('Leaf_ID_table_test.csv') 

# outputting sorted values for DC
#data.sort_values(['Region','Core_ID','Transition','Component']).to_csv('../Documents/scouse/leaf_ID_table_test_sorted.csv')
#sort = data.sort_values(['Region','Core_ID','Transition','Component', 'Detection','Number_of_comps','Amp','Amp_err','Vlsr','Vlsr_err','fwhm','fwhm_err','RMS'])
# Creating new column for unique leaf ID
data['Leaf_ID'] = data.apply(lambda row: row.Region+row.Core_ID, axis = 1)
data['Frequency'] = data.apply(lambda row: float(line_dict[row.Transition[:-3]][1])*(1.0-((row.Vlsr*1000)/c.c.value)), axis = 1)

ID = np.array(data[(data['Leaf_ID'] == 'G359.648-0.133a') & (data['Transition'] == 'H2CO.218.8GHz')].index)
data.drop(ID[1],inplace=True)
ID = np.array(data[(data['Leaf_ID'] == 'G359.648-0.133a') & (data['Transition'] == 'OCS.218.9GHz')].index)
data.drop(ID[1],inplace=True)
#Prepare and save data to latex table
#temp_trans =  np.array(sort['Transition'])
#temp_amp = np.array(sort['Amp'])
#temp_amp_err = np.array(sort['Amp_err'])
#temp_Vlsr = np.array(sort['Vlsr'])
#temp_Vlsr_err = np.array(sort['Vlsr_err'])
#temp_fwhm = np.array(sort['fwhm'])
#temp_fwhm_err = np.array(sort['fwhm'])

#new_trans = [str(temp_trans[i].split('.')[0]) + " (" + str(temp_trans[i].split('.')[1]) + "." + str(temp_trans[i].split('.')[2]) + ")" for i in range(len(sort['Transition']))]
#amp = ['$'+str(temp_amp[i])+'\pm'+str(temp_amp_err[i])+'$' for i in range(len(sort['Amp']))]
#vel = ['$'+str(temp_Vlsr[i])+'\pm'+str(temp_Vlsr_err[i])+'$' for i in range(len(sort['Vlsr']))]
#wid = ['$'+str(temp_fwhm[i])+'\pm'+str(temp_fwhm_err[i])+'$' for i in range(len(sort['fwhm']))]
#d = [sort['Leaf_ID'],new_trans,sort['Detection'],sort['Number_of_comps'],sort['Component'],amp,vel,wid,sort['RMS']]
#table_cols = ['Leaf','Transition','Detection','Number of Comps','Component','Amplitude','Velocity','Width','RMS']
#ascii.write(d, '../Documents/scouse/leaf_ID_table_test_sorted.tex', names=table_cols, format='latex', overwrite=True)

# -------------------------------
# Removing instrumental response from the FWHM by subtracting channel width
# in quadrature

# Velocity channel width --> check this is the same for all files!!
CRPIX3 = 1.057309424046252
data['fwhm_int']=data.apply(lambda row: row.fwhm**2 - CRPIX3**2, axis = 1)**0.5

# data['fwhm_int'].hist(log=True)

# -------------------------------
# quality controls use to exclude poor fits

vlsr_err_max = 1.5
fwhm_err_max = 1.5
amp_err_max  = 0.5
fwhm_max     = 15.0
fwhm_min     = 0.5

data_qc = data[
    (data['Detection'] == "Y") &
    (data['Vlsr_err']   < vlsr_err_max) &
    (data['fwhm_err']   < fwhm_err_max) &
    (data['Amp_err']    < amp_err_max) &
    (data['fwhm']       < fwhm_max) &
    (data['fwhm']       > fwhm_min) &
    (data['Transition'] != '12CO.230.5GHz') &
    (data['Transition'] != '13CO.220.4GHz')
     ]

# -------------------------------
#Identifying misidentified lines

walker_id, walker_vel = np.loadtxt('Walker_vels.txt', dtype=str, unpack=True)

data_ID = pd.DataFrame(columns=data.columns)

for line in np.unique(data_qc['Transition']):
    subcube = data_qc[data_qc['Transition'] == line]
    subcube['original_cube'] = np.nan
    plt.figure(figsize=(8,8))
    plt.suptitle(line)
    diffs = []
    for region in np.unique(subcube['Region']):
        vel = walker_vel[walker_id == region]
        if len(vel) > 0:
            vel = int(vel)
        elif len(vel) == 0 and region == 'G359.137+0.031':
            vel = 0
        else:
            print(region)
            continue
        core_vels = np.array(subcube[subcube['Region'] == region]['Vlsr'])
        [diffs.append(x) for x in core_vels-vel]
    plt.hist(diffs)
    for index in subcube.index:
        reg = subcube.loc[index,'Region']
        vel = walker_vel[walker_id == reg]
        if len(vel) > 0:
            vel = int(vel)
        else:
            continue
        if np.abs(subcube.loc[index,'Vlsr']-vel) > 25:
            subcube.loc[index,'original_cube'] = subcube.loc[index,'Transition']
            subcube.loc[index,'Transition'] = 'unknown'
    data_ID = data_ID.append(subcube)
    
# Correcting some peaks by hand
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058l') & (data_ID['Transition'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == 56.86)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058n') & (data_ID['Transition'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 68.96)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058n') & (data_ID['Transition'] == 'H2CO.218.8GHz') & (data_ID['Vlsr'] == 66.81)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058n') & (data_ID['Transition'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == 66.08)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058p') & (data_ID['Transition'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == 72.97)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.068-0.075c') & (data_ID['Transition'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == 63.53)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.068-0.075g') & (data_ID['Transition'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == 29.92)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['Transition'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == 29.79)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['Transition'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == 39.18)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'CH3OH'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['Transition'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 46.14)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['Transition'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == 40.69)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.106-0.082e') & (data_ID['Transition'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == 32.68)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'
ID = data_ID[(data_ID['Leaf_ID'] == 'G359.889-0.093a') & (data_ID['Transition'] == 'SiO.217.1GHz') & (data_ID['Vlsr'] == -3.86)].index
data_ID.loc[ID,'original_cube'] = data_ID.loc[ID,'Transition']
data_ID.loc[ID,'Transition'] = 'unknown'

ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058p') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 21.54)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.001-0.058p') & (data_ID['original_cube'] == 'H2CO.218.8GHz') & (data_ID['Vlsr'] == 22.8)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.8GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.068-0.075a') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 20.91)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.068-0.075a') & (data_ID['original_cube'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == 20.11)].index
data_ID.loc[ID,'Transition'] = 'C18O.219.6GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.068-0.075a') & (data_ID['original_cube'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == 21.37)].index
data_ID.loc[ID,'Transition'] = 'SO.219.9GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == -10.72)].index
data_ID.loc[ID,'Transition'] = 'C18O.219.6GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == -7.22)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == -9.17)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.5GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'H2CO.218.8GHz') & (data_ID['Vlsr'] == -9.5)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.8GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'OCS.218.9GHz') & (data_ID['Vlsr'] == -10.01)].index
data_ID.loc[ID,'Transition'] = 'OCS.218.9GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == 187.02)].index
data_ID.loc[ID,'Transition'] = 'CH3OH'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'OCS.218.9GHz') & (data_ID['Vlsr'] == 187.02)].index
data_ID.loc[ID,'Transition'] = 'H2CO'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'SiO.217.1GHz') & (data_ID['Vlsr'] == -8.02)].index
data_ID.loc[ID,'Transition'] = 'SiO.217.1GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035d') & (data_ID['original_cube'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == -9.07)].index
data_ID.loc[ID,'Transition'] = 'SO.219.9GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == -5.90)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == -8.50)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.5GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'H2CO.218.5GHz') & (data_ID['Vlsr'] == 40.69)].index
data_ID.loc[ID,'Transition'] = 'CH3OH'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'H2CO.218.8GHz') & (data_ID['Vlsr'] == -8.19)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.8GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'OCS.218.9GHz') & (data_ID['Vlsr'] == 188.32)].index
data_ID.loc[ID,'Transition'] = 'H2CO'
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'SiO.217.1GHz') & (data_ID['Vlsr'] == -7.68)].index
data_ID.loc[ID,'Transition'] = 'SiO.217.1GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.070-0.035e') & (data_ID['original_cube'] == 'SO.219.9GHz') & (data_ID['Vlsr'] == -7.74)].index
data_ID.loc[ID,'Transition'] = 'SO.219.9GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.212-0.001c') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 37.42)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.212-0.001d') & (data_ID['original_cube'] == 'C18O.219.6GHz') & (data_ID['Vlsr'] == 38.65)].index
data_ID.loc[ID,'Transition'] = 'C18O.219.6GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.326-0.085b') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 10.18)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G0.340+0.055a') & (data_ID['original_cube'] == 'H2CO.218.2GHz') & (data_ID['Vlsr'] == 23.3)].index
data_ID.loc[ID,'Transition'] = 'H2CO.218.2GHz'
data_ID.loc[ID,'original_cube'] = np.nan
ID = data_ID[(data_ID['Leaf_ID'] == 'G359.137+0.031a') & (data_ID['original_cube'] == 'OCS.218.9GHz') & (data_ID['Vlsr'] == 194.89)].index
data_ID.loc[ID,'Transition'] = 'H2CO'
#data_ID.drop(1529,inplace=True)
#data_ID.drop(1530,inplace=True)
#data_ID.loc[176,'original_cube'] = data_ID.loc[176,'Transition']
#data_ID.loc[176,'Transition'] = 'unknown'
#data_ID.loc[178,'original_cube'] = data_ID.loc[178,'Transition']
#data_ID.loc[178,'Transition'] = 'unknown'
#data_ID.loc[179,'original_cube'] = data_ID.loc[179,'Transition']
#data_ID.loc[179,'Transition'] = 'unknown'
#data_ID.loc[183,'original_cube'] = data_ID.loc[183,'Transition']
#data_ID.loc[183,'Transition'] = 'unknown'
#data_ID.loc[182,'original_cube'] = data_ID.loc[182,'Transition']
#data_ID.loc[182,'Transition'] = 'unknown'
#data_ID.loc[198,'original_cube'] = data_ID.loc[198,'Transition']
#data_ID.loc[198,'Transition'] = 'unknown'
#data_ID.loc[185,'original_cube'] = data_ID.loc[185,'Transition']
#data_ID.loc[185,'Transition'] = 'unknown'
#data_ID.loc[184,'original_cube'] = data_ID.loc[184,'Transition']
#data_ID.loc[184,'Transition'] = 'unknown'
#data_ID.loc[461,'original_cube'] = data_ID.loc[461,'Transition']
#data_ID.loc[461,'Transition'] = 'unknown'
#data_ID.loc[474,'original_cube'] = data_ID.loc[474,'Transition']
#data_ID.loc[474,'Transition'] = 'unknown'
#data_ID.loc[487,'original_cube'] = data_ID.loc[487,'Transition']
#data_ID.loc[487,'Transition'] = 'unknown'
#data_ID.loc[466,'original_cube'] = data_ID.loc[466,'Transition']
#data_ID.loc[466,'Transition'] = 'unknown'
#data_ID.loc[203,'original_cube'] = data_ID.loc[203,'Transition']
#data_ID.loc[203,'Transition'] = 'unknown'
#data_ID.loc[490,'original_cube'] = data_ID.loc[490,'Transition']
#data_ID.loc[490,'Transition'] = 'unknown'
#data_ID.loc[479,'original_cube'] = data_ID.loc[479,'Transition']
#data_ID.loc[479,'Transition'] = 'unknown'
#data_ID.loc[472,'original_cube'] = data_ID.loc[472,'Transition']
#data_ID.loc[472,'Transition'] = 'unknown'

# Separating unknowns and ID'ing them
unknowns = data_ID[data_ID['Transition'] == 'unknown']
IDs = data_ID[data_ID['Transition'] == 'unknown'].index

for id in np.unique(unknowns['Region']):
    vel = int(walker_vel[walker_id == id])*1000
    index = unknowns[unknowns['Region'] == id].index
    for ind in index:
        unknowns.loc[ind,'rest_freq']=unknowns.loc[ind,'Frequency']/(1.0 - (vel/float(c.c.value)))
  
from astroquery.splatalogue import Splatalogue
from astroquery.nist import Nist
import astropy.units as u

S = Splatalogue(energy_max=50,energy_type='eu_k',energy_levels=['el4'],line_strengths=['ls4'],only_NRAO_recommended=True,noHFS=True)
def trimmed_query(*args,**kwargs):
    columns = ('Species','Chemical Name','Resolved QNs','Log<sub>10</sub> (A<sub>ij</sub>)',
               'E_U (K)')
    table = S.query_lines(*args, **kwargs)[columns]
    table.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
    table.rename_column('E_U (K)','EU_K')
    table.rename_column('Resolved QNs','QNs')
    table.sort('EU_K')
    return table

for index,row in unknowns.iterrows():
    try:
        print(unknowns.loc[index,'Leaf_ID'],unknowns.loc[index,'Vlsr'])
        freq = unknowns.loc[index,'rest_freq']
        print(freq)
        columns = ('Species','Chemical Name','Freq-GHz(rest frame,redshifted)','Meas Freq-GHz(rest frame,redshifted)',
                   'Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)')
        lines = Splatalogue.query_lines((freq-0.02)*u.GHz, (freq+0.02)*u.GHz,top20='top20',
                                        energy_max=100,energy_type='eu_k',only_NRAO_recommended=True,noHFS=True,
                                        line_lists=['JPL','CDMS'])[columns]
        lines.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
        lines.rename_column('E_U (K)','EU_K')
        lines.rename_column('Freq-GHz(rest frame,redshifted)','Freq-GHz')
        lines.rename_column('Meas Freq-GHz(rest frame,redshifted)','Meas Freq-GHz')
        lines.sort('log10(Aij)')
        lines.reverse()
        if len(lines) > 0:
            print(lines)
            number = int(input('Enter preferred row: '))
            if number > 0:
                unknowns.loc[index,'Transition'] = lines[number-1]['Species']+'.'+str(lines[number-1]['Freq-GHz'])
            else:
                print('No line selected.')
    except ConnectionError:
        continue
    '''
    else:
        lines = Splatalogue.query_lines((freq-0.02)*u.GHz, (freq+0.02)*u.GHz,top20='top20',
                                    energy_max=50,energy_type='eu_k',energy_levels=['el4'],
                                    line_strengths=['ls4'],only_NRAO_recommended=True,noHFS=True,
                                    line_lists=['JPL','CDMS'])[columns]
        lines.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
        lines.rename_column('E_U (K)','EU_K')
        lines.sort('EU_K')
        print(lines)
    '''
    
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.001-0.058p') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == 21.54)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.001-0.058p') & (unknowns['Transition'] == 'H2CO.218.8GHz') & (unknowns['Vlsr'] == 21.54)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.068-0.075a') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == 20.91)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.068-0.075a') & (unknowns['Transition'] == 'C18O.219.6GHz') & (unknowns['Vlsr'] == 20.11)].index
unknowns.loc[ID,'Transition'] = 'C18O.219.6GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.068-0.075a') & (unknowns['Transition'] == 'SO.219.9GHz') & (unknowns['Vlsr'] == 21.37)].index
unknowns.loc[ID,'Transition'] = 'SO.219.9GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'C18O.219.6GHz') & (unknowns['Vlsr'] == -10.72)].index
unknowns.loc[ID,'Transition'] = 'C18O.219.6GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == -7.22)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'H2CO.218.5GHz') & (unknowns['Vlsr'] == -9.17)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.5GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'H2CO.218.8GHz') & (unknowns['Vlsr'] == -9.5)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'OCS.218.9GHz') & (unknowns['Vlsr'] == -10.01)].index
unknowns.loc[ID,'Transition'] = 'OCS.218.9GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'OCS.218.9GHz') & (unknowns['Vlsr'] == 187.02)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'SiO.217.1GHz') & (unknowns['Vlsr'] == -8.02)].index
unknowns.loc[ID,'Transition'] = 'SiO.217.1GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035d') & (unknowns['Transition'] == 'SO.219.9GHz') & (unknowns['Vlsr'] == -9.07)].index
unknowns.loc[ID,'Transition'] = 'SO.219.9GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == -5.90)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'H2CO.218.5GHz') & (unknowns['Vlsr'] == -8.50)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.5GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'H2CO.218.8GHz') & (unknowns['Vlsr'] == -8.19)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'OCS.218.9GHz') & (unknowns['Vlsr'] == 188.32)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'SiO.217.1GHz') & (unknowns['Vlsr'] == -7.68)].index
unknowns.loc[ID,'Transition'] = 'SiO.217.1GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.070-0.035e') & (unknowns['Transition'] == 'SO.219.9GHz') & (unknowns['Vlsr'] == -7.74)].index
unknowns.loc[ID,'Transition'] = 'SO.219.9GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.212-0.001c') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == 37.42)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.212-0.001d') & (unknowns['Transition'] == 'C18O.219.6GHz') & (unknowns['Vlsr'] == 38.65)].index
unknowns.loc[ID,'Transition'] = 'C18O.219.6GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.326-0.085b') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == 10.18)].index
unknowns.loc[ID,'Transition'] = 'C18O.219.6GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G0.340+0.055a') & (unknowns['Transition'] == 'H2CO.218.2GHz') & (unknowns['Vlsr'] == 23.3)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.2GHz'
unknowns.loc[ID,'original_cube'] = np.nan
unknowns.loc[ID,'rest_freq'] = np.nan
ID = unknowns[(unknowns['Leaf_ID'] == 'G359.137+0.031a') & (unknowns['Transition'] == 'OCS.218.9GHz') & (unknowns['Vlsr'] == 194.89)].index
unknowns.loc[ID,'Transition'] = 'H2CO.218.8GHz'

data_ID.drop(IDs, inplace=True)
data_ID = data_ID.append(unknowns)
#data_ID.to_csv('Documents/scouse/Leaf_spectra_ID.csv')


# -------------------------------
# Inspecting output

# data['Leaf_ID'].describe()

# print the column names

# data.columns
# 'Leaf_ID', 'Transition', 'Number_of_comps', 'Amp', 'Amp_err', 
# 'Vlsr','Vlsr_err', 'fwhm', 'fwhm_err', 'RMS'

# data['Number_of_comps'].hist()
# data['RMS'].hist(log=True)

# data['Transition'].unique()
# 'C18O.219.6GHz', '12CO.230.5GHz', '13CO.220.4GHz', 'H2CO.218.2GHz',
# 'H2CO.218.5GHz', 'H2CO.218.8GHz', 'OCS.218.9GHz', 'OCS.231.1GHz',
# 'SO.219.9GHz', 'SiO.217.1GHz'

# data['Leaf_ID']

# fwhm
# fwhm_lt10 = data[data['fwhm']<10.]
# fwhm_lt10['fwhm'].hist()
# data[data['fwhm']<10.]['fwhm'].hist()

# fwhm err
# data['fwhm_err'].hist(log=True)
# data[data['fwhm_err']<10.]['fwhm_err'].hist()

# ------------------------------
# Compare kinematic properties for all detections for every leaf

# getting unique leaf names
unique_all, counts_all       = np.unique(data['Leaf_ID'].values, return_counts=True)
unique_all_qc, counts_all_qc = np.unique(data_ID['Leaf_ID'].values, return_counts=True)

# creating variable and column names to hold kinematic info for each leaf
leaf_kin = pd.DataFrame(np.zeros([unique_all.size,10]))
leaf_kin.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max']

leaf_kin_qc = pd.DataFrame(np.zeros([unique_all_qc.size,10]))
leaf_kin_qc.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max']

# filling the leaf_kin variable with appropriate info for non qc
i=0
for id in unique_all:
    leaf_kin.iloc[i,0] = id
    leaf_kin.iloc[i,1] = data[data['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin.iloc[i,2] = data[data['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin.iloc[i,3] = data[data['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin.iloc[i,4] = data[data['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin.iloc[i,5] = data[data['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin.iloc[i,6] = data[data['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin.iloc[i,7] = data[data['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin.iloc[i,8] = data[data['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin.iloc[i,9] = data[data['Leaf_ID'] == id]['Vlsr'].max()
    i=i+1


#and for qc
i=0
for id in unique_all_qc:
    leaf_kin_qc.iloc[i,0] = id
    leaf_kin_qc.iloc[i,1] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin_qc.iloc[i,2] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin_qc.iloc[i,3] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin_qc.iloc[i,4] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin_qc.iloc[i,5] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin_qc.iloc[i,6] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin_qc.iloc[i,7] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin_qc.iloc[i,8] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin_qc.iloc[i,9] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].max()
    i=i+1

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['fwhm_std'].hist(ax=ax[0])
leaf_kin_qc['fwhm_std'].hist(ax=ax[1])
ax[0].set_ylabel('Number of leaves')
ax[0].set_xlabel('FWHM STD')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['fwhm_std']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('FWHM STD')
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['fwhm_std']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_fwhm_std.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['fwhm_av'].hist(ax=ax[0])
leaf_kin_qc['fwhm_av'].hist(ax=ax[1])
ax[0].set_ylabel('Number of leaves',fontsize=15)
ax[0].set_xlabel('FWHM Mean',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['fwhm_av']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('FWHM Mean',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['fwhm_av']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_fwhm_av.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['Vlsr_std'].hist(ax=ax[0])
leaf_kin_qc['Vlsr_std'].hist(ax=ax[1])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Vlsr STD',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_std']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('Vlsr STD',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_std']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_Vlsr_std.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
(leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']).hist(ax=ax[0])
(leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']).hist(ax=ax[1])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Vlsr difference',fontsize=15)
ax[0].title.set_text('Non-quality Controlled')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('Vlsr difference',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_Vlsr_diff.pdf')

plt.figure(figsize=(12,12))
plt.scatter((leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']),leaf_kin_qc['fwhm_std'])

#leaf_kin['fwhm_std'].hist()
#leaf_kin['fwhm_av'].hist()
#leaf_kin['Vlsr_std'].hist()
#(leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']).hist()

# CRUCIAL TO UNDERSTAND WHY THE VLSR VALUES ARE SO FAR OFF

# -----------------------
# Splitting data by transition and detection statistics for each transition
# 'C18O.219.6GHz', '12CO.230.5GHz', '13CO.220.4GHz', 'H2CO.218.2GHz',
# 'H2CO.218.5GHz', 'H2CO.218.8GHz', 'OCS.218.9GHz', 'OCS.231.1GHz',
# 'SO.219.9GHz', 'SiO.217.1GHz'
C18O       = data_qc[data_qc['Transition'] == "C18O.219.6GHz"]  # 78
ThCO       = data_qc[data_qc['Transition'] == "13CO.220.4GHz"]  #245, 247
H2CO_218_2 = data_qc[data_qc['Transition'] == "H2CO.218.2GHz"]  #120
H2CO_218_5 = data_qc[data_qc['Transition'] == "H2CO.218.5GHz"]  #128
H2CO_218_8 = data_qc[data_qc['Transition'] == "H2CO.218.8GHz"]  # 60
OCS_218_9  = data_qc[data_qc['Transition'] == "OCS.218.9GHz"]   # 25
OCS_231_1  = data_qc[data_qc['Transition'] == "OCS.231.1GHz"]   # 18
SO         = data_qc[data_qc['Transition'] == "SO.219.9GHz"]    # 82
SiO        = data_qc[data_qc['Transition'] == "SiO.217.1GHz"]   # 50, 49

fig,ax = plt.subplots(figsize=(16,16),nrows=3,ncols=3)
ax = ax.ravel()
#ax[0].hist(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])
#ax[0].title.set_text(r'C$^{18}$O')
#ax[0].set_aspect(1./ax[0].get_data_ratio())
#ax[0].text(0.77, 0.95, "N = "+str(len(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])), transform=ax[0].transAxes)
#fig.savefig('Documents/scouse/All_Vlsr_diff.pdf')

#C18O['fwhm_int'].hist()   # With instrumwent response removed
ax[1].hist(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])#ax[1].set_ylabel('Number of spectra',fontsize=15)
#ax[1].set_xlabel('FWHM',fontsize=15)
ax[1].title.set_text(r'$C^{18}$O')
ax[1].text(0.77, 0.95, "N = "+str(len(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])), transform=ax[1].transAxes)
#ax[1].set_aspect(1./ax[1].get_data_ratio())

ax[2].hist(H2CO_218_2['fwhm_int'][~np.isnan(H2CO_218_2['fwhm_int'])])
#ax[2].set_ylabel('Number of spectra',fontsize=15)
#ax[2].set_xlabel('FWHM',fontsize=15)
ax[2].title.set_text(r'H$_{2}$CO (218.2 GHz)')
ax[2].text(0.77, 0.95, "N = "+str(len(H2CO_218_2['fwhm_int'][~np.isnan(H2CO_218_2['fwhm_int'])])), transform=ax[2].transAxes)
#ax[2].set_aspect(1./ax[2].get_data_ratio())

ax[3].hist(H2CO_218_5['fwhm_int'][~np.isnan(H2CO_218_5['fwhm_int'])])
#ax[3].set_xlabel('FWHM',fontsize=15)
ax[3].title.set_text(r'H$_{2}$CO (218.5 GHz)')
ax[3].text(0.77, 0.95, "N = "+str(len(H2CO_218_5['fwhm_int'][~np.isnan(H2CO_218_5['fwhm_int'])])), transform=ax[3].transAxes)
#ax[3].set_aspect(1./ax[3].get_data_ratio())

ax[4].hist(H2CO_218_8['fwhm_int'][~np.isnan(H2CO_218_8['fwhm_int'])])#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_xlabel('FWHM',fontsize=15)
ax[4].title.set_text(r'H$_{2}$CO (218.8 GHz)')
ax[4].text(0.77, 0.95, "N = "+str(len(H2CO_218_8['fwhm_int'][~np.isnan(H2CO_218_8['fwhm_int'])])), transform=ax[4].transAxes)
#ax[4].set_aspect(1./ax[4].get_data_ratio())

ax[5].hist(OCS_218_9['fwhm_int'][~np.isnan(OCS_218_9['fwhm_int'])])#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_xlabel('FWHM',fontsize=15)
ax[5].title.set_text(r'OCS (218.9 GHz)')
ax[5].text(0.77, 0.95, "N = "+str(len(OCS_218_9['fwhm_int'][~np.isnan(OCS_218_9['fwhm_int'])])), transform=ax[5].transAxes)
#ax[5].set_aspect(1./ax[5].get_data_ratio())

ax[6].hist(OCS_231_1['fwhm_int'][~np.isnan(OCS_231_1['fwhm_int'])])#ax[3].set_ylabel('Number of spectra',fontsize=15)
ax[6].set_ylabel('Number of spectra',fontsize=15)
ax[6].set_xlabel('FWHM',fontsize=15)
ax[6].title.set_text(r'OCS (231.1 GHz)')
ax[6].text(0.77, 0.95, "N = "+str(len(OCS_231_1['fwhm_int'][~np.isnan(OCS_231_1['fwhm_int'])])), transform=ax[6].transAxes)
#ax[6].set_aspect(1./ax[6].get_data_ratio())

ax[7].hist(SiO['fwhm_int'][~np.isnan(SiO['fwhm_int'])])#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_xlabel('FWHM',fontsize=15)
ax[7].title.set_text('SiO')
ax[7].text(0.77, 0.95, "N = "+str(len(SiO['fwhm_int'][~np.isnan(SiO['fwhm_int'])])), transform=ax[7].transAxes)
#ax[7].set_aspect(1./ax[7].get_data_ratio())

ax[8].hist(SO['fwhm_int'][~np.isnan(SO['fwhm_int'])])#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_xlabel('FWHM',fontsize=15)
ax[8].title.set_text('SO')
ax[8].text(0.77, 0.95, "N = "+str(len(SO['fwhm_int'][~np.isnan(SO['fwhm_int'])])), transform=ax[8].transAxes)
#ax[8].set_aspect(1./ax[8].get_data_ratio())
fig.savefig('Documents/scouse/line_fwhms.pdf')

# =================================================
# Plots


# -----------------------
# Plotting RMS histogram

fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(16,12)) 
ax[0].hist(data['RMS'], log=True)
ax[1].hist(data_qc['RMS'],log=True) 
ax[0].set_xlabel('RMS')
ax[0].set_ylabel('Number of spectra')
ax[1].set_ylabel('Number of spectra')
ax[1].set_xlabel('RMS')
ax[0].grid()
ax[1].grid()
ax[0].title.set_text('Non-quality Controlled')
ax[1].title.set_text('Quality Controlled')
fig.savefig('Documents/scouse/RMS_all.pdf')

fig,ax = plt.subplots(figsize=(16,16),nrows=3,ncols=3)
ax = ax.ravel()
ax[0].hist(C18O['RMS'])
ax[0].title.set_text(r'C$^{18}$O')

ax[1].hist(ThCO['RMS'])
ax[1].title.set_text(r'$^{13}$CO')

ax[2].hist(H2CO_218_2['RMS'])
ax[2].title.set_text(r'H$_{2}$CO (218.2 GHz)')

ax[3].hist(H2CO_218_5['RMS'])
ax[3].title.set_text(r'H$_{2}$CO (218.5 GHz)')

ax[4].hist(H2CO_218_8['RMS'])
ax[4].title.set_text(r'H$_{2}$CO (218.8 GHz)')

ax[5].hist(OCS_218_9['RMS'])
ax[5].title.set_text(r'OCS (218.9 GHz)')

ax[6].hist(OCS_231_1['RMS'])
ax[6].title.set_text(r'OCS (231.1 GHz)')
ax[6].set_xlabel('RMS')
ax[6].set_ylabel('Number of Spectra')

ax[7].hist(SiO['RMS'])
ax[7].title.set_text('SiO')

ax[8].hist(SO['RMS'])
ax[8].title.set_text('SO')
fig.savefig('Documents/scouse/line_rms.pdf')

unique_all, counts_all       = np.unique(data['Leaf_ID'].values, return_counts=True)
unique_all_qc, counts_all_qc = np.unique(data_ID['Leaf_ID'].values, return_counts=True)

# creating variable and column names to hold kinematic info for each leaf
leaf_kin = pd.DataFrame(np.zeros([unique_all.size,10]))
leaf_kin.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max']

leaf_kin_qc = pd.DataFrame(np.zeros([unique_all_qc.size,10]))
leaf_kin_qc.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max']

# filling the leaf_kin variable with appropriate info for non qc
i=0
for id in unique_all:
    leaf_kin.iloc[i,0] = id
    leaf_kin.iloc[i,1] = data[data['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin.iloc[i,2] = data[data['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin.iloc[i,3] = data[data['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin.iloc[i,4] = data[data['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin.iloc[i,5] = data[data['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin.iloc[i,6] = data[data['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin.iloc[i,7] = data[data['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin.iloc[i,8] = data[data['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin.iloc[i,9] = data[data['Leaf_ID'] == id]['Vlsr'].max()
    i=i+1

i=0
for id in unique_all_qc:
    leaf_kin_qc.iloc[i,0] = id
    leaf_kin_qc.iloc[i,1] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin_qc.iloc[i,2] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin_qc.iloc[i,3] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin_qc.iloc[i,4] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin_qc.iloc[i,5] = data_qc[data_qc['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin_qc.iloc[i,6] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin_qc.iloc[i,7] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin_qc.iloc[i,8] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin_qc.iloc[i,9] = data_qc[data_qc['Leaf_ID'] == id]['Vlsr'].max()
    i=i+1
    
fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['Vlsr_std'].hist(ax=ax[0])
leaf_kin_qc['Vlsr_std'].hist(ax=ax[1])
ax[0].set_ylabel('Number of leaves',fontsize=15)
ax[0].set_xlabel('Vlsr STD',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_std']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('Vlsr STD',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_std']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_Vlsr_std_temp.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
(leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']).hist(ax=ax[0])
(leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']).hist(ax=ax[1])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Vlsr diff',fontsize=15)
ax[0].title.set_text('Non-quality Controlled')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']))),
           transform=ax[0].transAxes)
ax[1].set_xlabel('Vlsr diff',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']))),
           transform=ax[1].transAxes)
fig.savefig('Documents/scouse/All_Vlsr_diff_temp.pdf')

C18O       = data_ID[data_ID['Transition'] == "C18O.219.6GHz"]  # 78
H2CO_218_2 = data_ID[data_ID['Transition'] == "H2CO.218.2GHz"]  #120
H2CO_218_5 = data_ID[data_ID['Transition'] == "H2CO.218.5GHz"]  #128
H2CO_218_8 = data_ID[data_ID['Transition'] == "H2CO.218.8GHz"]  # 60
OCS_218_9  = data_ID[data_ID['Transition'] == "OCS.218.9GHz"]   # 25
OCS_231_1  = data_ID[data_ID['Transition'] == "OCS.231.1GHz"]   # 18
SO         = data_ID[data_ID['Transition'] == "SO.219.9GHz"]    # 82
SiO        = data_ID[data_ID['Transition'] == "SiO.217.1GHz"]   # 50, 49
CH3OH      = data_ID[data_ID['Transition'] == "CH3OH"]  # 78

lines = [C18O,H2CO_218_2,H2CO_218_5,H2CO_218_8,OCS_218_9,OCS_231_1,SO, SiO]
#lines = [H2CO_218_2]

plt.figure(figsize=(8,8))
x = leaf_kin_qc['fwhm_av']
y = (leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min'])
plt.scatter(y, x)

from astropy.io import fits
catalog = fits.open('catalog_acc-2.fits')
cat_cols = catalog[1].columns
cat_data = catalog[1].data

unique_l, unique_n = np.unique(data_qc['Leaf_ID'].values, return_counts=True)

vel_out = pd.DataFrame(columns=['id','l','b','v_peak','v_avg'])

#Method 1 - takes the velocity of the brightest component

ID_data = data_ID[data_ID['rest_freq'].isna() == True]
ID_data = ID_data[ID_data['Region'].isna() == False]
unique_l, unique_n = np.unique(ID_data['Leaf_ID'].values, return_counts=True)
i=0
empty_dict = []
for id in unique_l:
    print(cat_data[cat_data['leaf_ID'] == id],id)
    try:
        glat = cat_data[cat_data['leaf_ID'] == id]['glat'][0]
        glon = cat_data[cat_data['leaf_ID'] == id]['glon'][0]
        test_vels = np.array(ID_data[ID_data['Leaf_ID'] == id]['Vlsr'])
        test_peaks = np.array(ID_data[ID_data['Leaf_ID'] == id]['Amp'])
        print(id,test_peaks)
        if len(test_vels) > 0:
            peak_vel = np.nanmean(test_vels[test_peaks == np.amax(test_peaks)])
            print(peak_vel)
            test_errs = np.array(ID_data[ID_data['Leaf_ID'] == id]['Vlsr_err'])
            if len(test_vels) > 0:
                peak_avg = np.average(test_vels, weights=1./test_errs)
                temp_data = {'id': id, 'l': glat, 'b':glon, 'v_peak':peak_vel, 'v_avg':peak_avg}
                print(temp_data)
                empty_dict.append(temp_data)
    except IndexError:
        "No Leaf"
    #temp_df = pd.DataFrame(data=temp_data, index=[i])
    #i += 1
    #vel_out = vel_out.append(temp_df, ignore_index=True)

vel_out = pd.DataFrame(data=empty_dict)

IDs = data_ID[data_ID['Leaf_ID'].isna() == True].index
data_ID.drop(IDs,inplace=True)

for id in np.unique(data_ID['Leaf_ID']):
    if id[:-1] not in walker_id:
        subcube = data_ID[data_ID['Leaf_ID'] == id]
    else:
        continue
    for index,row in subcube.iterrows():
        print(row.Vlsr, np.array(vel_out[vel_out['id'] == id]['v_peak']))
        if np.abs(row.Vlsr-np.array(vel_out[vel_out['id'] == id]['v_peak'])) > 25:
            data_ID.loc[index,'original_cube'] = data_ID.loc[index,'Transition']
            data_ID.loc[index,'Transition'] = 'unknown'
        
unknowns = data_ID[data_ID['Transition'] == 'unknown']
IDs = data_ID[data_ID['Transition'] == 'unknown'].index

for index,row in unknowns.iterrows():
    vel = np.array(vel_out[vel_out['id'] == row.Leaf_ID]['v_peak'])
    unknowns.loc[index,'rest_freq']=unknowns.loc[index,'Frequency']/(1.0 - (vel/float(c.c.value))) 
    
for index,row in unknowns.iterrows():
    try:
        print(unknowns.loc[index,'Leaf_ID'],unknowns.loc[index,'Vlsr'])
        freq = unknowns.loc[index,'rest_freq']
        columns = ('Species','Chemical Name','Freq-GHz(rest frame,redshifted)','Meas Freq-GHz(rest frame,redshifted)',
                   'Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)')
        lines = Splatalogue.query_lines((freq-0.02)*u.GHz, (freq+0.02)*u.GHz,top20='top20',
                                        energy_max=100,energy_type='eu_k',only_NRAO_recommended=True,noHFS=True,
                                        line_lists=['JPL','CDMS'])[columns]
        lines.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
        lines.rename_column('E_U (K)','EU_K')
        lines.rename_column('Freq-GHz(rest frame,redshifted)','Freq-GHz')
        lines.rename_column('Meas Freq-GHz(rest frame,redshifted)','Meas Freq-GHz')
        lines.sort('log10(Aij)')
        lines.reverse()
        if len(lines) > 0:
            print(lines)
            number = int(input('Enter preferred row: '))
            if number > 0:
                unknowns.loc[index,'Transition'] = lines[number-1]['Species']+'.'+str(lines[number-1]['Freq-GHz'])
            else:
                print('No line selected.')
    except ConnectionError:
        continue
    
    
data_ID.drop(IDs, inplace=True)
data_ID = data_ID.append(unknowns)
data_ID.to_csv('Leaf_spectra_ID.csv')

'''
ONLY USE THIS FOR CORRELATION MATRIX
print(data_ID)
CO_data = data[(data['Transition'] == '12CO.230.5GHz')]
for index,row in CO_data.iterrows():
    data_ID = data_ID.append(row)
print(data_ID)
CO_data = data[(data['Transition'] == '13CO.220.4GHz')]
for index,row in CO_data.iterrows():
    data_ID = data_ID.append(row)
print(data_ID)

data_ID.to_csv('Corr_matrix.csv')
'''

from astropy.io import fits
catalog = fits.open('megacatalog_acc.fits')
cat_cols = catalog[1].columns
cat_data = catalog[1].data

unique_l, unique_n = np.unique(data_qc['Leaf_ID'].values, return_counts=True)

vel_out = pd.DataFrame(columns=['id','l','b','v_peak','v_avg'])

#Method 1 - takes the velocity of the brightest component
ID_data = pd.read_csv('Leaf_spectra_ID.csv')
ID_data['fwhm_int_err'] = ID_data.apply(lambda row: (row.fwhm/np.sqrt(row.fwhm**2 - CRPIX3**2)*row.fwhm_err), axis=1)
ID_data = ID_data[ID_data['rest_freq'].isna() == True]
ID_data = ID_data[ID_data['Region'].isna() == False]
unique_l, unique_n = np.unique(ID_data['Leaf_ID'].values, return_counts=True)
i=0
empty_dict = []
for id in unique_l:
    try:
        glat = cat_data[cat_data['leaf_ID'] == id]['glat'][0]
        glon = cat_data[cat_data['leaf_ID'] == id]['glon'][0]
        test_vels = np.array(ID_data[ID_data['Leaf_ID'] == id]['Vlsr'])
        test_peaks = np.array(ID_data[ID_data['Leaf_ID'] == id]['Amp'])
        if len(test_vels) > 0:
            peak_vel = np.nanmean(test_vels[test_peaks == np.amax(test_peaks)])
            test_errs = np.array(ID_data[ID_data['Leaf_ID'] == id]['Vlsr_err'])
        if len(test_vels) > 0:
            peak_avg = np.average(test_vels, weights=1./test_errs)
            vels = np.array(ID_data[(ID_data['Leaf_ID'] == id) & (ID_data['Transition'] == 'H2CO.218.2GHz')]['Vlsr'])
        if len(vels) > 0:
            vel = int(np.mean(vels))
        else:
            vel = np.nan
        fwhms = np.array(ID_data[(ID_data['Leaf_ID'] == id) & (ID_data['Transition'] == 'H2CO.218.2GHz')]['fwhm'])
        fwhm_ints = np.array(ID_data[(ID_data['Leaf_ID'] == id) & (ID_data['Transition'] == 'H2CO.218.2GHz')]['fwhm_int'])
        fwhm_errs = np.array(ID_data[(ID_data['Leaf_ID'] == id) & (ID_data['Transition'] == 'H2CO.218.2GHz')]['fwhm_err'])
        fwhm_int_errs = np.array(ID_data[(ID_data['Leaf_ID'] == id) & (ID_data['Transition'] == 'H2CO.218.2GHz')]['fwhm_int_err'])
        if len(fwhms) > 0:
            fwhm = np.around(np.mean(fwhms),2)
        else:
            fwhm = np.nan
        if len(fwhm_ints) > 0:
            fwhm_int = np.around(np.mean(fwhm_ints),2)
        else:
            fwhm_int = np.nan
        if len(fwhm_errs) > 0:
            fwhm_err = np.around(np.mean(fwhm_errs),2)
        else:
            fwhm_err = np.nan
        if len(fwhm_int_errs) > 0:
            fwhm_int_err = np.around(np.mean(fwhm_int_errs),2)
        else:
            fwhm_int_err = np.nan
        temp_data = {'id': id, 'l': glat, 'b':glon, 'v_peak':peak_vel, 'v_avg':peak_avg, 'h_vel':vel, 'h_fwhm':fwhm, 'h_fwhm_unc':fwhm_err,'h_fwhm_int':fwhm_int,'h_fwhm_int_unc':fwhm_int_err}
        empty_dict.append(temp_data)
    except IndexError:
        print(id,"No Leaf")
    #temp_df = pd.DataFrame(data=temp_data, index=[i])
    #i += 1
    #vel_out = vel_out.append(temp_df, ignore_index=True)

vel_out = pd.DataFrame(data=empty_dict)

vel_out.to_csv('Velocity_vals.csv')


unique_all, counts_all             = np.unique(data['Leaf_ID'].values,return_counts=True)
unique_all_qc, counts_all_qc       = np.unique(ID_data['Leaf_ID'].values, return_counts=True)

# creating variable and column names to hold kinematic info for each leaf
leaf_kin = pd.DataFrame(np.zeros([unique_all.size,14]))
leaf_kin_qc = pd.DataFrame(np.zeros([unique_all.size,14]))

leaf_kin.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max',
                     'amp_av', 'amp_std', 'amp_min', 'amp_max']
leaf_kin_qc.columns = ['Leaf_ID', 'count', 
                     'fwhm_av', 'fwhm_std', 'fwhm_min', 'fwhm_max',
                     'Vlsr_av', 'Vlsr_std', 'Vlsr_min', 'Vlsr_max',
                     'amp_av', 'amp_std', 'amp_min', 'amp_max']

i=0
for id in unique_all:
    leaf_kin.iloc[i,0] = id
    leaf_kin.iloc[i,1] = data[data['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin.iloc[i,2] = data[data['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin.iloc[i,3] = data[data['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin.iloc[i,4] = data[data['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin.iloc[i,5] = data[data['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin.iloc[i,6] = data[data['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin.iloc[i,7] = data[data['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin.iloc[i,8] = data[data['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin.iloc[i,9] = data[data['Leaf_ID'] == id]['Vlsr'].max()
    leaf_kin.iloc[i,10] = data[data['Leaf_ID'] == id]['Amp'].mean()
    leaf_kin.iloc[i,11] = data[data['Leaf_ID'] == id]['Amp'].std()
    leaf_kin.iloc[i,12] = data[data['Leaf_ID'] == id]['Amp'].min()
    leaf_kin.iloc[i,13] = data[data['Leaf_ID'] == id]['Amp'].max()
    i=i+1
    
i=0
for id in unique_all_qc:
    leaf_kin_qc.iloc[i,0] = id
    leaf_kin_qc.iloc[i,1] = ID_data[ID_data['Leaf_ID'] == id]['fwhm_int'].count()
    leaf_kin_qc.iloc[i,2] = ID_data[ID_data['Leaf_ID'] == id]['fwhm_int'].mean()
    leaf_kin_qc.iloc[i,3] = ID_data[ID_data['Leaf_ID'] == id]['fwhm_int'].std()
    leaf_kin_qc.iloc[i,4] = ID_data[ID_data['Leaf_ID'] == id]['fwhm_int'].min()
    leaf_kin_qc.iloc[i,5] = ID_data[ID_data['Leaf_ID'] == id]['fwhm_int'].max()
    leaf_kin_qc.iloc[i,6] = ID_data[ID_data['Leaf_ID'] == id]['Vlsr'].mean()
    leaf_kin_qc.iloc[i,7] = ID_data[ID_data['Leaf_ID'] == id]['Vlsr'].std()
    leaf_kin_qc.iloc[i,8] = ID_data[ID_data['Leaf_ID'] == id]['Vlsr'].min()
    leaf_kin_qc.iloc[i,9] = ID_data[ID_data['Leaf_ID'] == id]['Vlsr'].max()
    leaf_kin_qc.iloc[i,10] = ID_data[ID_data['Leaf_ID'] == id]['Amp'].mean()
    leaf_kin_qc.iloc[i,11] = ID_data[ID_data['Leaf_ID'] == id]['Amp'].std()
    leaf_kin_qc.iloc[i,12] = ID_data[ID_data['Leaf_ID'] == id]['Amp'].min()
    leaf_kin_qc.iloc[i,13] = ID_data[ID_data['Leaf_ID'] == id]['Amp'].max()
    i=i+1
    
leaf_kin_qc = leaf_kin_qc[leaf_kin_qc['Leaf_ID'] != 0]

# FWHM
fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
data['fwhm'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves')
ax[0].set_xlabel('Velocity Dispersion [km s$^{-1}$]')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len(data['fwhm'])),
           transform=ax[0].transAxes)
data_qc['fwhm'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves')
ax[1].set_xlabel('Velocity Dispersion [km s$^{-1}$]')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len(data_qc['fwhm'])),
           transform=ax[1].transAxes)
fig.savefig('All_fwhm.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['fwhm_std'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves')
ax[0].set_xlabel('Velocity Dispersion STD [km s$^{-1}$]')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len(leaf_kin['fwhm_std'])),
           transform=ax[0].transAxes)
leaf_kin_qc['fwhm_std'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves')
ax[1].set_xlabel('Velocity Dispersion STD [km s$^{-1}$]')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len(leaf_kin_qc['fwhm_std'])),
           transform=ax[1].transAxes)
fig.savefig('All_fwhm_std.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['fwhm_av'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves',fontsize=15)
ax[0].set_xlabel('Mean Velocity Dispersion [km s$^{-1}$]',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['fwhm_av']))),
           transform=ax[0].transAxes)
leaf_kin_qc['fwhm_av'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves',fontsize=15)
ax[1].set_xlabel('Mean Velocity Dispersion [km s$^{-1}$]',fontsize=15)
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['fwhm_av']))),
           transform=ax[1].transAxes)
fig.savefig('All_fwhm_av.pdf')

# AMP

fig,ax = plt.subplots(figsize=(16,8),nrows=1,ncols=2)
data['Amp'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves')
ax[0].set_xlabel('Peak Intensity [Jy/beam]')
#ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len(data['Amp'])),
           transform=ax[0].transAxes)
ax[0].set_yscale('log')
data_qc['Amp'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves')
ax[1].set_xlabel('Peak Intensity [Jy/beam]')
#ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len(data_qc['Amp'])),
           transform=ax[1].transAxes)
ax[1].set_yscale('log')
fig.savefig('All_amp.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['amp_std'].hist(ax=ax[0])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Peak Intensity STD [Jy/beam]',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['amp_std']))),
           transform=ax[0].transAxes)
leaf_kin_qc['amp_std'].hist(ax=ax[1])
ax[1].set_ylabel('Number of spectra',fontsize=15)
ax[1].set_xlabel('Peak Intensity STD [Jy/beam]',fontsize=15)
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['amp_std']))),
           transform=ax[1].transAxes)
fig.savefig('All_amp_std.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['amp_av'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves',fontsize=15)
ax[0].set_xlabel('Mean Peak Intensity [Jy/beam]',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['amp_av']))),
           transform=ax[0].transAxes)
leaf_kin_qc['amp_av'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves',fontsize=15)
ax[1].set_xlabel('Mean Peak Intensity [Jy/beam]',fontsize=15)
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['amp_av']))),
           transform=ax[1].transAxes)
fig.savefig('All_amp_av.pdf')

#Vlsr

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
data['Vlsr'].hist(ax=ax[0])
ax[0].set_ylabel('Number of leaves')
ax[0].set_xlabel('Line of sight velocity [km s$^{-1}$]')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len(data['Vlsr'])),
           transform=ax[0].transAxes)
data_qc['Vlsr'].hist(ax=ax[1])
ax[1].set_ylabel('Number of leaves')
ax[1].set_xlabel('Line of sight velocity [km s$^{-1}$]')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len(data_qc['Vlsr'])),
           transform=ax[1].transAxes)
fig.savefig('All_vel.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['Vlsr_std'].hist(ax=ax[0])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Line of sight velocity STD [km s$^{-1}$]',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_std']))),
           transform=ax[0].transAxes)
leaf_kin_qc['Vlsr_std'].hist(ax=ax[1])
ax[1].set_ylabel('Number of spectra',fontsize=15)
ax[1].set_xlabel('Line of sight velocity STD [km s$^{-1}$]',fontsize=15)
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_std']))),
           transform=ax[1].transAxes)
fig.savefig('All_Vlsr_std.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
leaf_kin['Vlsr_av'].hist(ax=ax[0])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Mean Line of sight velocity [km s$^{-1}$]',fontsize=15)
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].title.set_text('Non-quality Controlled')
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_av']))),
           transform=ax[0].transAxes)
leaf_kin_qc['Vlsr_av'].hist(ax=ax[1])
ax[1].set_ylabel('Number of spectra',fontsize=15)
ax[1].set_xlabel('Mean Line of sight velocity [km s$^{-1}$]',fontsize=15)
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].title.set_text('Quality Controlled')
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_av']))),
           transform=ax[1].transAxes)
fig.savefig('All_Vlsr_av.pdf')

fig,ax = plt.subplots(figsize=(16,12),nrows=1,ncols=2)
(leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']).hist(ax=ax[0])
ax[0].set_ylabel('Number of spectra',fontsize=15)
ax[0].set_xlabel('Maximum line of sight velocity difference [km s$^{-1}$]',fontsize=15)
ax[0].title.set_text('Non-quality Controlled')
ax[0].set_aspect(1./ax[0].get_data_ratio())
ax[0].text(0.8,0.95,"N = "+str(len((leaf_kin['Vlsr_max']-leaf_kin['Vlsr_min']))),
           transform=ax[0].transAxes)
(leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']).hist(ax=ax[1])
ax[1].set_ylabel('Number of spectra',fontsize=15)
ax[1].set_xlabel('Maximum line of sight velocity difference [km s$^{-1}$]',fontsize=15)
ax[1].title.set_text('Quality Controlled')
ax[1].set_aspect(1./ax[1].get_data_ratio())
ax[1].text(0.8,0.95,"N = "+str(len((leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']))),
           transform=ax[1].transAxes)
fig.savefig('All_Vlsr_diff.pdf')

plt.figure(figsize=(12,12))
plt.scatter((leaf_kin_qc['Vlsr_max']-leaf_kin_qc['Vlsr_min']),leaf_kin_qc['fwhm_std'])
#fig.savefig('Documents/scouse/All_fwhm_std.pdf')

C18O       = ID_data[ID_data['Transition'] == "C18O.219.6GHz"]  # 78
ThCO       = ID_data[ID_data['Transition'] == "13CO.220.4GHz"]  #245, 247
H2CO_218_2 = ID_data[ID_data['Transition'] == "H2CO.218.2GHz"]  #120
H2CO_218_5 = ID_data[ID_data['Transition'] == "H2CO.218.5GHz"]  #128
H2CO_218_8 = ID_data[ID_data['Transition'] == "H2CO.218.8GHz"]  # 60
OCS_218_9  = ID_data[ID_data['Transition'] == "OCS.218.9GHz"]   # 25
OCS_231_1  = ID_data[ID_data['Transition'] == "OCS.231.1GHz"]   # 18
SO         = ID_data[ID_data['Transition'] == "SO.219.9GHz"]    # 82
SiO        = ID_data[ID_data['Transition'] == "SiO.217.1GHz"]   # 50, 49

fig,ax = plt.subplots(figsize=(16,10),nrows=2,ncols=4)
ax = ax.ravel()
#ax[0].hist(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])
#ax[0].title.set_text(r'C$^{18}$O')
#ax[0].set_aspect(1./ax[0].get_data_ratio())
#ax[0].text(0.77, 0.95, "N = "+str(len(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])), transform=ax[0].transAxes)
#fig.savefig('Documents/scouse/All_Vlsr_diff.pdf')
#[~np.isnan(C18O['fwhm']
#C18O['fwhm_int'].hist()   # With instrumwent response removed
ax[0].hist(C18O['fwhm'])
#ax[1].set_ylabel('Number of spectra',fontsize=15)
#ax[1].set_xlabel('FWHM',fontsize=15)
ax[0].title.set_text(r'$C^{18}$O')
ax[0].text(0.7, 0.90, "N = "+str(len(C18O['fwhm'])), transform=ax[0].transAxes)
#ax[1].set_aspect(1./ax[1].get_data_ratio())

ax[1].hist(H2CO_218_2['fwhm'])
#ax[2].set_ylabel('Number of spectra',fontsize=15)
#ax[2].set_xlabel('FWHM',fontsize=15)
ax[1].title.set_text(r'H$_{2}$CO (218.2 GHz)')
ax[1].text(0.7, 0.90, "N = "+str(len(H2CO_218_2['fwhm'])), transform=ax[1].transAxes)
#ax[2].set_aspect(1./ax[2].get_data_ratio())

ax[2].hist(H2CO_218_5['fwhm'])
#ax[3].set_xlabel('FWHM',fontsize=15)
ax[2].title.set_text(r'H$_{2}$CO (218.5 GHz)')
ax[2].text(0.7, 0.90, "N = "+str(len(H2CO_218_5['fwhm'])), transform=ax[2].transAxes)
#ax[3].set_aspect(1./ax[3].get_data_ratio())

ax[3].hist(H2CO_218_8['fwhm'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_xlabel('FWHM',fontsize=15)
ax[3].title.set_text(r'H$_{2}$CO (218.8 GHz)')
ax[3].text(0.7, 0.90, "N = "+str(len(H2CO_218_8['fwhm'])), transform=ax[3].transAxes)
#ax[4].set_aspect(1./ax[4].get_data_ratio())

ax[4].hist(OCS_218_9['fwhm'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_xlabel('Velocity Dispersion [km s$^{-1}$]',fontsize=15)
ax[4].title.set_text(r'OCS (218.9 GHz)')
ax[4].text(0.7, 0.90, "N = "+str(len(OCS_218_9['fwhm'])), transform=ax[4].transAxes)
#ax[5].set_aspect(1./ax[5].get_data_ratio())

ax[5].hist(OCS_231_1['fwhm'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_xlabel('FWHM',fontsize=15)
ax[5].title.set_text(r'OCS (231.1 GHz)')
ax[5].text(0.7, 0.90, "N = "+str(len(OCS_231_1['fwhm'])), transform=ax[5].transAxes)
#ax[6].set_aspect(1./ax[6].get_data_ratio())

ax[6].hist(SiO['fwhm'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_xlabel('FWHM',fontsize=15)
ax[6].title.set_text('SiO')
ax[6].text(0.7, 0.90, "N = "+str(len(SiO['fwhm'])), transform=ax[6].transAxes)
#ax[7].set_aspect(1./ax[7].get_data_ratio())

ax[7].hist(SO['fwhm'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_xlabel('FWHM',fontsize=15)
ax[7].title.set_text('SO')
ax[7].text(0.7, 0.90, "N = "+str(len(SO['fwhm'])), transform=ax[7].transAxes)
#ax[8].set_aspect(1./ax[8].get_data_ratio())
fig.savefig('line_fwhms.pdf')

fig,ax = plt.subplots(figsize=(16,10),nrows=2,ncols=4)
ax = ax.ravel()
#ax[0].hist(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])
#ax[0].title.set_text(r'C$^{18}$O')
#ax[0].set_aspect(1./ax[0].get_data_ratio())
#ax[0].text(0.77, 0.95, "N = "+str(len(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])), transform=ax[0].transAxes)
#fig.savefig('Documents/scouse/All_Vlsr_diff.pdf')
#[~np.isnan(C18O['fwhm']
#C18O['fwhm_int'].hist()   # With instrumwent response removed
ax[0].hist(C18O['Amp'])
#ax[1].set_ylabel('Number of spectra',fontsize=15)
#ax[1].set_xlabel('FWHM',fontsize=15)
ax[0].title.set_text(r'$C^{18}$O')
ax[0].text(0.7, 0.90, "N = "+str(len(C18O['Amp'])), transform=ax[0].transAxes)
#ax[1].set_aspect(1./ax[1].get_data_ratio())

ax[1].hist(H2CO_218_2['Amp'])
#ax[2].set_ylabel('Number of spectra',fontsize=15)
#ax[2].set_xlabel('FWHM',fontsize=15)
ax[1].title.set_text(r'H$_{2}$CO (218.2 GHz)')
ax[1].text(0.7, 0.90, "N = "+str(len(H2CO_218_2['Amp'])), transform=ax[1].transAxes)
#ax[2].set_aspect(1./ax[2].get_data_ratio())

ax[2].hist(H2CO_218_5['Amp'])
#ax[3].set_xlabel('FWHM',fontsize=15)
ax[2].title.set_text(r'H$_{2}$CO (218.5 GHz)')
ax[2].text(0.7, 0.90, "N = "+str(len(H2CO_218_5['Amp'])), transform=ax[2].transAxes)
#ax[3].set_aspect(1./ax[3].get_data_ratio())

ax[3].hist(H2CO_218_8['Amp'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_xlabel('FWHM',fontsize=15)
ax[3].title.set_text(r'H$_{2}$CO (218.8 GHz)')
ax[3].text(0.7, 0.90, "N = "+str(len(H2CO_218_8['Amp'])), transform=ax[3].transAxes)
#ax[4].set_aspect(1./ax[4].get_data_ratio())

ax[4].hist(OCS_218_9['Amp'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_xlabel('Peak Intensity [Jy/beam]',fontsize=15)
ax[4].title.set_text(r'OCS (218.9 GHz)')
ax[4].text(0.7, 0.90, "N = "+str(len(OCS_218_9['Amp'])), transform=ax[4].transAxes)
#ax[5].set_aspect(1./ax[5].get_data_ratio())

ax[5].hist(OCS_231_1['Amp'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_xlabel('FWHM',fontsize=15)
ax[5].title.set_text(r'OCS (231.1 GHz)')
ax[5].text(0.7, 0.90, "N = "+str(len(OCS_231_1['Amp'])), transform=ax[5].transAxes)
#ax[6].set_aspect(1./ax[6].get_data_ratio())

ax[6].hist(SiO['Amp'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_xlabel('FWHM',fontsize=15)
ax[6].title.set_text('SiO')
ax[6].text(0.7, 0.90, "N = "+str(len(SiO['Amp'])), transform=ax[6].transAxes)
#ax[7].set_aspect(1./ax[7].get_data_ratio())

ax[7].hist(SO['Amp'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_xlabel('FWHM',fontsize=15)
ax[7].title.set_text('SO')
ax[7].text(0.7, 0.90, "N = "+str(len(SO['Amp'])), transform=ax[7].transAxes)
#ax[8].set_aspect(1./ax[8].get_data_ratio())
fig.savefig('line_amps.pdf')

fig,ax = plt.subplots(figsize=(16,10),nrows=2,ncols=4)
ax = ax.ravel()
#ax[0].hist(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])
#ax[0].title.set_text(r'C$^{18}$O')
#ax[0].set_aspect(1./ax[0].get_data_ratio())
#ax[0].text(0.77, 0.95, "N = "+str(len(C18O['fwhm_int'][~np.isnan(C18O['fwhm_int'])])), transform=ax[0].transAxes)
#fig.savefig('Documents/scouse/All_Vlsr_diff.pdf')
#[~np.isnan(C18O['fwhm']
#C18O['fwhm_int'].hist()   # With instrumwent response removed
ax[0].hist(C18O['Vlsr'])
#ax[1].set_ylabel('Number of spectra',fontsize=15)
#ax[1].set_xlabel('FWHM',fontsize=15)
ax[0].title.set_text(r'$C^{18}$O')
ax[0].text(0.05, 0.90, "N = "+str(len(C18O['Vlsr'])), transform=ax[0].transAxes)
#ax[1].set_aspect(1./ax[1].get_data_ratio())

ax[1].hist(H2CO_218_2['Vlsr'])
#ax[2].set_ylabel('Number of spectra',fontsize=15)
#ax[2].set_xlabel('FWHM',fontsize=15)
ax[1].title.set_text(r'H$_{2}$CO (218.2 GHz)')
ax[1].text(0.05, 0.90, "N = "+str(len(H2CO_218_2['Vlsr'])), transform=ax[1].transAxes)
#ax[2].set_aspect(1./ax[2].get_data_ratio())

ax[2].hist(H2CO_218_5['Vlsr'])
#ax[3].set_xlabel('FWHM',fontsize=15)
ax[2].title.set_text(r'H$_{2}$CO (218.5 GHz)')
ax[2].text(0.05, 0.90, "N = "+str(len(H2CO_218_5['Vlsr'])), transform=ax[2].transAxes)
#ax[3].set_aspect(1./ax[3].get_data_ratio())

ax[3].hist(H2CO_218_8['Vlsr'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_ylabel('Number of spectra',fontsize=15)
#ax[4].set_xlabel('FWHM',fontsize=15)
ax[3].title.set_text(r'H$_{2}$CO (218.8 GHz)')
ax[3].text(0.05, 0.90, "N = "+str(len(H2CO_218_8['Vlsr'])), transform=ax[3].transAxes)
#ax[4].set_aspect(1./ax[4].get_data_ratio())

ax[4].hist(OCS_218_9['Vlsr'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_ylabel('Number of spectra',fontsize=15)
ax[4].set_xlabel('Line of sight velocity [km s$^{-1}$]',fontsize=15)
ax[4].title.set_text(r'OCS (218.9 GHz)')
ax[4].text(0.05, 0.90, "N = "+str(len(OCS_218_9['Vlsr'])), transform=ax[4].transAxes)
#ax[5].set_aspect(1./ax[5].get_data_ratio())

ax[5].hist(OCS_231_1['Vlsr'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_ylabel('Number of spectra',fontsize=15)
#ax[5].set_xlabel('FWHM',fontsize=15)
ax[5].title.set_text(r'OCS (231.1 GHz)')
ax[5].text(0.05, 0.90, "N = "+str(len(OCS_231_1['Vlsr'])), transform=ax[5].transAxes)
#ax[6].set_aspect(1./ax[6].get_data_ratio())

ax[6].hist(SiO['Vlsr'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_ylabel('Number of spectra',fontsize=15)
#ax[7].set_xlabel('FWHM',fontsize=15)
ax[6].title.set_text('SiO')
ax[6].text(0.05, 0.90, "N = "+str(len(SiO['Vlsr'])), transform=ax[6].transAxes)
#ax[7].set_aspect(1./ax[7].get_data_ratio())

ax[7].hist(SO['Vlsr'])
#ax[3].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_ylabel('Number of spectra',fontsize=15)
#ax[8].set_xlabel('FWHM',fontsize=15)
ax[7].title.set_text('SO')
ax[7].text(0.05, 0.90, "N = "+str(len(SO['Vlsr'])), transform=ax[7].transAxes)
#ax[8].set_aspect(1./ax[8].get_data_ratio())
fig.savefig('line_vlsrs.pdf')
# =================================================
# Plots


# -----------------------
# Plotting RMS histogram

def Jy_to_K(val):
    bmaj = 3*u.arcsec
    bmin = 3*u.arcsec
    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
    freq = 230*u.GHz
    equiv = u.brightness_temperature(freq)
    return (val*u.Jy/beam_area).to(u.K, equivalencies=equiv)
    
    
fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(16,8)) 
ax[0].hist(data['RMS'], log=True)
ax[1].hist(ID_data['RMS'],log=True) 
ax[0].set_xlabel('RMS [Jy / beam]')
ax[0].set_ylabel('Number of spectra')
ax[1].set_ylabel('Number of spectra')
ax[1].set_xlabel('RMS [Jy / beam]')
ax[0].grid()
ax[1].grid()
ax[0].title.set_text('Non-quality Controlled')
ax[1].title.set_text('Quality Controlled')
fig.savefig('RMS_all.pdf')

fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(16,8)) 
ax[0].hist(Jy_to_K(data['RMS'].to_numpy()).value, log=True)
ax[1].hist(Jy_to_K(ID_data['RMS'].to_numpy()).value,log=True) 
ax[0].set_xlabel('RMS [K]')
ax[0].set_ylabel('Number of spectra')
ax[1].set_ylabel('Number of spectra')
ax[1].set_xlabel('RMS [K]')
ax[0].grid()
ax[1].grid()
ax[0].title.set_text('Non-quality Controlled')
ax[1].title.set_text('Quality Controlled')
fig.savefig('RMS_all_K.pdf')

fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(16,8)) 
ax[0].hist(data['SNR'], log=True)
ax[1].hist(ID_data['SNR'],log=True) 
ax[0].set_xlabel('RMS')
ax[0].set_ylabel('Number of spectra')
ax[1].set_ylabel('Number of spectra')
ax[1].set_xlabel('RMS')
ax[0].grid()
ax[1].grid()
ax[0].title.set_text('Non-quality Controlled')
ax[1].title.set_text('Quality Controlled')
fig.savefig('RMS_all.pdf')

fig,ax = plt.subplots(figsize=(16,10),nrows=2,ncols=4)
ax = ax.ravel()
ax[0].hist(C18O['RMS'])
ax[0].title.set_text(r'C$^{18}$O')

#ax[1].hist(ThCO['RMS'])
#ax[1].title.set_text(r'$^{13}$CO')

ax[1].hist(H2CO_218_2['RMS'])
ax[1].title.set_text(r'H$_{2}$CO (218.2 GHz)')

ax[2].hist(H2CO_218_5['RMS'])
ax[2].title.set_text(r'H$_{2}$CO (218.5 GHz)')

ax[3].hist(H2CO_218_8['RMS'])
ax[3].title.set_text(r'H$_{2}$CO (218.8 GHz)')

ax[4].hist(OCS_218_9['RMS'])
ax[4].title.set_text(r'OCS (218.9 GHz)')
ax[4].set_xlabel('RMS [Jy / beam]')
ax[4].set_ylabel('Number of Spectra')

ax[5].hist(OCS_231_1['RMS'])
ax[5].title.set_text(r'OCS (231.1 GHz)')


ax[6].hist(SiO['RMS'])
ax[6].title.set_text('SiO')

ax[7].hist(SO['RMS'])
ax[7].title.set_text('SO')
fig.savefig('line_rms.pdf')
# xs = 8             # Fig size x-axis.
# ys = 6             # Fig size y-axis.
# edgepad = 0.15     # Size of padding at the edges of the sub-plots 
# textpad = 0.03     # Size of padding around text
# msize = 5          # Marker size

# dx = 0.02                                        # Slight x-shift making sure vertical error bars don't overlap
# xmin = 0.5                                       # Setting x and y ranges
# xmax = np.max(counts_all)+0.5
# ymin = 0.3
# ymax = 10.**np.ceil(np.log10(np.max(nsys_all)))

# fig, ax = plt.subplots(1, 1, figsize=(xs, ys))
# ax.errorbar(np_all,     nsys_all,  yerr=np.sqrt(nsys_all),  fmt='o', c='k', ls='-', ms=msize, label='All systems')
# ax.errorbar(np_low-dx,  nsys_low,  yerr=np.sqrt(nsys_low),  fmt='o', c='b', ls='-', ms=msize, label='Field')
# ax.errorbar(np_high+dx, nsys_high, yerr=np.sqrt(nsys_high), fmt='o', c='r', ls='-', ms=msize, label='Overdensities')
# ax.set_yscale('log')
# ax.set_xlim((xmin,xmax))
# ax.set_ylim((ymin,ymax))
# ax.set_xlabel('Planets per system')
# ax.set_ylabel('Observed systems')
# ax.legend(loc=1, frameon=0)
# fig.subplots_adjust(left=edgepad, right=1-edgepad, top=1-edgepad, bottom=edgepad)
# figsave = plt.gcf()
# figsave.set_size_inches(xs,ys)
# plt.savefig('figures/planets_per_system.pdf',bbox_inches='tight',dpi=300)

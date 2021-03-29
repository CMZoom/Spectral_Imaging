import pandas as pd
import os
import astropy.constants as c
import astropy.units as u
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import csv
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

try:
    import Image
except ImportError:
    from PIL import Image
import imf

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

onesix_deg = ['G1.683-0.089', 'G1.670-0.130', 'G1.651-0.050', 'G1.602+0.018']
oneone_deg = ['G1.085-0.027', 'G1.038-0.074', 'G0.891-0.048']
Sgr_B2 = ['G0.714-0.100', 'G0.699-0.028', 'G0.619+0.012']
cloudef = ['G0.489+0.010']
cloudd = ['G0.412+0.052']
cloudc = ['G0.380+0.050']
cloudb = ['G0.340+0.055']
#Dust_Ridge = ['G0.489+0.010','G0.412+0.052','G0.380+0.050','G0.340+0.055']
three_pigs = ['G0.253+0.016','G0.145-0.086','G0.106-0.082','G0.068-0.075']
Arches = ['G0.054+0.027','G0.014+0.021']
fifty = ['G0.001-0.058']
twenty = ['G359.889-0.093']
stream = ['G359.648-0.133']
#fifty_twenty_cloud = ['G0.001-0.058','G359.889-0.093']
iso_HMSF = ['G0.393-0.034','G0.316-0.201','G0.212-0.001','G359.615-0.243','G359.137+0.031']
FSC = ['G0.326-0.085','G359.865+0.022','G359.734+0.002','G359.611+0.018']
CND = ['G359.948-0.052']
Apex = ['G0.070-0.035']

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def ID_region(sourcename):
    if sourcename in onesix_deg:
        region = '1.6 Degree Cloud'
        color='k'
        marker='x'
    elif sourcename in oneone_deg:
        region = '1.1 Degree Cloud'
        color='b'
        marker='+'
    elif sourcename in Sgr_B2:
        region = 'Sagittarius B2'
        color='r'
        marker='o'
    elif sourcename in cloudef:
        region = 'Cloud E/F'
        color='c'
        marker='^'
    elif sourcename in cloudd:
        region = 'Cloud D'
        color='steelblue'
        marker='^'
    elif sourcename in cloudc:
        region = 'Cloud C'
        color='red'
        marker='^'
    elif sourcename in cloudb:
        region = 'Cloud B'
        color='purple'
        marker='^'
    elif sourcename in three_pigs:
        region = 'Three Little Pigs'
        color='m'
        marker='x'
    elif sourcename in Arches:
        region = 'Arches'
        color='olive'
        marker='x'
    elif sourcename in fifty:
        region = '50 km s$^{-1}$ cloud'
        color='olive'
        marker='v'
    elif sourcename in twenty:
        region = '20 km s$^{-1}$ cloud'
        color='orange'
        marker='.'
    elif sourcename in iso_HMSF:
        region = 'Isolated HMSF Region'
        color='lime'
        marker='D'
    elif sourcename in FSC:
        region = 'Far-side Stream Candidate'
        color='brown'
        marker='+'
    elif sourcename in CND:
        region = 'Circumnuclear Disk'
        color='pink'
        marker='+'
    elif sourcename in Apex:
        region='Apex H$_{2}$CO Bridge'
        color='grey'
        marker='x'
    else:
        print(sourcename)
        region = 'None'
        marker='x'
        color='magenta'
    return color, marker, region

def SF_rate(id_arr):
    sfe = 0.1
    stars_per_year = 0
    starmass_per_year = 0
    freefreefrac = 0.5
    
    masses = []
    tff = []
    mask_num = []
    for id in id_arr:
        SF_t = mast_data[mast_data['leaf_ID'] == id]['SF_any_certain']
        masses.append(SF_t*mast_data[mast_data['leaf_ID'] == id]['mass'])
        tff.append(mast_data[mast_data['leaf_ID'] == id]['tff'])
        mask_num.append(mast_data[mast_data['leaf_ID'] == id]['mask_num'])
        
    for i in range(len(masses)):
        if mask_num[i]==9:
            cluster_for_leafi = imf.make_cluster(masses[i]*sfe*freefreefrac,silent=True,massfunc='kroupa')
        else:
            cluster_for_leafi = imf.make_cluster(masses[i]*sfe,silent=True,massfunc='kroupa')
            num_stars_for_leafi = len(cluster_for_leafi)
            starmass_for_leafi = sum(cluster_for_leafi)
            stars_per_yeari = num_stars_for_leafi/tff[i]
            starmass_per_yeari = starmass_for_leafi/tff[i]
            stars_per_year += stars_per_yeari
            starmass_per_year += starmass_per_yeari
    print("SFR over freefall time for each leaf, with SFE%1.1f"%sfe)
    print("Stars per year: %2.3f"%stars_per_year)
    print("SFR: %2.3f msun per year"%starmass_per_year)
    return(sfe,stars_per_year,starmass_per_year)

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique),fontsize=10,ncol=2,loc='lower center')
    
os.chdir('../Documents/scouse/')
vals = pd.read_csv('velocity_vals.csv')
vals = vals[vals['h_fwhm'].isna() == False]

catalog = fits.open('catalog_acc-2.fits')
cat_cols = catalog[1].columns
cat_data = catalog[1].data

new_catalog = fits.open('catalog_acc-2.fits')
new_cat_cols = new_catalog[1].columns
new_cat_data = new_catalog[1].data

master_cat = fits.open('megacatalog_acc.fits')
mast_cols = master_cat[1].columns
mast_data = master_cat[1].data

ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093a')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093b')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093c')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093d')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093h')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093i')[0][0]
mast_data[ID]['SF_any_certain'] = 1
ID = np.where(mast_data['leaf_ID'] == 'G359.889-0.093m')[0][0]
mast_data[ID]['SF_any_certain'] = 1


G       = c.G.value
k       = c.k_B.value
pc      = c.pc.value
Gamma   = 0.73

    
CRPIX = 1.057309424046252
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8))

x_vals = []
y_vals = []
ids = []
SF_id_all = []
SF_id_CMZ = []
iso_cores = []
for id in np.unique(vals['id']):
    #print(id[:-1])
    #c,m,r = ID_region(id[:-1])
    fwhm = vals[vals['id'] == id]['h_fwhm'].to_numpy()
    fwhm = fwhm*(u.km/u.s)
    R = new_cat_data[new_cat_data['leaf_ID'] == id]['r_eff_pc']*(u.pc)
    sigma = new_cat_data[new_cat_data['leaf_ID'] == id]['Sigma']*(u.cm**-2*u.g)
    y = np.log10((fwhm**2/R).value)
    x = np.log10(sigma.value)
    #print(x,y,id)
    #ax2.scatter(np.log10(R.value),np.log10(fwhm.value),marker='x')
    #ax.scatter(x,y,marker=m,color=c,label=r)
    ids.append(id)
    if id[:-1] in iso_HMSF:
        iso_cores.append(id)
    x_vals.append(x[0])
    y_vals.append(y[0])
    if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1 and id[:-1] not in iso_HMSF:
        marker = 'o'
        SF_id_all.append(id)
        SF_id_CMZ.append(id)
    elif mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1 and id[:-1] in iso_HMSF:
        marker = 'o'
        SF_id_all.append(id)
    elif mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1 and id[:-1] not in iso_HMSF:
        marker='s'
        SF_id_all.append(id)
        SF_id_CMZ.append(id)
    elif mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1 and id[:-1] in iso_HMSF:
        marker='s'
        SF_id_all.append(id)
    else:
        if id[:-1] in iso_HMSF:
            print(id)
        marker='x'
    c,m,r = ID_region(id[:-1])
    #ax.scatter(x,y,marker=marker,color=c,label=r)
    #if id[:-1] not in iso_HMSF:
        #ax.scatter(x,y,marker=marker,color=c,label=r
        
    if fwhm.value < CRPIX:
        ax.scatter(x,y,marker='x',color='gray',label='CMZoom Cores ($\sigma < \sigma_{int})$')
    elif fwhm.value < 1.5*CRPIX and fwhm.value > CRPIX:
        ax.scatter(x,y,marker='x',color='red',label='CMZoom Cores ($X < \sigma < 1.5\sigma_{int})$')
    else:
        ax.scatter(x,y,marker='x',color='black',label='CMZoom Cores ($\sigma > 1.5\sigma_{int}$)')
'''
unique_ids = [np.unique([x[:-1] for x in ids])][0]

patches = []
for id in unique_ids:
    print(id)
    c,m,r = ID_region(id)
    cores = [x for x in ids if x[:-1] == id]
    loc = [np.where(np.array(ids) == x)[0][0] for x in cores]
    X = [x_vals[x] for x in loc]
    X = X[(~np.isnan(X)]
    Y = [y_vals[x] for x in loc]
    if ~np.all(np.isnan(X)) and ~np.all(np.isnan(Y)):
        minx = np.nanmin(X)
        minxy = Y[np.where(X == minx)[0][0]]
        maxx = np.nanmax(X)
        maxxy = Y[np.where(X == maxx)[0][0]]
        miny = np.nanmin(Y)
        minyx = X[np.where(Y == miny)[0][0]]
        maxy = np.nanmax(Y)
        maxyx = X[np.where(Y == maxy)[0][0]]
        if ~np.isnan(minx) and ~np.isnan(maxx) and ~np.isnan(miny) and ~np.isnan(maxy):
            print(c)
            points = [[minx,minxy],[minyx,miny],[maxx,maxxy],[maxyx,maxy]]
            if id not in iso_HMSF:
                polygon = Polygon(points,closed=True,color=c,edgecolor='None',alpha=0.25)
            if id == 'G359.889-0.093':
                print(X,Y)
            ax.add_patch(polygon)
'''

N = 3
d = 0.05
x = int((N/d)+1)

fig2,ax2 = plt.subplots(nrows=1,ncols=2,figsize=(16,8))
for j in range(x):
    n = 0
    Y = np.zeros(shape=[5000])
    Sig = np.zeros(shape=[5000])
    Pe = 0
    for sigma in np.logspace(-2.5,1.5,num=5000):
        Y[n] = (1/3)*((np.pi*Gamma*G*sigma*10) + ((4*(Pe*(k*1e6)))/(sigma*10)))
        Sig[n] = sigma
        n = n+1
    if j == 0:
        ax.plot(np.log10(Sig), np.log10((Y*pc)/1e6)+(j*d), color='k', ls='--', alpha=0.4,zorder=1)
    SF_id = []
    for i in range(len(x_vals)):
        sig_id = np.where(np.log10(Sig) == find_nearest(np.log10(Sig),x_vals[i]))[0]
        if y_vals[i] < (np.log10((Y*pc)/1e6)+(j*d))[sig_id]:
            SF_id.append(ids[i])
            if j == 0:
                print(ids[i])
    # isolated only
    D = 0
    Dn = 0
    SF_ids = [x for x in SF_id if x[:-1] in iso_HMSF]
    for id in SF_id:
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1:
            D += 1
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1:
            Dn += 1
    tot = D + Dn
    ax2[1].scatter((j*d),D/len(SF_id),marker='x',color='red',label='All - Definite SF')
    ax2[1].scatter((j*d),tot/len(SF_id),marker='x',color='black',label='All - Definite + Possible SF')
    
    # CMZ Only
    D = 0
    Dn = 0
    SF_ids = [x for x in SF_id if x[:-1] not in iso_HMSF]
    for id in SF_ids:
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1:
            D += 1
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1:
            Dn += 1
    tot = D + Dn
    ax2[1].scatter((j*d),D/len(SF_ids),marker='x',color='purple',label='CMZ - Definite SF')
    ax2[1].scatter((j*d),tot/len(SF_ids),marker='x',color='blue',label='CMZ - Definite + Possible SF')
    
ax2[1].set_xlabel('Maximum Distance above P$_{e} = 0$ line')
#ax2[0].set_ylabel('Fraction of star forming cores')
ax2[1].set_ylim([0,1])
#fig2.savefig('F_SF.pdf')
e,spy,smpy=SF_rate(SF_id_CMZ)
#ax.text(-0.3,3.2,'SFR = %2.3f M$_{\odot}$ yr$^{-1}$'%smpy[0])
legend_without_duplicate_labels(ax2[1])  


print('Pressure =',Pe)
SF_rate(SF_id)

n = 0
Pe = [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]
for j in range(0,len(Pe),1):
    for sigma in np.logspace(-2.5,1.5,num=5000):
        Y[n] = (1/3)*((np.pi*Gamma*G*sigma*10) + ((4*(Pe[j]*(k*1e6)))/(sigma*10)))
        Sig[n] = sigma
        n = n+1
        
    ax.plot(np.log10(Sig), np.log10((Y*pc)/1e6), color='k', alpha=0.4,zorder=1)
    SF_id = []
    for i in range(len(x_vals)):
        sig_id = np.where(np.log10(Sig) == find_nearest(np.log10(Sig),x_vals[i]))[0]
        if y_vals[i] < np.log10((Y*pc)/1e6)[sig_id]:
            SF_id.append(ids[i])
    e,spy,smpy = SF_rate(SF_id)
    
    #isolated only
    
    D = 0
    Dn = 0
    SF_ids = [x for x in SF_id if x[:-1] in iso_HMSF]
    for id in SF_id:
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1:
            D += 1
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1:
            Dn += 1
    tot = D + Dn
    ax2[0].scatter(Pe[j],D/len(SF_id),marker='x',color='red')
    ax2[0].scatter(Pe[j],tot/len(SF_id),marker='x',color='black')
    
    #CMZ only
    D = 0
    Dn = 0
    SF_ids = [x for x in SF_id if x[:-1] not in iso_HMSF]
    for id in SF_ids:
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 1:
            D += 1
        if mast_data[mast_data['leaf_ID'] == id]['SF_any_certain'] == 0 and mast_data[mast_data['leaf_ID'] == id]['SF_any_all'] == 1:
            Dn += 1
    tot = D + Dn
    ax2[0].scatter(Pe[j],D/len(SF_ids),marker='x',color='purple')
    ax2[0].scatter(Pe[j],tot/len(SF_ids),marker='x',color='blue')
    n=0
ax2[0].set_xscale('log')   
ax2[0].set_xlabel('Upper Limit on External Pressure')
ax2[0].set_ylabel('Fraction of star forming cores')
ax2[0].set_ylim([0,1])
fig2.savefig('P_f_SF_old.pdf')
ax.set_xlim([-2.3,0.8])
ax.set_ylim([-2,3.5])
ax.set_xlabel('log $\Sigma$ (g cm$^{-2}$)')
ax.set_ylabel('log $\sigma^{2}$/R (km$^{2}$ s$^{-2}$ pc$^{-1}$)')

ax.text(-2.28,1.1,'P$_{e}$/k = 10$^{6}$ K cm$^{-3}$',fontsize=10)
ax.text(-2.2,2,'P$_{e}$/k = 10$^{7}$ K cm$^{-3}$',fontsize=10)
ax.text(-2.12,2.9,'P$_{e}$/k = 10$^{8}$ K cm$^{-3}$',fontsize=10)


with open('disc_clouds.txt') as inf:
    gal = np.loadtxt(inf)
    gal_x = gal[:,0]
    gal_y = gal[:,1]
    #gal_x = list(zip(*reader))[0]
    #gal_y = list(zip(*reader))[1]

ax.scatter(gal_x,gal_y,marker='+',color='black',s=10,linewidths=0.5,label='FBK (2011)')


X = [-0.62087078, 0.11148861, 0.08395594]
Y = [1.32944709, 1.25647033, 1.38995607]
ax.scatter(X,Y,marker='*',color='green',label=r'Dust Ridge Clouds')

X = [-0.02655597]
Y = [1.15157407]
ax.scatter(X,Y,marker='*',color='green',label='Dust Ridge Clouds',zorder=2)

X = [-0.3067933]
Y = [1.17524341]
ax.scatter(X,Y,marker='*',color='green',label='Dust Ridge Clouds')

X = [-0.47954163]
Y = [1.32944709]
ax.scatter(X,Y,marker='*',color='green',label='Dust Ridge Clouds')

legend_without_duplicate_labels(ax)
'''
ax.text(-0.8,-0.5,'1.6 Degree Cloud',color='k',fontsize=10)
ax.text(-0.8,-0.7,'1.1 Degree Cloud',color='b',fontsize=10)
ax.text(-0.8,-0.9,'Sagittarius B2',color='r',fontsize=10)
ax.text(-0.8,-1.1,'Cloud E/F',color='c',fontsize=10)
ax.text(-0.8,-1.3,'Cloud D',color='steelblue',fontsize=10)
ax.text(-0.8,-1.5,'Cloud C',color='red',fontsize=10)
ax.text(-0.8,-1.7,'Cloud B',color='purple',fontsize=10)
ax.text(-0.8,-1.9,'Three Little Pigs',color='m',fontsize=10)
#ax.text(-0.1,-0.5,'Arches',color='olive',fontsize=10)
ax.text(-0.1,-0.5,'50 km s$^{-1}$ cloud',color='olive',fontsize=10)
ax.text(-0.1,-0.7,'20 km s$^{-1}$ cloud',color='orange',fontsize=10)
#ax.text(-0.1,-0.9,'Isolated HMSF Region',color='lime',fontsize=10)
ax.text(-0.1,-1.1,'Far-side Stream Candidate',color='brown',fontsize=10)
ax.text(-0.1,-1.3,'Circumnuclear Disk',color='pink',fontsize=10)
ax.text(-0.1,-1.5,'Apex H$_{2}$CO Bridge',color='gray',fontsize=10)
'''

fig.savefig('Pressure.pdf')

version = '2.0' 

''' 
Welcome to CASTOR! This software aims at estimating general parameters of Core Collapse Supernovae by applying 
standard physical assumptions. We make use of Gaussian Processes to interpolate light curves and create 
synthetic spectra. We call "supernova of study" the object you want to retrieve data of, and "reference-supernova"
the most resembling supernova out of the training set of CASTOR. 

Steps:
		- Extraction of data 
		- Comparison with the training set and election of the reference supernova
		- Building synthetic spectra 
		- Estimation of parameters 
		- Plotting and printing results 


Data employed: 
		- light curves of the supernova-of-study 
		- (optional) spectra of the supernova-of-study
		- light curves of the reference-supernova
		- (optional) spectra of the reference supernova

Plots:
		- comparison.png    : plot of the comparison 
		- light_curves.png	: plot of the interpolated light curves
		- spectra.png       : plot of the synthetic spectra 
		- velocity.png      : plot of the velocity of the ejecta 
		- luminosity.png    : plot of bol and pseudo luminosity 
        - photosphere.png   : plot of photospheric parameters
''' 


import numpy as np 
import pandas as pd 
import os
import sys
from   collections     import defaultdict 
from   astropy.io      import ascii
from   scipy.stats     import chisquare, linregress
from   scipy.optimize  import curve_fit
from   scipy.signal    import find_peaks
from   scipy.integrate import trapezoid  as trapz
import matplotlib.pyplot as plt 
import matplotlib.cm     as cm
import matplotlib.colors as mcolors
import warnings
from   sklearn.gaussian_process         import GaussianProcessRegressor
from   sklearn.gaussian_process.kernels import Matern, ConstantKernel as C, WhiteKernel
from   joblib         import Parallel, delayed
from   threadpoolctl  import ThreadpoolController
from   functools       import partial
import time   as TT
warnings.filterwarnings('ignore')


# print logo 

print(r'##############################################')
print(r'#   ____    _    ____ _____ ___  ____         #')
print(r'#  / ___|  / \  / ___|_   _/ _ \|  _ \        #')
print(r'# | |     / _ \ \___ \ | || | | | |_) |       #')
print(r'# | |___ / ___ \ ___) || || |_| |  _ <        #')
print(r'#  \____/_/   \_\____/ |_| \___/|_| \_\       #')
print(r'#                                             #')
print(f'#     Welcome to CASTOR v{version}               #')
print(r'#     	Simongini+2024, MNRAS                 #')
print(r'#                                             #')
print(r'#   In case of issues please contact          #')
print(r'#      andrea.simongini@inaf.it               #')
print(r'###############################################')
print('\n\n')


# Interactive plotting 

plt.ion()


''' 
	Filters parameters and constants: 

	limit_sp_points = limit of spectral points to build synthetic spectra with 
	n_templates     = how many templates to build
	alpha_templ     = noise parameter for GP
	n_cores         = how many cores to use 
	neutrino_weight = relative weight of neutrinos 
	cal_parms       = calibration parameters
	bandpass        = dictionary containing filter names, width and effective wavelength 
	lines           = dictionary containing effective wavelength of common lines 
	extco           = dictionary containining the extinction of different filters from York et al. 2004 
	
	Note that every filter has to be in the AB system!
''' 

start_time = TT.time()  
	
			 
n_templates     = 50				
alpha_templ     = 1e-12				 
cal_parms       = ([0.9063, 1.4732 ])	
neutrino_weight = 100 / 0.1   
light_vel_A     = 3e+18         
light_vel_km    = 3e+5          
Hubble_70       = 70   
Mag_sun         = 4.74    
L_sun           = 3.828e33      
M_sun_g         = 1.98892e+33  
R_sun_mpc       = 2.25e-14      
h_planck        = 6.63e-27      
k_boltz         = 1.38e-16      
mpc_to_cm       = 3.0856775813e+24
days_to_sec     = 86400
km_to_cm        = 100000
n_cores         = 5 


bandpass               = {}

filterlist  = ['u', 'g', 'r', 'i', 'z', 'y', 
	           'U', 'B', 'V', 'R', 'I', 'Y',
	           'w2', 'm2', 'w1', 
	           'J', 'H', 'K', 'Ks']

bandpass = { 'u': (3546, 457), 'g': (4767, 928), 'r': (6156, 813), 'i': (7472, 891), 'z': (8917, 1183), 'y': (10100, 600), 
	     'U': (3571, 524), 'B': (4344, 797), 'V': (5456, 807), 'R': (6442, 1377), 'I': (7994, 1079), 'Y': (10380, 2130), 
	     'w2': (1928, 657), 'm2': (2246, 498), 'w1': (2600, 693),
	     'J':(12500, 3110), 'H':(16300, 4020), 'K':(21900, 3930)
	     }

lines = {'Nitrogen'     : 3485 , 'Calcium II'   : 3706 , 'Calcium H&K'  : 3932 , 
         'Hydrogen d'   : 4102 , 'Hydrogen g'   : 4341 , 'Helium Ia'    : 4471 ,
         'Magnesium I]' : 4571 , 'Hydrogen b'   : 4861 , 'Iron IIa'     : 4900 ,
         'Iron IIb'     : 5300 , 'Iron IIc'     : 5500 , 'Oxygen Ia'    : 5577 ,
         'Helium Ib'    : 5876 , 'Silicon II'   : 6150 , 'Oxygen Ib'    : 6300 ,
         'Hydrogen a'   : 6563 , 'Helium Ic'    : 6678 , 'Helium Id'    : 7065 ,
         'Iron IId'     : 7155 , '[Calcium II]' : 7300 , 'Oxygen Ic'    : 7330 , 
         'Carbon IV'    : 7724 , 'Oxygen Id'    : 7776 , 'Calcium nir'  : 8500 
        }

colors = {'u': 'royalblue', 'g': 'springgreen', 'r': 'firebrick', 'i': 'darkorange', 'z': 'brown', 'y': 'pink', 
		  'U': 'mediumblue', 'B': 'deepskyblue', 'V': 'limegreen', 'R': 'red', 'I': 'orange',
		  'm2': 'darkviolet','w1': 'mediumpurple','w2': 'indigo',
          'J': 'red',  'H': 'lightcoral', 'Ks': 'darkred'
		  }


#Extinction coefficients in A_lam / E(B-V).
extco = {'u': 4.786,  'g': 3.587, 'r': 2.471, 'i': 1.798,  'z': 1.403, 'y': 1.228, 
         'U': 4.744,  'B': 4.016, 'V': 3.011, 'R': 2.386, 'I': 1.684,  'Y': 1.228,
         'w2': 8.795, 'm2': 9.270, 'w1': 6.432,
         'J': 0.813, 'H': 0.516,  'Ks':0.337 }

# Solar magnitudes in different filters.
solar_magnitudes = {'u': 6.27, 'g': 5.06, 'r': 4.64, 'i': 4.52, 'z': 4.51, 'y': 4.50,
                    'U': 6.4981, 'B': 5.31, 'V': 4.80, 'R': 4.60, 'I': 4.51, 
                    'J': 4.54, 'H': 4.66, 'Ks': 5.09, 
                    'w2': 10.2868, 'm2': 6.3952, 'w21': 8.4705}

# Some useful functions 

def sampling_step_templ(x):
    return np.diff(np.sort(x))

def sampling_step(x):
    dx = np.diff(np.sort(x))
    return dx[dx != 0]


####

# First step: enter the name of the supernova-of-study 

print('*** Step 1: supernova-of-study ***') 

sn_name = input('\n> Enter SN name:      ') 

# Create the output directory 

localpath   = os.path.dirname(os.path.abspath(__file__))
data_path   = os.path.join(localpath, 'Analysis', sn_name)
outdir      = os.path.join(data_path, "output_" + sn_name) 
if not os.path.exists(outdir): os.makedirs(outdir)

# Find light curves 

lc_path     = os.path.join(data_path, "light_curves", sn_name + '.dat') 

if os.path.exists(lc_path): 
	print(f'\n> light curves file found: {lc_path}') 
else:
	print(f'\n> light curves file not found') 
	print(f'\n *** Stopping ***')
	sys.exit(1) 

	
# Collecting light curves data 
	
def collect_lightcurves(lc_path):
	lc_set    = defaultdict(list)
	band_dict = defaultdict(list)
	with open(lc_path, 'r') as file:
		for line in file:
			parts = line.strip().split()
			time  = float(parts[0])
			mag   = float(parts[1])
			emag  = float(parts[2])
			if emag == 0.0: 
				emag = mag * 0.1 
			band  = parts[3]
			band_dict[band].append((time, mag, emag))
	for band, values in band_dict.items():
		lc_set[f'time_{band}'], lc_set[f'mag_{band}'], lc_set[f'emag_{band}'] = zip(*values)
	lc_filters         = list(band_dict.keys())
	ordered_lc_filters = [band for band in filterlist if band in lc_filters]
	return ordered_lc_filters, lc_set

def time_of_explosion(lc_filters, lc_set):
	t0, t1, dt = [], [], []
	for filtro in lc_filters:
		time = np.array(lc_set['time_' + filtro])
		t1.append(time[0])
		dt.append(abs(np.mean(np.diff(time[:10]))))
	t0        = np.array(t1) - np.array(dt)
	min_value = np.min(t0)
	err_t0    = np.std(dt)
	return min_value, err_t0

def scaling_light_curves(lc_filters, lc_set, t0): 
	for filtro in lc_filters: 
		lc_set['time_' + filtro] = np.array(lc_set['time_' + filtro]) - t0
	return lc_set 

lc_filters, lc_set = collect_lightcurves(lc_path)

# Estimate or enter the time of explosion 
	
q1 = 0 

while q1!= 'y' and q1!= 'n': 

	q1 = input('\n> Do you want to estimate automatically the time of explosion? [Y/n]:      ').strip().lower() 

	if q1 == 'y': 
		t0, err_t0         = time_of_explosion(lc_filters, lc_set)
		print(f'\n> Time of explosion: {round(t0, 2)} +/- {round(err_t0, 2)}')
	if q1 == 'n': 
		t0     = float(input('\n> Enter time of explosion in MJD:      '))
		err_t0 = t0 * 0.05 

lc_set = scaling_light_curves(lc_filters, lc_set, t0)
		
# Interpolation by Gaussian Processes
		

def process_filter(filtro, lc_set, ref_set, is_bol = False):
	warnings.simplefilter("ignore")

	x = np.array(lc_set['time_' + filtro])
	y = np.array(lc_set['mag_' + filtro])
	yerr = np.array(lc_set['emag_' + filtro])

	lengthscale1 = np.min(sampling_step_templ(x))
	lengthscale2 = np.max(sampling_step_templ(x))
	kernel = (C(1.0, (1e-2, 1e2)) * Matern(length_scale=lengthscale1, nu=1.5) +
			C(1.0, (1e-2, 1e2))   * Matern(length_scale=lengthscale2, nu=1.5) +
			C(1.0, (1e-2, 1e2))   * Matern(length_scale=1.0, nu=1.5)) + WhiteKernel(noise_level=1e-1)

	gp = GaussianProcessRegressor(kernel=kernel, alpha=yerr**2, n_restarts_optimizer=10, normalize_y=True)
	gp.fit(x[:, np.newaxis], y)

	ref_time = np.array(ref_set['time_' + filtro])

	if is_bol: 
		ref_filter = max(lc_filters, key=lambda f: lc_set['time_' + f][-1] - lc_set['time_' + f][0])
		ref_time   = lc_set['time_' + ref_filter]

	t       = np.linspace(np.min(ref_time), np.max(ref_time), 250)
	mu, std = gp.predict(t[:, np.newaxis], return_std=True)

	centroid = bandpass[filtro][0] 
	wave     = np.repeat(centroid, len(mu))
	flux     = light_vel_A * 10**(- 0.4* (((mu - cal_parms[1]) / cal_parms[0]) + 48.6)) / (centroid **2)

	return filtro, t, mu, std, wave, flux

def gaussian_process_on_lcs(lc_set, lc_filters, ref_set=lc_set, is_bol = False):
    
    results = Parallel(n_jobs=n_cores)(delayed(process_filter)(f, lc_set, ref_set, is_bol) for f in lc_filters)
    
    gp_set = defaultdict(list)
    for filtro, t, mu, std, wave, flux in results:
        gp_set['time_' + filtro] = t
        gp_set['emag_' + filtro] = std
        gp_set['mag_' + filtro]  = mu
        gp_set['wave_' + filtro] = wave 
        gp_set['flux_' + filtro] = flux 
		
    return gp_set
		
gp_set = gaussian_process_on_lcs(lc_set, lc_filters)

###

# Plot the light curves and the GP on them 

plt.ion()  
plt.figure(figsize=(8, 6))  

for filtro in lc_filters: 

	time, mag, emag = np.array(lc_set['time_' + filtro]), lc_set['mag_' + filtro], lc_set['emag_' + filtro]
	gp_time, gp_mag, gp_err = gp_set['time_' + filtro], gp_set['mag_' + filtro], gp_set['emag_' + filtro]

	plt.errorbar(time, mag, yerr=emag, fmt='.', color=colors[filtro], label=filtro) 
	plt.plot(gp_time, gp_mag, '--', color=colors[filtro], alpha = 0.5) 

plt.xlabel("Days after the explosion")
plt.ylabel("Apparent magnitude")
plt.title('Gaussian Process Interpolated Light Curves')
plt.legend(loc='best')
plt.grid(True)
plt.gca().invert_yaxis()
plt.savefig(outdir + "/light_curves.png", dpi = 300)
plt.show(block=True) 


####

# Second step: compare the light curves with those from the training set 
		
print('\n*** Step 2: comparison with the training set  ***') 
print('\n> This process may require some time ...')


training_path = os.path.join(localpath, 'Training_Set')
training_set  = pd.read_excel(training_path + '/Training_Set.xlsx' )['Name'].tolist()

def chi_squared(gp_set, ref_gps, common_filters):
	chi2_total_norm         = ([])
	chi2_total              = ([])
	for filtro in common_filters:
		mag1                = ref_gps['mag_' + filtro]
		mag2                = np.array(gp_set['mag_' + filtro])
		mag2                = np.sum(mag1)/np.sum(mag2) * mag2
		chi2_single, _      = chisquare(f_obs=mag1, f_exp=mag2)
		chi2_total_norm     = np.append(chi2_total_norm, chi2_single/len(mag2))
		chi2_total          = np.append(chi2_total, chi2_single)
	total_chi2              = np.sum(chi2_total_norm)
	total_chi2_norm         = total_chi2 / len(common_filters)
	return (total_chi2_norm)



total_chi2    = []
for ref_sn in training_set: 
	ref_path              = os.path.join(training_path, 'data_lightcurves', ref_sn + '.dat') 
	ref_filters, ref_lcs  = collect_lightcurves(ref_path)
	common_filters        = set(ref_filters) & set(lc_filters)
	if len(common_filters) <=2 or sn_name == ref_sn:
		chi_norm = np.inf
	else: 
		ref_gps           = gaussian_process_on_lcs(ref_lcs, common_filters, lc_set)
		chi_norm          = chi_squared(gp_set, ref_gps, common_filters)
	total_chi2.append(chi_norm) 

final_name  = training_set[np.argmin(total_chi2)]
final_chi   = total_chi2[np.argmin(total_chi2)]

print(f'\n> The reference supernova is {final_name} with a chi-square of {round(final_chi, 5)}')

###

# Plot the comparison

ref_path              = os.path.join(training_path, 'data_lightcurves', final_name + '.dat') 
ref_filters, ref_lcs  = collect_lightcurves(ref_path)
common_filters        = set(ref_filters) & set(lc_filters)
ref_gps               = gaussian_process_on_lcs(ref_lcs, common_filters, lc_set)

plt.ion() 
num_rows   = (len(common_filters) + 1) // 2
fig, axes  = plt.subplots(num_rows, 2, figsize=(8, 3 * num_rows))
axes       = axes.flatten()

for i, filtro in enumerate(common_filters):
    ax = axes[i]  

    sn_time, sn_mag       = lc_set['time_' + filtro], np.array(lc_set['mag_' + filtro])
    sn_gp_time, sn_gp_mag = gp_set['time_' + filtro], np.array(gp_set['mag_' + filtro])
    ref_time              = ref_lcs['time_' + filtro]
    ref_mag               = np.array(ref_lcs['mag_' + filtro])
    gp_mag                = np.array(ref_gps['mag_' + filtro])
    scale                 = np.min(sn_gp_mag) / np.min(gp_mag)
    ref_mag_scaled        = ref_mag * scale
    ax.scatter(ref_time, ref_mag_scaled, label=final_name, color='red', marker = '1', alpha=0.7)
    ax.plot(sn_gp_time, sn_gp_mag, label=sn_name, color='blue', linewidth=2)
    ax.set_title(f"Filter: {filtro}")
    ax.set_xlabel("Days after the explosion")
    ax.set_ylabel("Magnitude (free-scale)")
    ax.invert_yaxis()  
    ax.grid()
    ax.legend()
for i in range(len(common_filters), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig(outdir + "/comparison.png", dpi = 300)
plt.show(block=True) 


####

# Third step: collect spectral data
		
print('\n*** Step 3: collect spectral data  ***') 
		
q2 = 0 
while q2!= 'y' and q2!='n': 

	q2 = input('\n> Do you have spectral data [y/n]:      ').strip().lower() 

	if q2 == 'y': spectral_availability = 1 
	if q2 == 'n': spectral_availability = 0 

# Collect spectral data 

def collect_spectra(name, sp_path, t0=0):

    sp_set       = defaultdict(list)
    all_files    = os.listdir(sp_path)
    all_epochs   = set()
    min_wave     = [] 
    max_wave     = [] 
    
    for file_path in all_files:
        file_name = os.path.join(sp_path, file_path)
        epoch     = float(file_path.split(name + '_')[1].split('.dat')[0]) - t0
        all_epochs.add(epoch)
        data      = ascii.read(file_name)
        wave = np.array(data['col1'])
        flux = np.array(data['col2'])
        sp_set[f'wave_{epoch}'] = wave
        sp_set[f'flux_{epoch}'] = flux
        min_wave.append(np.min(wave))
        max_wave.append(np.max(wave)) 

    files    = np.array(sorted(all_epochs))
    min_wave = np.mean(min_wave) 
    max_wave = np.mean(max_wave) 
     
    return files, sp_set, min_wave, max_wave 

if spectral_availability: 

	sp_path     = os.path.join(data_path,  "spectra") 
	if os.path.exists(sp_path): 
		print(f'\n> spectral file found: {sp_path}') 
		files, sp_set, min_wave, max_wave  = collect_spectra(sn_name, sp_path, t0)
	else:
		print(f'\n> spectral file not found: applying standard CASTOR')
		spectral_availability = 0 

# If you don't have spectra, we have to use CASTOR in the standard way ;) 

if not spectral_availability: 

	sp_path                           = os.path.join(training_path,  "data_spectra", final_name) 
	files, sp_set, min_wave, max_wave = collect_spectra(final_name, sp_path)


# Let's define the maximum epoch

maximum_epoch = np.max(files) + np.mean(np.diff(files)) 
time_series   = np.linspace(0, maximum_epoch, n_templates)


print(f'\n> There are {len(files)} available spectra, with a maximum epoch of t0+{round(maximum_epoch, 2)}')


####

# Fourth step: collect spectral data
		
print('\n*** Step 4: building synthetic spectra  ***') 

# The first thing to do is to collect ALL flux points from
# real spectra and from interpolated light curves. 
# In order not to exceed on computation time, we set 
# a maximum combined sp_points of 5000. You can change.  

limit_sp_points = 5000	


def total_spectral_points(gp_set, lc_filters, files, sp_set, spectral_availability):

    total_flux, total_wave, total_time = [], [], []  
    
    if not spectral_availability:
        lc_time = np.concatenate([gp_set['time_' + filtro] for filtro in lc_filters])
        lc_flux = np.concatenate([gp_set['flux_' + filtro] for filtro in lc_filters])
        lc_wave = np.concatenate([gp_set['wave_' + filtro] for filtro in lc_filters])
    else:
        lc_flux, lc_wave, lc_time = np.array([]), np.array([]), np.array([])
    
    cut_epoch = 1
    
    while True:
        cut_points = 1  
        while cut_points <= 15:  
            cut_files = files[::cut_epoch]
            all_points = np.concatenate([sp_set[f'flux_{epoch}'] for epoch in cut_files])
            if len(all_points) + len(lc_flux) * (not spectral_availability) <= limit_sp_points:
                break  
            cut_points += 1
            total_points = all_points[::cut_points]
            if len(total_points) + len(lc_flux) * (not spectral_availability) <= limit_sp_points:
                break  
        if cut_points > 15:
            cut_epoch += 1
        else:
            break  
    
    if not spectral_availability:
        lc_order = int(np.log10(np.max(lc_flux)))
        sp_order = int(np.log10(np.max(sp_set['flux_' + str(files[-1])])))
        diff_order = lc_order - sp_order if lc_order != sp_order else 0
    else:
        diff_order = 0
    
    for epoch in list(files)[0::cut_epoch]:
        sp_wave = sp_set['wave_' + str(epoch)][0::cut_points]
        sp_flux = sp_set['flux_' + str(epoch)][0::cut_points] * 10 ** diff_order
        sp_time = np.repeat(epoch, len(sp_flux))
        
        total_flux.extend(sp_flux)  
        total_wave.extend(sp_wave)
        total_time.extend(sp_time)

    total_flux = np.array(total_flux)
    total_wave = np.array(total_wave)
    total_time = np.array(total_time)

    # If spectral data is unavailable, include light curve data
    if not spectral_availability:
        total_flux = np.concatenate((total_flux, lc_flux))
        total_wave = np.concatenate((total_wave, lc_wave))
        total_time = np.concatenate((total_time, lc_time))
    
    return total_wave, total_flux, total_time



tot_wave, tot_flux, tot_time = total_spectral_points(gp_set, lc_filters, files, sp_set, spectral_availability)



def synthetic_spectra(tot_wave, tot_flux, tot_time, time_series):
    
    warnings.simplefilter("ignore")

    x    = tot_time
    y    = tot_flux
    z    = tot_wave

    y_scale   = np.mean(y)
    z_scale   = 70
    x_scale_1 = np.min(sampling_step(x))
    x_scale_2 = np.max(sampling_step(x))

    k1 = C(1.0, constant_value_bounds='fixed') * Matern(length_scale=[x_scale_1, z_scale], nu = 1.5)
    k2 = C(1.0, constant_value_bounds='fixed') * Matern(length_scale=[x_scale_2, z_scale], nu = 1.5)

    kernel = y_scale* (k1 + k2)

    gp = GaussianProcessRegressor(kernel=kernel, alpha= alpha_templ)
    gp.fit(np.column_stack([x.flatten(), z.flatten()]), y)

    pred_set = defaultdict(list)
    for epoch in time_series:
        x_pred = epoch
        z_pred = np.linspace(min_wave, max_wave, len(y))
        tot_pred = np.vstack([np.full_like(z_pred, x_pred), z_pred]).T

        pred_mean, _               = gp.predict(tot_pred, return_std=True)
        pred_set[f'flux_{epoch}']  = pred_mean
        pred_set[f'wave_{epoch}']  = z_pred
    return pred_set

controller = ThreadpoolController()

with controller.limit(limits=n_cores, user_api='blas'):
    pred_set = synthetic_spectra(tot_wave, tot_flux, tot_time, time_series)


###
	
# Plot the result

plt.ion()
plt.figure(figsize=(8, 6))

norm     = mcolors.Normalize(vmin=time_series[0], vmax=time_series[-1])
colormap = cm.viridis  

for epoch in time_series:
    wave, flux = pred_set['wave_' + str(epoch)], pred_set['flux_' + str(epoch)] * 10**15
    plt.plot(wave, flux, color=colormap(norm(epoch)))

plt.xlabel("Wavelength [A]")
plt.ylabel("Normalized Flux")
plt.title('Synthetic Spectra')
plt.grid(True)
plt.savefig(outdir + "/spectra.png", dpi=300)
plt.show(block=True)


####

#  Fourth step: collect spectral data
		
print('\n*** Step 5: estimate parameters  ***') 

def automatic_redshift_procedure(time_series, pred_set, diagnostic=0):
    def find_closest_peak(value, peaks, wave):
        peaks_right = peaks[wave[peaks] > value]
        return peaks_right[0] if peaks_right.size > 0 else None
    
    total_redshift, single_redshift, distances = [], [], []
    
    for epoch in time_series:
        wave, flux = pred_set[f'wave_{epoch}'], pred_set[f'flux_{epoch}']
        peaks, _ = find_peaks(flux, height=0)
        closest_peaks = {line: (wave[closest_peak], wave[closest_peak] - value)
                         for line, value in lines.items()
                         if (closest_peak := find_closest_peak(value, peaks, wave)) is not None}
        
        redshifts = [delta / peak_wave for peak_wave, delta in closest_peaks.values()]
        total_redshift.extend(redshifts)
        single_redshift.extend(redshifts)
        distances.extend(delta for _, delta in closest_peaks.values())
    
    if not total_redshift:
        return None, None
    
    if diagnostic == 0:
        final_redshift = np.nanmedian(total_redshift)
    else:
        median_dist = np.nanmedian(distances)
        filtered_redshifts = np.array([z for i, z in enumerate(single_redshift) if distances[i] < median_dist])
        
        if diagnostic == 1:
            final_redshift = np.nanmedian(filtered_redshifts[filtered_redshifts != 0])
        elif diagnostic == 2:
            sorted_redshifts = np.sort(filtered_redshifts[(filtered_redshifts != 0) & (~np.isnan(filtered_redshifts))])
            final_redshift = next((z for z in sorted_redshifts if int(z * light_vel_km / Hubble_70) > 3), 0)
    
    return final_redshift, final_redshift * 0.1


def redshift_estimate(time_series, pred_set):
    diagnostic = 0
    while True:
        z, err_z = automatic_redshift_procedure(time_series, pred_set, diagnostic)
        if z is None:
            print("\n> Unable to estimate redshift.")
            return None, None
        
        print(f'\n> Estimated Redshift: {round(z, 4)} +/- {round(err_z, 4)}')
        response = input('\n> Accept this redshift? [Y/n]: ').strip().lower()
        if response == 'y' or diagnostic == 2:
            return z, err_z
        diagnostic += 1


def extinction_estimate(lc_filters, gp_set): 

    filter_1 = 'B' if 'B' in lc_filters else 'g' 
    filter_2 = 'V' if 'V' in lc_filters else 'i' 

    time_1, mag_1, emag_1 = gp_set['time_' + filter_1], gp_set['mag_' + filter_1], gp_set['emag_' + filter_1] 
    time_2, mag_2, emag_2 = gp_set['time_' + filter_2], gp_set['mag_' + filter_2], gp_set['emag_' + filter_2] 

    index_2 = np.argmin(mag_2) 
    index_1 = np.argmin(np.abs(time_1 - time_2[index_2]))

    Ebv     = mag_1[index_1] - mag_2[index_2] 
    err_Ebv = np.sqrt(emag_1[index_1]**2 + emag_2[index_2]**2)

    return (Ebv, err_Ebv)

def distance_estimate(z, err_z): 

    dist        = z     * light_vel_km / Hubble_70 
    err_dist    = err_z * light_vel_km / Hubble_70 

    return dist, err_dist
    


while True:
        z_input = input('\n> Please enter the redshift (Enter a number or 0.0 if unknown): ').strip()
        try:
            z = float(z_input)
            if z == 0:
                z, err_z = redshift_estimate(time_series, pred_set)
            else:
                err_z = z * 0.05
            
            print(f'\n> Final Redshift: {round(z, 4)} +/- {round(err_z, 4)}')
            break
        except ValueError:
            print("Invalid input. Please enter a numeric value.")


while True:
    dist = input('\n> Please enter the distance (Enter a number or 0.0 if unknown)     ').strip()
    try:
        dist = float(dist)  
        if not dist: 
            dist, err_dist = distance_estimate(z, err_z)
            print(f'\n> Distance: {round(dist, 2)} +/- {round(err_dist, 2)}')

        else: 
            err_dist = dist * 0.05 
        break  
    except ValueError:
        print("Invalid input. Please enter a numeric value.")

while True:
    Ebv = input('\n> Please enter the extinction (Enter a number or 0.0 if unknown)     ').strip()
    try:
        Ebv = float(Ebv)
        if not Ebv: 
            Ebv, err_Ebv = extinction_estimate(lc_filters, gp_set) 
            print(f'\n> Extinction: {round(Ebv, 2)} +/- {round(err_Ebv, 2)}')
        else: 
            err_Ebv = Ebv * 0.05   

        break  
    except ValueError:
        print("Invalid input. Please enter a numeric value.")


###
        
#   First thing to do is to estimate bolometric luminosity 
    

print(f'\n> Estimating bolometric luminosity...')

L_sun   = 3.828e33

solar_magnitudes = {'u': 6.27, 'g': 5.06, 'r': 4.64, 'i': 4.52, 'z': 4.51, 'y': 4.50,
                    'U': 6.4981, 'B': 5.31, 'V': 4.80, 'R': 4.60, 'I': 4.51, 
                    'J': 4.54, 'H': 4.66, 'Ks': 5.09, 
                    'w2': 10.2868, 'm2': 6.3952, 'w21': 8.4705}


def absolute_magnitudes(Ebv, err_Ebv, dist, err_dist, lc_filters): 

    gp_int = gaussian_process_on_lcs(lc_set, lc_filters, lc_set, True)
      
    dist     = dist     * 10**6 
    err_dist = err_dist * 10**6 
    mu       = 5 * np.log10(dist) - 5
    err_mu   = 5 * err_dist / (dist * np.log(10))
    lc_int   = {}
    for filtro in lc_filters:
        time, mag, emag = gp_int['time_' + filtro], gp_int['mag_' + filtro], gp_int['emag_' + filtro]
        Mag             = mag - mu - Ebv*extco[filtro] 
        Emag            = np.sqrt(emag**2 + err_mu**2 + err_Ebv**2) 
        lc_int['time_' + filtro] = time 
        lc_int['mag_'  + filtro] = Mag 
        lc_int['emag_' + filtro] = Emag 

    return lc_int 

lc_int = absolute_magnitudes(Ebv, err_Ebv, dist, err_dist, lc_filters)

def bolometric_light_curve_from_magnitudes(lc_int):

    pseudo_luminosity = {} 
    ref_time          = lc_int['time_' + lc_filters[0]]
    bol_lum           = np.zeros(len(ref_time))
    bol_err           = np.zeros(len(ref_time))

    for filtro in lc_filters:

        time, mag, emag = lc_int['time_' + filtro], lc_int['mag_' + filtro], lc_int['emag_' + filtro]
        Mag_sun         = solar_magnitudes[filtro]
        L_filtro        = L_sun * 10 ** (-0.4 * (mag - Mag_sun))
        err_L           = 0.4 * L_sun * 10 ** (-0.4 * (mag - Mag_sun)) * np.log(10) * emag

        pseudo_luminosity['lum_'  + filtro]  = L_filtro
        pseudo_luminosity['time_' + filtro]  = time
        pseudo_luminosity['err_'  + filtro]  = err_L

        for i in range(len(bol_lum)): 

            bol_lum[i] += L_filtro[i] 
            bol_err[i] += err_L[i] ** 2 

    bol_err     = np.sqrt(bol_err)
    max_index   = np.argmax(bol_lum)
    lum_bol_max = bol_lum[max_index]
    lum_bol_err = bol_err[max_index]

    return ref_time, bol_lum, bol_err, lum_bol_max, lum_bol_err

###
        
#   Second thing to do is to estimate pseudo-bolometric luminosity 


def bolometric_light_curve_from_sed(dist, err_dist, Ebv, z, lc_set, lc_filters): 
    
    dist_in_cm               = dist     * mpc_to_cm
    err_dist_in_cm           = err_dist * mpc_to_cm
    integrated_flux          = np.array([])
    ref_time                 = np.array([])
    gp_int                   = gaussian_process_on_lcs(lc_set, lc_filters, lc_set, True)
    total_time               = max((gp_int['time_' + filtro] for filtro in lc_filters), key=len)
    

    for t in total_time:     

        total_flux  = [] 
        total_wave  = [] 

        for filtro in lc_filters:   

            time, mag   = gp_int['time_' + filtro], gp_int['mag_' + filtro] 
            eff_wav     = bandpass[filtro][0]
            
            wave        = np.repeat(eff_wav, len(mag))
            flux        = light_vel_A * 10**(- 0.4* (mag   + 48.6 - Ebv*extco[filtro])) / (eff_wav **2)

            closest_time     = min(time, key = lambda x: abs(x-(t)))
            idx              = list(time).index(closest_time)

            if abs(t-closest_time) < 0.01: 

                total_flux.append(flux[idx])
                total_wave = np.append(total_wave, wave[idx])

            total_wave = np.array(total_wave) / (1+z) 

        if len(total_wave) > 1: 

            integrated_flux = np.append(integrated_flux, trapz(total_flux, total_wave))  
            ref_time        = np.append(ref_time, t) 
        else:

            pass
            
    luminosity  = 4 * np.pi * integrated_flux * dist_in_cm**2
    err_lum     = 8 * np.pi * integrated_flux * dist_in_cm * err_dist_in_cm 
    index_max   = np.argmax(luminosity) 
    lum_max     = luminosity[index_max]
    err_lum_max = err_lum[index_max]

    return ref_time, luminosity, err_lum, lum_max, err_lum_max


###

#  Plot luminosity 

bol_time, bol_lum, bol_err, lum_bol_max, lum_bol_err             = bolometric_light_curve_from_magnitudes(lc_int)
pseudo_time, pseudo_lum, pseudo_err, lum_pseudo_max, lum_err_max = bolometric_light_curve_from_sed(dist, err_dist, Ebv, z, gp_set, lc_filters)


plt.figure(figsize=(8, 6))  
plt.fill_between(bol_time, bol_lum - bol_err, bol_lum + bol_err, color='orange', alpha=0.5)
plt.plot(bol_time, bol_lum, '--', color='red', label = 'Bol Lum') 
plt.fill_between(pseudo_time, pseudo_lum - pseudo_err, pseudo_lum + pseudo_err, color='cyan', alpha=0.5)
plt.plot(pseudo_time, pseudo_lum, '--', color='blue', label = 'Pseudo Lum') 
plt.yscale('log')
plt.xlabel("Days after the explosion")
plt.ylabel("Luminosity [erg/s]")
plt.title('Bolometric luminosity curve')
plt.legend(loc='best')
plt.grid(True)
plt.savefig(outdir + "/luminosity.png", dpi=300)
plt.show(block=True) 


###
        
#   Third thing to do is the estimate of velocity 

print(f'\n> Estimating velocity...')


def time_at_maximum(lc_filters, gp_set, err_t0): 
    
    if 'V' in lc_filters:
        filtro = 'V' 
    elif 'i' in lc_filters:
        filtro = 'i' 
    else:
        filtro = 'Ks'

    mag = gp_set['mag_' + filtro]
    max_index = np.argmin(mag) 
    time = gp_set['time_' + filtro]
    err_time = np.mean(np.diff(time))
    tmax = time[max_index]  
    err_tmax = np.sqrt(err_time**2 + err_t0**2) 

    return tmax, err_tmax 

tmax, err_tmax       = time_at_maximum(lc_filters, gp_set, err_t0)

class velocity_class(): 

    def __init__(self, pred_set, time_series): 

        self.pred_set = pred_set 
        self.time_series = time_series 

    def find_closest_peak(self, value, peaks, wave): 

        peaks_right = peaks[wave[peaks] > value]
        if peaks_right.size == 0:
            return None
        closest_peak = peaks_right[0]

        return closest_peak
    
    def find_the_peaks(self, epoch): 

        wave     = np.array(self.pred_set['wave_' + str(epoch)])
        flux     = np.array(self.pred_set['flux_'  + str(epoch)])
        peaks, _ = find_peaks(flux, height=0)  

        closest_peaks = {}
        for line, value in lines.items(): 
            closest_peak = self.find_closest_peak(value, peaks, wave)
            if closest_peak is not None:
                closest_peaks[line] = (wave[closest_peak], wave[closest_peak] - value, value)

        distances = [] 
        for value in closest_peaks.values(): 
            distances.append(value[1]) 


        median_value = np.median(distances)
        saved_lines = {} 
        aux = []

        for line, (peak_wave, distance, value) in closest_peaks.items():
            
                if distance < median_value:
                    if peak_wave not in aux:
                        saved_lines[line] = (peak_wave, distance, value)
                        aux.append(peak_wave)
                    else:
                        for existing_line, (existing_peak_wave, existing_distance, existing_value) in saved_lines.items():
                            if existing_peak_wave == peak_wave:
                                if distance < existing_distance:
                                    saved_lines[line] = (peak_wave, distance, value)
                                    del saved_lines[existing_line]
                                break
        return saved_lines 


    def velocity(self): 

        save_files = []
        err_vel  = []
        doppler    = []
        err_dop    = []
        velocity_lis   = [] 
        diagnostic = False

        wave, flux = self.pred_set['wave_' + str(self.time_series[0])], self.pred_set['flux_' + str(self.time_series[0])]
        saved_lines = self.find_the_peaks(self.time_series[0])
        best_index_1 = None
        best_mean_r_squared = -np.inf  
        index_1_values = [0.1, 0.5, 1, 2, 5]

        for index_1 in index_1_values:
            rr = []  

            for line, (peak_wave, distance, value) in saved_lines.items():
                mask = np.logical_and(wave >= value - index_1 * distance, wave <= peak_wave + 100)
                in_wave, in_flux = wave[mask], flux[mask]
                params = np.polyfit(in_wave, in_flux, deg=5)
                fitted_curve = np.poly1d(params)(in_wave)
                total_sum_squares = np.sum((in_flux - np.mean(in_flux))**2)
                residuals = in_flux - fitted_curve
                sum_squared_residuals = np.sum(residuals**2)
                r_squared = 1 - (sum_squared_residuals / total_sum_squares)
                rr.append(r_squared)

            mean_r_squared = np.mean(rr)

            if mean_r_squared > best_mean_r_squared:
                best_mean_r_squared = mean_r_squared
                best_index_1 = index_1
                
        for rr_limit in [0.9, 0.87, 0.85, 0.83, 0.8, 0.7]:

            for epoch in self.time_series: 
                
                wave, flux = self.pred_set['wave_' + str(epoch)], self.pred_set['flux_' + str(epoch)]
                saved_lines = self.find_the_peaks(epoch)
                for line, (peak_wave, distance, value) in saved_lines.items():
                    mask = np.logical_and(wave >= value - best_index_1 * distance, wave <= peak_wave + 100)
                    in_wave, in_flux = wave[mask], flux[mask]
                    params = np.polyfit(in_wave, in_flux, deg=5)
                    fitted_curve = np.poly1d(params)(in_wave)
                    total_sum_squares = np.sum((in_flux - np.mean(in_flux))**2)
                    residuals = in_flux - fitted_curve
                    sum_squared_residuals = np.sum(residuals**2)
                    r_squared = 1 - (sum_squared_residuals / total_sum_squares)

                    if r_squared > rr_limit: 
                        
                        maxi_wave    = peak_wave
                        reduced_wave = fitted_curve[in_wave < peak_wave]
                        if len(reduced_wave) == 0: continue 
                        mini_wave    = in_wave[np.argmin(reduced_wave)]

                        if maxi_wave > mini_wave: 

                            diagnostic    = True 
                            observed_line = (maxi_wave + mini_wave) / 2 
                            doppler       = np.append(doppler, abs(observed_line - value)/ value)
                            err_w         = np.mean(np.diff(in_wave))
                            err_obs       = 1/np.sqrt(2) * err_w
                            err_dop       = np.append(err_dop, err_obs / value)

                if diagnostic: 

                    final_error    = np.mean(err_dop) * light_vel_km
                    final_velocity = np.mean(doppler)* light_vel_km
                    err_vel.append(final_error)
                    velocity_lis.append(final_velocity)
                    save_files.append(epoch) 

            if len(velocity_lis) >  2: 
                break

        return np.array(velocity_lis), np.array(err_vel), save_files

    def expansion_velocity(self, velocity_lis, err_vel, save_files, tmax): 
            
        save_tmax     = save_files.index(min(save_files, key = lambda x: abs(x-(tmax))))  
        vel_at_max    = velocity_lis[save_tmax]
        err_vel_max   = err_vel[save_tmax] 

        return vel_at_max, err_vel_max
    
velocity_lis, err_vel, save_files = velocity_class(pred_set, time_series).velocity()
vel_at_max, err_vel_max           = velocity_class(pred_set, time_series).expansion_velocity(velocity_lis, err_vel, save_files, tmax)



###

# Plot 


plt.figure(figsize=(8, 6))  
plt.fill_between(save_files, velocity_lis - err_vel, velocity_lis + err_vel, color='cyan', alpha=0.5)
plt.plot(save_files, velocity_lis, '--', color='blue', label = 'Vel') 
plt.xlabel("Days after the explosion")
plt.ylabel("Velocity [km/s]")
plt.title('Evolution of expansion velocity')
plt.legend(loc='best')
plt.grid(True)
plt.savefig(outdir + "/velocity.png", dpi=300)
plt.show(block=True) 


###

# Now the rest of the parameters


def kinetic_energy(luminosity, err_lum, tmax, err_tmax): 

    trise_in_s     = tmax     * days_to_sec 
    err_trise_in_s = err_tmax * days_to_sec

    energy         = luminosity * trise_in_s * neutrino_weight 
    err_energy     = np.sqrt( (trise_in_s * err_lum)**2 + (luminosity * err_trise_in_s) **2 ) * neutrino_weight 

    return energy, err_energy

def mass_of_the_ejecta(vel_at_max, err_vel_max, energy, err_energy): 

    vel_in_cm     = vel_at_max  * km_to_cm  
    err_vel_in_cm = err_vel_max * km_to_cm

    mass_ejecta   = (10/3) * energy / vel_in_cm**2 
    err_mass_ej   = (10/3/vel_in_cm**2) * np.sqrt(err_energy**2 + (2*energy*err_vel_in_cm/vel_in_cm)**2) 

    mass_ejecta   = mass_ejecta / M_sun_g 
    err_mass_ej   = err_mass_ej / M_sun_g 

    return mass_ejecta, err_mass_ej

def progenitor_mass(mass_ejecta, err_mass_ej): 

    mass_ns     = 1.2 
    mass_bh     = 10 
    mass_pr     = (round(mass_ejecta + mass_ns, 2), round(mass_ejecta + mass_bh, 2)) 
    err_mass_pr = err_mass_ej 

    return mass_pr, err_mass_pr 

def mass_of_nichel(ref_time, luminosity, tmax): 

    def detect_linear_decay(x, y, tmax): 

        def r_squared_estimate(x, y):
            slope, intercept, _, _, _ = linregress(x, y)
            y_fit = slope * x + intercept
            ss_res = np.sum((y - y_fit) ** 2)      
            ss_tot = np.sum((y - np.mean(y)) ** 2)  
            r_squared = 1 - (ss_res / ss_tot)
            return r_squared

        starting_point = int(tmax + 20)
        decay_time     = [] 
        decay_lum      = [] 

        for index in range(starting_point, int(np.max(x)) - 2):

            x1 = x[x > index] 
            y1 = y[x > index]
            r_squared  = r_squared_estimate(x1, y1)

            if r_squared > 0.95: 

                decay_time = x1 
                decay_lum  = y1 
                break

        return decay_time, decay_lum

    decay_time, decay_lum = detect_linear_decay(ref_time, luminosity, tmax)

    if len(decay_time) > 0: 

        gamma_1       = 1.32e-6 
        gamma_2       = 1.02e-7
        time_in_s     = decay_time * days_to_sec 
        s             = 3.90e+10 * np.exp(-gamma_1*time_in_s) + 6.78e+9 * (np.exp(-gamma_2*time_in_s) - np.exp(-gamma_1*time_in_s))
        pars          = np.polyfit(s, decay_lum, 1)
        mass_ni       = pars[0] / M_sun_g 
        residuals     = decay_lum - (pars[0] * s + pars[1]) 
        rss           = np.sum(residuals**2)
        std_err       = np.sqrt(rss / (len(s) - 2)) / np.sqrt(np.sum((s - np.mean(s))**2))
        err_mass_ni   = std_err/M_sun_g

        if mass_ni <=0: 
            return 0, 0 
        else: 
            return mass_ni, err_mass_ni

    else:
        return 0, 0 
    


energy, err_energy   = kinetic_energy(lum_pseudo_max, lum_bol_err, tmax, err_tmax) 
mej, err_mej         = mass_of_the_ejecta(vel_at_max, err_vel_max, energy, err_energy)
mass_pr, err_mass_pr = progenitor_mass(mej, err_mej)
mni, err_mni         = mass_of_nichel(pseudo_time, pseudo_lum, tmax)


###
# Photospheric fit 

print(f'\n> Estimating photospheric parameters...')


def blackbody(wave, temperature, radius, zeta, distance=dist):

  black_func_cost = ( 2 * h_planck * light_vel_A**2 ) * 10e+16
  black_func_var  = ( 1 / (np.exp(h_planck * light_vel_A / (wave * k_boltz * temperature)) -1 )) / wave**5
  modified_flux   = ( zeta * radius / distance ) ** 2 * np.pi * black_func_cost * black_func_var * wave

  return modified_flux



def photosphere_fit(lc_filters, time_series, pred_set, dist): 

    temp_min, temp_max     = 30000, 50000
    radius_min, radius_max = 10000, 30000
    temp_steps             = 5
    radius_steps           = 5

    temperature_ranges = np.linspace(temp_min, temp_max, temp_steps)
    radius_ranges      = np.linspace(radius_min, radius_max, radius_steps) 


    epoch_temp_results = [[] for _ in time_series]
    epoch_rad_results  = [[] for _ in time_series]

    for temp_bound in temperature_ranges:
        for radius_bound in radius_ranges:

            bounds_temperature = np.array([1000, temp_bound])
            bounds_radius      = np.array([1000, radius_bound]) * R_sun_mpc
            bounds_zeta        = np.array([0, 1])
            bounds             = ([bounds_temperature[0], bounds_radius[0], bounds_zeta[0]], 
                                  [bounds_temperature[1], bounds_radius[1], bounds_zeta[1]])
            
            initial_guess = [14000, 2000*R_sun_mpc, 0.2]

            iteration_temp = []
            iteration_rad  = []

            for epoch_idx, epoch in enumerate(time_series):               
                wave, flux = pred_set[f'wave_{epoch}'], pred_set[f'flux_{epoch}']
                sed = []
                wave_center = []
                epoch_flux = np.array([])

                for filtro in lc_filters:
                    center   = bandpass[filtro][0]
                    width    = bandpass[filtro][1]
                    low_wave = center - width / 2
                    up_wave  = center + width / 2
                    if np.min(wave) < low_wave and np.max(wave) > up_wave:
                        mask = np.logical_and(wave >= low_wave, wave <= up_wave)
                        in_wave  = wave[mask]
                        in_flux  = flux[mask]
                        int_flux = trapz(in_flux, in_wave)
                        if int_flux != 0: 
                            sed         = np.append(sed, int_flux)
                            wave_center = np.append(wave_center, center)

                fit_func   = partial(blackbody, distance=dist)
                err        = sed * 0.1
                pars, cov  = curve_fit(fit_func, wave_center, sed, p0=initial_guess, sigma=err, maxfev=50000, bounds=bounds)
                variances  = np.diag(cov)
                errors     = np.sqrt(variances)
                ft         = 1 / np.sqrt(pars[2])
                final_temp = pars[0] / ft  
                final_rad  = pars[1] / R_sun_mpc 

                iteration_temp.append(final_temp)
                iteration_rad.append(final_rad)

            # Append this iteration's results for each epoch
            for epoch_idx in range(len(time_series)):
                epoch_temp_results[epoch_idx].append(iteration_temp[epoch_idx])
                epoch_rad_results[epoch_idx].append(iteration_rad[epoch_idx])

    # Average results for each epoch
    average_temp = [np.mean(temp_vals) for temp_vals in epoch_temp_results]
    average_rad  = [np.mean(rad_vals) for rad_vals in epoch_rad_results]

    return average_temp, average_rad

average_temp, average_rad  = photosphere_fit(lc_filters, time_series, pred_set, dist)

###
# Plot 


plt.figure(figsize=(8, 6))
plt.scatter(time_series, average_temp, color='red', marker = '*', label='Temperature (K)')
plt.scatter(time_series, average_rad,  color='green', marker = 'v', label='Radius (Rsun)')
plt.xlabel("Days after the explosion")
plt.ylabel('Photospheric parameters')
plt.grid(True)
plt.tight_layout()
plt.legend(loc='best')
plt.savefig(outdir + "/photosphere.png", dpi=300)
plt.show(block=True) 



### 
# Final results 


print('\n*** Step 6: Print results  ***') 

filepath = os.path.join(outdir, "results.txt")

# Apri il file in modalità scrittura
with open(filepath, "w") as f:
    f.write(f'Results of the analysis of {sn_name}\n\n')
    f.write(f"> Time of explosion [MJD]         : {int(t0)} +/- {int(err_t0)}\n")
    f.write(f"> Time of maximum [MJD]           : {int(tmax+t0)} +/- {int(err_tmax)}\n")
    f.write(f"> Redshift [-]                    : {float(z)} +/- {float(err_z)}\n")
    f.write(f"> Distance [Mpc]                  : {round(dist,2)} +/- {round(err_dist,2)}\n")
    f.write(f"> Extinction [-]                  : {float(Ebv)} +/- {round(err_Ebv,4)}\n")
    f.write(f"> Velocity of the ejecta [km/s]   : {int(vel_at_max)} +/- {int(err_vel_max)}\n")
    f.write(f"> Bolometric maximum [erg/s]      : {lum_bol_max:.2e} +/- {lum_bol_err:.2e}\n")
    f.write(f"> Pseudo-bol maximum [erg/s]      : {lum_pseudo_max:.2e} +/- {lum_err_max:.2e}\n")
    f.write(f"> Total energy [erg]              : {energy:.2e} +/- {err_energy:.2e}\n")
    f.write(f"> Mass of the ejecta [Msun]       : {round(mej,2)} +/- {round(err_mej,2)}\n")
    f.write(f"> Mass of nichel [Msun]           : {round(mni,2)} +/- {round(err_mni,2)}\n")
    f.write(f"> Max temperature [K]             : {int(np.max(average_temp))} +/- {int(np.max(average_temp)*0.05)}\n")
    f.write(f"> Max radius [Rsun]               : {int(np.max(average_rad))} +/- {int(np.max(average_rad)*0.05)}\n")
    f.write(f"> Max Radius of progenitor [Rsun] : {int(average_rad[0])} +/- {int(average_rad[0]*0.05)}\n")
    f.write(f"> Mass of progenitor [Msun]       : {float(mass_pr[0]), float(mass_pr[1])} +/- {round(err_mej,2)}\n")

print(f"Results saved in {filepath}")

end_time = TT.time() 
final_time = end_time - start_time 

minutes, seconds = divmod(final_time, 60)

print("\n*** Analysis completed in {:02.0f} min, {:.2f} sec! ***".format(minutes, seconds))

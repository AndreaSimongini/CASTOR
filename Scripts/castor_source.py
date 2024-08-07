######################################################################################################                   
#                                               ######################################################
#   ____    _    ____ _____ ___  ____           ################################. ####################
#  / ___|  / \  / ___|_   _/ _ \|  _ \          ###############################         ##############
# | |     / _ \ \___ \ | || | | | |_) |         #######################              #     ###########
# | |___ / ___ \ ___) || || |_| |  _ <          ####################                       ########### 
#  \____/_/   \_\____/ |_| \___/|_| \_\         ##################.                       ############ 
#                                               ##################                        ############
#                                               #################                     ################
#                                               ################                      ################
#  Welcome to CASTOR v1.0                       ###############                         .#############
#                                               ###############                     ##   #############
#  In case of problems please                   ###        .-                     +###################      
#  contact                                      #                                  ###################   
#  andrea.simongini@inaf.it                     #              #                    ################## 
#                                               ##          ##########################################                                                 
#                                               ######################################################  
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#
#   Define the needed libraries  
#
#
import  numpy                               as np
import  pandas                              as pd 
import  os 
import  math 
import  time                                as TT
import  george
import  warnings
import  matplotlib.pyplot                   as plt 
from    tabulate                            import tabulate
from    collections                         import defaultdict
from    joblib                              import Parallel, delayed
from    astropy.io                          import ascii
from    scipy                               import optimize   as op
from    scipy.stats                         import chisquare
from    scipy.optimize                      import curve_fit
from    scipy.integrate                     import trapezoid  as trapz
from    functools                           import partial
from    george.kernels                      import Matern32Kernel
from    multiprocessing                     import Pool
from    sklearn.gaussian_process            import GaussianProcessRegressor
from    sklearn.gaussian_process.kernels    import ConstantKernel, Matern
from    sklearn.exceptions                  import ConvergenceWarning
#
#   Some parameters you may want to change  
#
limit_sp_points = 5000
n_templates     = 50
alpha_templ     = 1e-12
n_cores         = 5 
acc             = 0.8 
cal_parms       = m, q = ([0.9186, 1.2233])
#
#
#   Other ancillary definitions 
#
#
bandpass               = {}
bandpass['filterlist'] = filterlist = ["UVW2", "UVM2", "UVW1", "u'", "u", "U", "B", "g", 
                                        "g'", "V", "r", "r'", "R", "i", "i'", "I", "z", "z'", 
                                        "Y", "J", "H", "Ks", "K"]
bandpass['centroid']   = centroid   = ([1928, 2246, 2600, 3500, 3546, 3571, 4344, 4767, 
                                        4825, 5456, 6156, 6261, 6442, 7472, 7672, 7994, 8917, 9097,
                                        10380, 12500, 16300, 21450, 21900])                  
bandpass['fwhm']       = fwhm       = ([657, 498, 693, 80, 457, 524, 797, 928, 1379, 
                                        807, 813, 1382, 1377, 891, 1535, 1079, 1183, 1370, 
                                        2130, 3110, 4020, 3770, 3930])

lines = {'Nitrogen'     : 3485 , 'Calcium II'   : 3706 , 'Calcium H&K'  : 3932 , 
         'Hydrogen d'   : 4102 , 'Hydrogen g'   : 4341 , 'Helium Ia'    : 4471 ,
         'Magnesium I]' : 4571 , 'Hydrogen b'   : 4861 , 'Iron IIa'     : 4900 ,
         'Iron IIb'     : 5300 , 'Iron IIc'     : 5500 , 'Oxygen Ia'    : 5577 ,
         'Helium Ib'    : 5876 , 'Silicon II'   : 6150 , 'Oxygen Ib'    : 6300 ,
         'Hydrogen a'   : 6563 , 'Helium Ic'    : 6678 , 'Helium Id'    : 7065 ,
         'Iron IId'     : 7155 , '[Calcium II]' : 7300 , 'Oxygen Ic'    : 7330 , 
         'Carbon IV'    : 7724 , 'Oxygen Id'    : 7776 , 'Calcium nir'  : 8500 
        }

light_vel_A     = 3e+18         
light_vel_km    = 3e+5          
Hubble_70       = 70            
M_sun_g         = 1.98892e+33  
R_sun_mpc       = 2.25e-14      
h_planck        = 6.63e-27      
k_boltz         = 1.38e-16      
mpc_to_cm       = 3.0856775813e+24
days_to_sec     = 86400
km_to_cm        = 100000
neutrino_weight = 100 / 0.1 


def sampling_step(x):
    dx = np.diff(np.sort(x))
    return dx[dx != 0]

def nll(p, y, gp, computed):
  gp.set_parameter_vector(p)
  ll = gp.log_likelihood(y, quiet=True)
  return -ll if np.isfinite(ll) else 1e25

def grad_nll(p, y, gp, computed):
  gp.set_parameter_vector(p)
  return -gp.grad_log_likelihood(y, quiet=True)

def blackbody(wave, temperature, radius, zeta, distance):

  black_func_cost = ( 2 * h_planck * light_vel_A**2 ) * 10e+16
  black_func_var  = ( 1 / (np.exp(h_planck * light_vel_A / (wave * k_boltz * temperature)) -1 )) / wave**5

  modified_flux   = ( zeta * radius / distance ) ** 2 * np.pi * black_func_cost * black_func_var * wave

  return modified_flux

def read_excel(excel_path): 
	
    name_file = 'Training_Set.xlsx' 
    list_of_names = pd.read_excel(excel_path + name_file)['Name'].tolist()
    list_of_types = pd.read_excel(excel_path + name_file)['Type'].tolist()
    list_of_reds  = pd.read_excel(excel_path + name_file)['Redshift'].tolist()
   
    return list_of_names, list_of_types, list_of_reds
#
#
#
#   
#   Class objects 
#
#
#
#
#

class comparison(): 

    def __init__(self, sn_name, sn_path, training_set, training_path):

        self.sn_name        = sn_name  
        self.sn_path        = sn_path 
        self.training_set   = training_set 
        self.training_path  = training_path


    def collect_lightcurves(self, lc_path):
        lc_set    = defaultdict(list)
        band_dict = defaultdict(list)
        with open(lc_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                time  = float(parts[0])
                mag   = float(parts[1])
                emag  = float(parts[2])
                band  = parts[3]
                band_dict[band].append((time, mag, emag))
        for band, values in band_dict.items():
            lc_set[f'time_{band}'], lc_set[f'mag_{band}'], lc_set[f'emag_{band}'] = zip(*values)
        lc_filters         = list(band_dict.keys())
        ordered_lc_filters = [band for band in filterlist if band in lc_filters]
        return ordered_lc_filters, lc_set

    def time_of_explosion(self, lc_set, lc_filters):
        t0, t1, dt = [], [], []
        for filtro in lc_filters:
            time = np.array(lc_set['time_' + filtro])
            t1.append(time[0])
            dt.append(abs(np.mean(np.diff(time[:10]))))
        t0        = np.array(t1) - np.array(dt)
        min_value = np.min(t0)
        err_t0    = np.std(dt)
        return min_value, err_t0


    def rescale_time(self, lc_filters, lc_set, t0):
        new_set = {}
        for filtro in lc_filters:
            time                        = np.array(lc_set['time_' + filtro])
            mag                         = np.array(lc_set['mag_' + filtro])
            emag                        = np.array(lc_set['emag_' + filtro])
            scale_time                  = time - t0
            sorted_indices              = np.argsort(scale_time)
            scale_time                  = scale_time[sorted_indices]
            mag                         = mag[sorted_indices]
            emag                        = emag[sorted_indices]
            new_set['time_%s' % filtro] = scale_time
            new_set['mag_%s' % filtro]  = mag
            new_set['emag_%s' % filtro] = emag
        return new_set


    def process_filter(self, filtro, ref_new_set, ref_time):
        x               = np.array(ref_new_set['time_' + filtro])
        y               = np.array(ref_new_set['mag_' + filtro])
        yerr            = np.array(ref_new_set['emag_' + filtro])
        amplitude       = np.mean(y)
        lengthscale0    = np.mean(x)
        lengthscale1    = np.min(sampling_step(x))
        lengthscale2    = np.max(sampling_step(x))
        k0              = amplitude * Matern32Kernel(lengthscale0)
        k1              = amplitude * Matern32Kernel(lengthscale1)
        k2              = amplitude * Matern32Kernel(lengthscale2)
        kernel          = k1 + k2 + k0 
        gp              = george.GP(kernel)
        star            = gp.compute(x, yerr)
        p0              = gp.get_parameter_vector()
        results         = op.minimize(nll, p0, args = (y, gp, star), jac=grad_nll, method="L-BFGS-B")
        gp.set_parameter_vector(results.x)
        mu, cov         = gp.predict(y, ref_time)
        std             = np.sqrt(np.diag(cov))
        return {'time_' + filtro: ref_time, 'mag_' + filtro: mu, 'std_' + filtro: std}

    def gaussian_process(self, common_filters, ref_new_set, sn_new_set):
        gp_set = defaultdict(list)
        results = Parallel(n_jobs=n_cores)(delayed(self.process_filter)(filtro, ref_new_set, np.array(sn_new_set['time_' + filtro])) for filtro in common_filters)
        for result in results:
            gp_set.update(result)
        return gp_set

    def chi_squared(self, gp_set, sn_new_set, common_filters):
        chi2_total_norm         = ([])
        chi2_total              = ([])
        for filtro in common_filters:
            mag1                = gp_set['mag_' + filtro]
            mag2                = np.array(sn_new_set['mag_' + filtro])
            mag2                = np.sum(mag1)/np.sum(mag2) * mag2
            chi2_single, _      = chisquare(f_obs=mag1, f_exp=mag2)
            chi2_total_norm     = np.append(chi2_total_norm, chi2_single/len(mag2))
            chi2_total          = np.append(chi2_total, chi2_single)
        total_chi2              = np.sum(chi2_total_norm)
        total_chi2_norm         = total_chi2 / len(common_filters)
        return (total_chi2_norm)
    

    def plot_comparison(self, final_name, sn_filters, sn_new_set, out_path): 

        ref_filters, ref_lc_set = self.collect_lightcurves(self.training_path + "data_lightcurves/" + final_name + ".dat")
        ref_t0, _               = self.time_of_explosion(ref_lc_set, ref_filters)
        ref_new_set             = self.rescale_time(ref_filters, ref_lc_set, ref_t0)
        common_filters          = [filtro for filtro in sn_filters if filtro in ref_filters]
        gp_set                  = self.gaussian_process(common_filters, ref_new_set, sn_new_set)

        num_rows = (len(common_filters) + 1) // 2
        fig, ax  = plt.subplots(num_rows, 2, figsize=(10, 4 * num_rows))

        if num_rows == 1:
            ax = ax.reshape(1, -1)

        for i, filtro in enumerate(common_filters):
            
            row = i // 2
            col = i % 2
            
            ref_time, ref_mag   = gp_set['time_' + filtro], np.array(gp_set['mag_' + filtro])
            sn_time, sn_mag     = sn_new_set['time_' + filtro], np.array(sn_new_set['mag_' + filtro])
            sn_mag              = np.sum(ref_mag)/np.sum(sn_mag) * sn_mag

            ax[row, col].plot(ref_time, ref_mag, '--', color='red', label=final_name)
            ax[row, col].scatter(sn_time, sn_mag, color='blue', label=self.sn_name)

            ax[row, col].invert_yaxis()
            ax[row, col].set_xlim(-1, 100)
            ax[row, col].set_xlabel('Time')
            ax[row, col].set_ylabel('Magnitude')
            ax[row, col].set_title(f'Filter: {filtro}')
            ax[row, col].grid()
            ax[row, col].legend()

        plt.tight_layout()
        figname = out_path + "comparison.png"
        plt.savefig(figname, dpi=300) 
        plt.close()

    

    def total_comparison(self, out_path, to_save=True):

        sn_filters, sn_lc_set   = self.collect_lightcurves(self.sn_path)
        sn_t0, err_t0           = self.time_of_explosion(sn_lc_set, sn_filters)
        sn_new_set              = self.rescale_time(sn_filters, sn_lc_set, sn_t0)
        total_chi2              = []
        total_t0                = []
        for ref_name in self.training_set: 
            data_path               = self.training_path + "data_lightcurves/" + ref_name + ".dat"
            ref_filters, ref_lc_set = self.collect_lightcurves(data_path)
            ref_t0, _               = self.time_of_explosion(ref_lc_set, ref_filters)
            ref_new_set             = self.rescale_time(ref_filters, ref_lc_set, ref_t0)
            common_filters          = set(ref_filters) & set(sn_filters)
            if len(ref_filters) == 0 or self.sn_name == ref_name: 
                chi_norm = np.inf 
            else: 
                gp_set                  = self.gaussian_process(common_filters, ref_new_set, sn_new_set)
                chi_norm                = self.chi_squared(gp_set, sn_new_set, common_filters) 
            total_chi2.append(chi_norm)
            total_t0.append(ref_t0)
        best_guess  = np.argmin(total_chi2)
        ref_t0      = total_t0[best_guess]
        final_name  = self.training_set[best_guess]
        final_chi   = total_chi2[best_guess]

        if to_save:
            self.plot_comparison(final_name, sn_filters, sn_new_set, out_path)

        print('----------------------------')
        print('The best match for %s'%self.sn_name + ' is %s'%final_name + ' with a chi-squared of %2f.'%final_chi)

        return sn_new_set, sn_filters, sn_t0, err_t0, final_name, final_chi, ref_t0
    
#
#
#
#

class spectra_selection(): 

    def __init__(self, final_name, ref_t0, training_path): 

        self.final_name    = final_name 
        self.ref_t0        = ref_t0 
        self.training_path = training_path 

    def read_single_file(self, file_path):

        file_name = os.path.basename(file_path)
        epoch     = float(file_name.split('_')[1].rsplit('.', 1)[0]) 
        epoch     = epoch - self.ref_t0 
        data      = ascii.read(file_path)

        return epoch, data

    def collect_spectra(self):

        sp_set       = defaultdict(list)
        spectra_path = os.path.join(self.training_path, 'data_spectra', self.final_name)
        all_files    = os.listdir(spectra_path)
        valid_files  = [os.path.join(spectra_path, file) for file in all_files if not (file.startswith('.') or file.lower() == 'desktop.ini')]
        all_epochs   = set()

        with Pool() as pool:
            results = pool.map(self.read_single_file, valid_files)

        for epoch, data in results:

            all_epochs.add(epoch)
            sp_set[f'wave_{epoch}']  = data[0][:]
            sp_set[f'flux_{epoch}']  = data[1][:]

        return sorted(all_epochs), sp_set

    def wavelength_coverage(self, ref_files, ref_sp_set): 

        min_wave = [] 
        max_wave = [] 

        for epoch in ref_files: 
            wave = np.array(ref_sp_set['wave_' + str(epoch)]) 
            min_wave.append(np.min(wave)) 
            max_wave.append(np.max(wave)) 
        
        return ([int(np.median(min_wave)), int(np.median(max_wave))])

    def filters_coverage(self, wave_range): 

        spectral_filters = []
        for i in range(len(centroid)): 
            min_wave = centroid[i] - fwhm[i] / 2
            max_wave = centroid[i] + fwhm[i] / 2
            if wave_range[0] <= min_wave and wave_range[1] >= max_wave: 
                spectral_filters.append(filterlist[i])

        return spectral_filters


    def time_coverage(self, ref_files): 

        maximum_epoch = np.max(ref_files) + np.mean(np.diff(ref_files)) 
        
        return maximum_epoch 


    def final_spectra_selection(self): 

        ref_files, ref_sp_set = self.collect_spectra()
        wave_range            = self.wavelength_coverage(ref_files, ref_sp_set)
        spectral_filters      = self.filters_coverage(wave_range)
        maximum_epoch         = self.time_coverage(ref_files) 

        print('----------------------------')
        print('The reference supernova has %i'%len(ref_files), 'spectra')
        print('----------------------------')
        print('The minimum epoch is at %2f'%min(ref_files), 'days after the explosion')
        print('The maximum epoch is at %2f'%max(ref_files), 'days after the explosion')
        print('----------------------------')
        print('The time cut is set at %i'%maximum_epoch, 'days after the explosion')
        print('----------------------------')
        print('The spectra cover the range between ', (wave_range[0]), '-', (wave_range[1]), 'A')
        print('corresponding to filters: ', spectral_filters)

        return ref_files, ref_sp_set, wave_range, spectral_filters, maximum_epoch

#
#
#
#
    

class build_templates(): 

    def __init__(self, sn_name, sn_new_set, sn_filters, sn_t0,  wave_range, 
                 ref_sp_set, ref_files, ref_t0, maximum_epoch): 
                       
        self.sn_name        = sn_name
        self.lc_set         = sn_new_set
        self.lc_filters     = sn_filters 
        self.t0             = sn_t0
        self.ref_files      = ref_files 
        self.ref_sp_set     = ref_sp_set 
        self.ref_t0         = ref_t0 
        self.min_w          = wave_range[0]
        self.max_w          = wave_range[1]
        self.maximum_epoch  = maximum_epoch


    def process_filter(self, filtro):

        x            = np.array(self.lc_set['time_' + filtro])
        y            = np.array(self.lc_set['mag_'  + filtro])
        yerr         = np.array(self.lc_set['emag_' + filtro])
        amplitude    = np.mean(y)
        lengthscale0 = np.mean(x)
        lengthscale1 = np.min(sampling_step(x))
        lengthscale2 = np.max(sampling_step(x))
        k0           = amplitude * Matern32Kernel(lengthscale0)
        k1           = amplitude * Matern32Kernel(lengthscale1)
        k2           = amplitude * Matern32Kernel(lengthscale2)
        kernel       = k1 + k2 + k0
        gp           = george.GP(kernel)
        star         = gp.compute(x, yerr)
        p0           = gp.get_parameter_vector()
        results      = op.minimize(nll, p0, args = (y, gp, star), jac=grad_nll, method="L-BFGS-B")
        gp.set_parameter_vector(results.x)
        t            = np.linspace(np.min(x), np.max(x), len(y))
        mu, cov      = gp.predict(y, t)
        std          = np.sqrt(cov.diagonal())
        eff_wav      = [centroid[i] for i in range(len(filterlist))  if filterlist[i] == filtro]
        wav_c        = np.repeat(eff_wav[0], len(mu))
        new_mu       = (mu - cal_parms[1]) / cal_parms[0]
        flux         = light_vel_A * 10**(- 0.4* (new_mu + 48.6)) / (eff_wav[0] **2)

        return {'time_' + filtro: t, 'flux_' + filtro: flux, 'wave_' + filtro: wav_c, 'std_' + filtro: std, 'mag_' + filtro: new_mu}

    
    def gaussian_process_on_lcs(self):

        gp_set  = defaultdict(list)
        results = Parallel(n_jobs=n_cores)(delayed(self.process_filter)(filtro) for filtro in self.lc_filters)
        for result in results:
            gp_set.update(result)
        return gp_set
    

    def injection(self, gp_set):
        
        total_flux = np.concatenate([gp_set['flux_'+filtro] for filtro in self.lc_filters])
        total_wave = np.concatenate([gp_set['wave_'+filtro] for filtro in self.lc_filters])
        total_time = np.concatenate([gp_set['time_'+filtro] for filtro in self.lc_filters])

        return total_flux, total_wave, total_time 
    
   
    def cut_real_spectra(self, total_flux): 
        
        cut_epoch = 1
        go        = True

        while go:
            cut_files  = list(self.ref_files)[0::cut_epoch]
            all_points = np.concatenate([self.ref_sp_set['flux_' + str(epoch)] for epoch in cut_files])

            if len(all_points)+len(total_flux) > limit_sp_points:
                medium_points = len(all_points) / len(cut_files)
                medium_desire = limit_sp_points / len(cut_files)
                cut_points    = math.ceil(medium_points / medium_desire)
                total_points  = all_points[0::cut_points]

                if len(total_points)+len(total_flux) > limit_sp_points:
                    cut_epoch += 1
                else:
                    go = False
            else:
                go = False


        return cut_epoch, cut_points 
    
    
    def total_spectra(self, total_flux, total_wave, total_time, cut_epoch, cut_points):

        for epoch in list(self.ref_files)[0::cut_epoch]:
            obs_wave      = self.ref_sp_set['wave_'+str(epoch)]
            obs_flux      = self.ref_sp_set['flux_'+str(epoch)]
            obs_time      = np.repeat(epoch, len(obs_flux))
            total_flux    = np.append(total_flux, obs_flux[0::cut_points])
            total_wave    = np.append(total_wave, obs_wave[0::cut_points])
            total_time    = np.append(total_time, obs_time[0::cut_points])
        return total_flux, total_wave, total_time
    
    
    def safe_region(self, total_flux, total_wave, total_time): 

        save_flux = total_flux[(total_wave>self.min_w)&(total_wave<self.max_w)] 
        save_time = total_time[(total_wave>self.min_w)&(total_wave<self.max_w)]  
        save_wave = total_wave[(total_wave>self.min_w)&(total_wave<self.max_w)] 
        return save_flux, save_wave, save_time



    def time_series(self):

        min_value   = 0 
        max_value   = self.maximum_epoch  
        time_series = np.linspace(min_value, max_value, n_templates)
        return time_series 
    

    def bi_gaussian_process(self, total_flux, total_wave, total_time, time_series):

        x    = total_time
        y    = total_flux
        z    = total_wave

        y_scale   = np.mean(y)
        z_scale   = 70
        x_scale_1 = np.min(sampling_step(x))
        x_scale_2 = np.max(sampling_step(x))

        k1 = ConstantKernel(1.0, constant_value_bounds='fixed') * Matern(length_scale=[x_scale_1, z_scale], nu = 1.5)
        k2 = ConstantKernel(1.0, constant_value_bounds='fixed') * Matern(length_scale=[x_scale_2, z_scale], nu = 1.5)

        kernel = y_scale* (k1 + k2)

        gp = GaussianProcessRegressor(kernel=kernel, alpha= alpha_templ)

        gp.fit(np.column_stack([x.flatten(), z.flatten()]), y)

        def process_epoch(epoch):
            x_pred   = epoch
            z_pred   = np.linspace(self.min_w, self.max_w, len(y))
            tot_pred = np.vstack([np.full_like(z_pred, x_pred), z_pred]).T
            pred_mean, pred_std = gp.predict(tot_pred, return_std=True)
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            return {
                'flux_%s' % str(epoch): pred_mean,
                'eflux_%s' % str(epoch): pred_std,
                'wave_%s' % str(epoch): z_pred
            }

        pred_set_list = Parallel(n_jobs=n_cores)(delayed(process_epoch)(epoch) for epoch in time_series)
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        

        pred_set = defaultdict(list)
        for result in pred_set_list:
            pred_set.update(result)

        return pred_set
        

    def total_templates(self, out_path, to_save=True): 

        start_templates_time = TT.time()

        templates_path = out_path + "templates/" 
        if not os.path.exists(templates_path): 
            os.makedirs(templates_path) 
        
        gp_set                                  = self.gaussian_process_on_lcs()
        total_flux, total_wave, total_time      = self.injection(gp_set)
        cut_epoch, cut_points                   = self.cut_real_spectra(total_flux)
        total_flux, total_wave, total_time      = self.total_spectra(total_flux, total_wave, total_time, cut_epoch, cut_points)
        total_flux, total_wave, total_time      = self.safe_region(total_flux, total_wave, total_time)
        time_series				                = self.time_series()
        pred_set				                = self.bi_gaussian_process(total_flux, total_wave, total_time, time_series)

        if to_save: 
            
            plt.figure()
            for epoch in time_series: 
                wave, flux = pred_set['wave_' + str(epoch)], pred_set['flux_' + str(epoch )]
                plt.plot(wave, flux)
                plt.grid()
                plt.title('Synthetic spectra')
                plt.xlabel('Wavelenght [$\AA$]')
                plt.ylabel('Flux [erg/s/cm$^2$/$\AA$]')

                file_name = f"{self.sn_name}_{epoch}.txt"
                file_path = templates_path + file_name
                with open(file_path, 'w') as file:
                    for w, f in zip(wave, flux): 
                        file.write(f"{w}\t{f}\n")
            plt.grid()
            plt.tight_layout()
            plt.savefig(out_path + 'templates.png')
            plt.close()

            plt.figure() 
            for filtro in self.lc_filters:
               time, mag       = self.lc_set['time_' + filtro], self.lc_set['mag_' + filtro]
               gp_time, gp_mag = gp_set['time_' + filtro], gp_set['mag_' + filtro]
               gp_mag = gp_mag * cal_parms[0] + cal_parms[1]
               plt.scatter(time, mag, label=filtro)
               plt.plot(gp_time, gp_mag, '--') 
            plt.legend()
            plt.grid()
            plt.gca().invert_yaxis()
            plt.tight_layout()        
            plt.title('Synthetic light curves')
            plt.xlabel('Days after the explosion')
            plt.ylabel('Apparent magnitude') 
            plt.savefig(out_path + 'gp.png')
            plt.close()

        end_templates_time = TT.time() 
        final_templates_time = end_templates_time - start_templates_time 
        print("Operation took {:.6f} seconds".format(final_templates_time))

        return (gp_set, pred_set, time_series)

#
#
#
#
    
class line_fitting():

  def __init__(self, sn_t0, pred_set, time_series, lines):

    self.t0 = sn_t0 
    self.pred_set = pred_set 
    self.time_series = time_series 
    self.lines = lines 


  def starting_epoch(self):

    YES = False
    i   = 0
    while YES == False:

      epoch = self.time_series[i]
      wave  = np.array(self.pred_set['wave_'  + str(epoch)])
      flux  = np.array(self.pred_set['flux_'  + str(epoch)])
      print('Can you see at least one P-Cygni and one emission line?')
      plt.figure(figsize=(8,6))
      plt.title('Spectrum at +' + str(int(epoch)) + ' days after explosion')
      for value in self.lines.values():
        if value < np.max(wave) and value > np.min(wave):
          plt.plot(wave, flux, '-', color='grey')
          plt.axvline(value, 0, 1, color = 'blue')
      plt.show()
      answer = input().lower()
      if answer == 'yes':
        YES = True
      elif answer == 'no':
        YES = False
        i+=1
      else:
        print("Please enter only 'yes' or 'no'")



    return epoch

  def find_lines(self, epoch):

    wave  = np.array(self.pred_set['wave_'  + str(epoch)])
    flux  = np.array(self.pred_set['flux_'  + str(epoch)])
    print('Select the P-Cygni and the emission lines')
    print('The answer should be comma-space separated i.e. Helium Ia, Hydrogen a')
    for line, value in self.lines.items():
      if value < np.max(wave) and value > np.min(wave):
        plt.figure(figsize=(8,6))
        plt.plot(wave, flux, '-', color = 'grey')
        plt.axvline(value, 0, 1, label=line, color = 'blue')
        plt.legend()
        plt.show()

    OK = False 
    while OK == False: 
      answer          = input().split(', ')
      if answer[0] in lines.keys() and answer[1] in lines.keys():
        p_cygni_line  = answer[0]
        emission_line = answer[1]
        OK  = True
      else: 
        print('Please select an existing line')
        OK = False 

    return p_cygni_line, emission_line
  

  def p_cygni_interval(self, epoch, p_cygni_line):

    wave  = np.array(self.pred_set['wave_'  + str(epoch)])
    flux  = np.array(self.pred_set['flux_'  + str(epoch)])

    red, blu = 100, 100  

    stop = False

    while stop == False:

      for line, value in self.lines.items():
        if line == p_cygni_line:
          blue_limit        = value - blu
          red_limit         = value + red
          mask              = np.logical_and(wave >= blue_limit, wave <= red_limit)
          in_wave, in_flux  = wave[mask], flux[mask]
          pars              = np.poly1d(np.polyfit(in_wave, in_flux, 8))
          computed_flux     = pars(in_wave)
          print('Please insert values to enlarge or reduce the interval i.e. "10, -35"')
          print('The first value will add on the left, while the second on the right')
          print('If the interval is okay, just write 0, 0')
          plt.figure(figsize=(8,6))
          plt.plot(in_wave, computed_flux, '--', color = 'red', label = 'fitted flux')
          plt.scatter(in_wave, in_flux, color = 'blue', label = 'obs. flux')
          plt.axvline(value, 0, 1, color = 'green', label = 'line')
          plt.show()

          answer     = input()
          values     = answer.split(', ')
          blue_value = int(values[0])
          red_value  = int(values[1])

          red += red_value
          blu += blue_value
          if red_value == blue_value == 0.0:
            stop = True
          else:
            stop = False



    return blu, red

  def emission_interval(self, epoch, emission_line):

    wave  = np.array(self.pred_set['wave_'  + str(epoch)])
    flux  = np.array(self.pred_set['flux_'  + str(epoch)])

    red_e, blu_e = 20, 50  #this are just hypotetical values to start with

    stop = False

    while stop == False:

      for line, value in self.lines.items():
        if line == emission_line:
          blue_limit        = value - blu_e
          red_limit         = value + red_e
          mask              = np.logical_and(wave >= blue_limit, wave <= red_limit)
          in_wave, in_flux  = wave[mask], flux[mask]

          pars              = np.poly1d(np.polyfit(in_wave, in_flux, 8))
          computed_flux     = pars(in_wave)

          print('Please insert values to enlarge or reduce the interval i.e. "10, -35"')
          print('The first value will add on the left, while the second on the right')
          print('If the interval is okay, just write 0, 0')

          plt.figure(figsize=(8,6))
          plt.plot(in_wave, computed_flux, '--', color = 'red', label = 'fitted flux')
          plt.scatter(in_wave, in_flux, color = 'blue', label = 'obs. flux')
          plt.axvline(value, 0, 1, color = 'green', label = 'line')
          plt.show()

          answer     = input()
          values     = answer.split(', ')
          blue_value = int(values[0])
          red_value  = int(values[1])

          red_e += red_value
          blu_e += blue_value
          if red_value == blue_value == 0.0:
            stop = True
          else:
            stop = False
    return blu_e, red_e


  def p_cygni_fitter(self, p_cygni_line, blu, red):

    save_files  = []
    doppler     = ([])
    err_dop     = ([])

    for epoch in self.time_series:
      wave    = np.array(self.pred_set['wave_'  + str(epoch)])
      flux    = np.array(self.pred_set['flux_'  + str(epoch)])

      value   = self.lines[p_cygni_line]
      mask    = np.logical_and(wave >= value - blu, wave <= value + red)
      in_wave = wave[mask]
      in_flux = flux[mask] 

      params                = np.polyfit(in_wave, in_flux, deg=5)
      fitted_curve          = np.poly1d(params)(in_wave)
      total_sum_squares     = np.sum((in_flux - np.mean(in_flux))**2) 
      residuals             = in_flux - fitted_curve
      sum_squared_residuals = np.sum(residuals**2)
      r_squared             = 1 - (sum_squared_residuals / total_sum_squares)

      if r_squared > 0.9: 
      
        maxi_wave = in_wave[np.argmax(fitted_curve)]
        mini_wave = in_wave[np.argmin(fitted_curve)]

        if maxi_wave > mini_wave: 
                  
          observed_line = (maxi_wave + mini_wave) / 2 
          doppler       = np.append(doppler, abs(observed_line - value)/ value)
          err_w         = np.mean(np.diff(in_wave))
          err_obs       = 1/np.sqrt(2) * err_w
          err_dop       = np.append(err_dop, err_obs / value)
          save_files.append(epoch) 
              
            
    return (save_files, doppler, err_dop)


  def emission_fitter(self, epoch, emission_line, blu_e, red_e):

    wave                  = np.array(self.pred_set['wave_'  + str(epoch)])
    flux                  = np.array(self.pred_set['flux_'  + str(epoch)])
    value                 = self.lines[emission_line]
    mask                  = np.logical_and(wave >= value - blu_e, wave <= value + red_e)
    in_wave               = wave[mask]
    in_flux               = flux[mask]
    params                = np.polyfit(in_wave, in_flux, deg=5)
    fitted_curve          = np.poly1d(params)(in_wave)

    maxi_wave     = in_wave[np.argmax(fitted_curve)]
    observed_line = (maxi_wave)
    if observed_line > value: 
      redshift = (observed_line-value) / value
      err_red  = np.mean(np.diff(in_wave)) / value
        
    return redshift, err_red 


  def classification(self):


    cmap = plt.get_cmap('viridis')

    for i in range(len(self.time_series)): 
        epoch = self.time_series[i] 
        color = cmap(i/len(self.time_series))
        wave, flux = self.pred_set['wave_' + str(epoch)], self.pred_set['flux_' + str(epoch)]
        new_flux = flux * 10**(14)
        plt.plot(wave, new_flux,  color = color)
    for line in lines: 
      plt.axvline(lines[line], 0, 1, color='grey') 
      if 'Helium' in line: 
        plt.axvline(lines[line], 0, 1, color='blue') 
      elif 'Hydrogen' in line: 
        plt.axvline(lines[line], 0, 1, color='red') 

    plt.xlabel('Wavelength [$\AA$]', fontweight = 'demibold', fontsize = 14)
    plt.ylabel('Normalized flux', fontweight = 'demibold', fontsize = 14)
    cax = plt.axes([0.91, 0.11, 0.03, 0.77])
    cmap = 'viridis' 
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=np.max(self.time_series)))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax, ticks=np.linspace(0, np.max(self.time_series), 5))
    cbar.set_label('Days after the explosion', fontsize=14, fontweight='bold')
    cax.invert_yaxis()
    
    plt.show() 
    print('Which class is this?')
    print('Choose between II, IIP, IIn, IIb, Ib, Ib/c, Ic')
    print('\n')
    print('-------------------------------------------------')
    YES = False
    while YES == False:
      answer = input()
      if answer in ['II', 'Ib', 'Ic', 'IIP', 'IIn', 'IIb', 'Ib/c']: 
        YES = True 
      else: 
        YES = False
        print("Please select only a real class")
    classe = answer 
    return classe


  def doppler_and_redshift(self):

    epoch                          = self.starting_epoch()
    p_cygni_line, emission_line    = self.find_lines(epoch)
    blu, red                       = self.p_cygni_interval(epoch, p_cygni_line)
    blu_e, red_e                   = self.emission_interval(epoch, emission_line)
    save_files, doppler, err_dop   = self.p_cygni_fitter(p_cygni_line, blu, red)
    redshift, err_red              = self.emission_fitter(epoch, emission_line, blu_e, red_e)
    classe                         = self.classification()

    return save_files, doppler, err_dop, redshift, err_red, classe
  
#
#
#
#
  
class parameter_estimation():

    def __init__(self, sn_t0, err_t0, sn_new_set, sn_filters, gp_set, spectral_filters, time_series, pred_set): 

        self.t0                 = sn_t0 
        self.err_t0             = err_t0 
        self.lc_filters         = sn_filters
        self.lc_set             = sn_new_set
        self.gp_set             = gp_set 
        self.spectral_filters   = spectral_filters
        self.time_series        = time_series 
        self.pred_set           = pred_set 


    def time_at_maximum(self): 
    
        if 'V' in self.lc_filters:
            filtro = 'V' 
        elif 'i' in self.lc_filters:
            filtro = 'i' 
        else:
            filtro = 'Ks'

        mag = self.gp_set['mag_' + filtro]
        max_index = np.argmin(mag) 
        time = self.gp_set['time_' + filtro]
        err_time = np.mean(np.diff(time))
        tmax = time[max_index] + self.t0 
        err_tmax = np.sqrt(err_time**2 + self.err_t0**2) 

        return tmax, err_tmax 
    
    def time_of_rising(self, tmax, err_tmax): 

        trise     = tmax - self.t0 
        err_trise = np.sqrt(err_tmax**2 + self.err_t0**2)

        return (trise, err_trise)
        
    def absorption(self, trise): 
       
        filter_B = 'B' if 'B' in self.lc_filters else 'g'
        filter_V = 'V' if 'V' in self.lc_filters else 'i'

        mag_V, emag_V = self.gp_set['mag_' + filter_V], self.gp_set['std_' + filter_V] 
        mag_B, emag_B = self.gp_set['mag_' + filter_B], self.gp_set['std_' + filter_B] 
        time_B        = self.gp_set['time_' + filter_B ]

        max_V_index = np.argmin(mag_V) 
        max_V       = mag_V[max_V_index]
        err_max_V   = emag_V[max_V_index]

        max_B_index = 0

        for i in range(len(time_B)):
            if abs(time_B[i] - trise) < abs(time_B[max_B_index] - trise): 
                max_B_index = i 

        max_B     = mag_B[max_B_index]
        err_max_B = emag_B[max_B_index]

        extinction  = (max_B-max_V)
        err_ext     = np.sqrt(err_max_B**2 + err_max_V**2)

        absorption  = 3.1 * extinction
        err_abs     = 3.1 * err_ext

        return absorption, err_abs


    def expansion_velocity(self, doppler, err_dop, save_files, tmax): 
            
        velocity_lis  = doppler * light_vel_km 
        err_vel       = err_dop * light_vel_km 

        save_tmax     = save_files.index(min(save_files, key = lambda x: abs(x-(tmax-self.t0)))) 
        vel_at_max    = velocity_lis[save_tmax]
        err_vel_max   = err_vel[save_tmax] 

        return velocity_lis, err_vel, vel_at_max, err_vel_max
    

    def distance(self, redshift, err_red):

        dist        = redshift * light_vel_km / Hubble_70 
        err_dist    = err_red * light_vel_km / Hubble_70 

        return dist, err_dist
    

    def bolometric_luminosity(self, trise, dist, err_dist): 

        total_flux, total_wave = ([]), ([])

        for filtro in self.lc_filters: 

            flux, wave, time = self.gp_set['flux_' + filtro], self.gp_set['wave_' + filtro], self.gp_set['time_' + filtro] 
            closest_time     = min(time, key = lambda x: abs(x-(trise)))

            if closest_time - trise < 10: 

                flux_at_max  = flux[time==closest_time]
                wave_at_max  = wave[time==closest_time]
                total_flux   = np.append(total_flux, flux_at_max) 
                total_wave   = np.append(total_wave, wave_at_max) 

        int_flux       = trapz(total_flux, total_wave) 
        dist_in_cm     = dist      * mpc_to_cm 
        err_dist_in_cm = err_dist  * mpc_to_cm
        luminosity     = 4 * np.pi * dist_in_cm**2 * int_flux 
        err_lum        = 8 * np.pi * dist_in_cm * err_dist_in_cm * int_flux

        return luminosity, err_lum 
    

    def kinetic_energy(self, luminosity, err_lum, trise, err_trise): 

        trise_in_s     = trise     * days_to_sec 
        err_trise_in_s = err_trise * days_to_sec

        energy         = luminosity * trise_in_s * neutrino_weight 
        err_energy     = np.sqrt( (trise_in_s * err_lum)**2 + (luminosity * err_trise_in_s) **2 ) * neutrino_weight 

        return energy, err_energy
    
    def mass_of_the_ejecta(self, vel_at_max, err_vel_max, energy, err_energy): 

        vel_in_cm     = vel_at_max  * km_to_cm  
        err_vel_in_cm = err_vel_max * km_to_cm

        mass_ejecta   = (10/3) * energy / vel_in_cm**2 
        err_mass_ej   = (10/3/vel_in_cm**2) * np.sqrt(err_energy**2 + (2*energy*err_vel_in_cm/vel_in_cm)**2) 

        mass_ejecta   = mass_ejecta / M_sun_g 
        err_mass_ej   = err_mass_ej / M_sun_g 

        return mass_ejecta, err_mass_ej
    
    def fitting_lcs(self, trise): 

        fit_set = defaultdict(int)
        upper_time = []

        for filtro in self.lc_filters: 
            time = self.lc_set['time_' + filtro] 
            upper_time.append(np.max(time)) 

        decay_time = np.arange(trise+20, np.mean(upper_time), dtype=int, step = 2) 


        for filtro in self.lc_filters: 

            x            = np.array(self.lc_set['time_' + filtro])
            y            = np.array(self.lc_set['mag_'  + filtro])
            yerr         = np.array(self.lc_set['emag_' + filtro])

            amplitude    = np.mean(y)
            lengthscale0 = np.mean(x)
            lengthscale1 = np.min(sampling_step(x))
            lengthscale2 = np.max(sampling_step(x))


            k0           = amplitude * Matern32Kernel(lengthscale0)
            k1           = amplitude * Matern32Kernel(lengthscale1)
            k2           = amplitude * Matern32Kernel(lengthscale2)
            kernel       = k1 + k2 + k0

            gp          = george.GP(kernel)
            star        = gp.compute(x, yerr)

            p0          = gp.get_parameter_vector()
            results     = op.minimize(nll, p0, args = (y, gp, star), jac=grad_nll, method="L-BFGS-B")
            gp.set_parameter_vector(results.x)
            t           = decay_time
            mu, _     = gp.predict(y, t)

            eff_wav     = [centroid[i] for i in range(len(filterlist))  if filterlist[i] == filtro]
            wav_c       = np.repeat(eff_wav[0], len(mu))
            new_mu      = (mu - cal_parms[1]) / cal_parms[0]
            flux        = light_vel_A * 10**(- 0.4* (new_mu + 48.6)) / (eff_wav[0] **2)
            fit_set['time_' + filtro]  = decay_time 
            fit_set['mag_'  + filtro]  = new_mu 
            fit_set['flux_' + filtro]  = flux  
            fit_set['wave_' + filtro]  = wav_c  

        return fit_set, decay_time 

    def mass_of_nichel(self, dist, trise):
        
        fit_set, decay_time = self.fitting_lcs(trise)
        dist_cm             = dist * mpc_to_cm
        gamma_1             = 1.32e-6 
        gamma_2             = 1.02e-7
        lum_array           = []

        for i in range(len(decay_time)):
            flux = ([])
            wave = ([])
            for filtro in self.lc_filters:

                flux = np.append(flux, fit_set['flux_' + filtro][i]) 
                wave = np.append(wave, fit_set['wave_' + filtro][i]) 

            int_flux   = trapz(flux, wave)                             
            Lum        = 4*np.pi * dist_cm**2 * int_flux              
            lum_array  = np.append(lum_array, Lum)

        time_in_s = decay_time * 86400 
        s         = 3.90e+10 * np.exp(-gamma_1*time_in_s) + 6.78e+9 * (np.exp(-gamma_2*time_in_s) - np.exp(-gamma_1*time_in_s))
        pars      = np.polyfit(s,lum_array,1)

        mass_nikel = pars[0] / M_sun_g 

        residuals     = lum_array - (pars[0] * s + pars[1]) 
        rss           = np.sum(residuals**2)
        std_err       = np.sqrt(rss / (len(s) - 2)) / np.sqrt(np.sum((s - np.mean(s))**2))
        err_mass_ni   = std_err/M_sun_g

        return mass_nikel, err_mass_ni
    

    def photosphere_fit(self, redshift, dist, trise): 
            
        black_body_set     = defaultdict(int) 
        bounds_temperature = np.array([1000, 30000])
        bounds_radius      = np.array([1000, 20000]) * R_sun_mpc
        bounds_zeta        = np.array([0, 1])
        bounds             = ([bounds_temperature[0], bounds_radius[0], bounds_zeta[0]], 
                              [bounds_temperature[1], bounds_radius[1], bounds_zeta[1]])
        closest_time       = min(self.time_series, key = lambda x: abs(x-trise)) 

        wave_masks  = []
        wave_center = ([])

        for sp_filtro in self.spectral_filters: 

            wave          = self.pred_set['wave_' + str(self.time_series[0])]
            filter_index  = filterlist.index(sp_filtro)
            low_wave      = centroid[filter_index] - fwhm[filter_index] / 2 
            up_wave       = centroid[filter_index] + fwhm[filter_index] / 2 
            mask          = np.logical_and(wave >= low_wave, wave<= up_wave)    
            wave_center   = np.append(wave_center, centroid[filter_index]/(1+redshift))
            wave_masks.append(mask)


        for epoch in self.time_series:               

            wave, flux  = self.pred_set['wave_' + str(epoch)], self.pred_set['flux_' + str(epoch)]
            sed         = ([])

            for mask in wave_masks: 

                masked_wave = (wave[mask])
                masked_flux = (flux[mask])
                sed         = np.append(sed, trapz(masked_flux, masked_wave))
                err         = sed * 0.1

            fit_func  = partial(blackbody, distance=dist)
            pars, cov = curve_fit(fit_func, wave_center, sed, sigma = err, maxfev=50000, bounds=bounds)
            variances = np.diag(cov)
            errors    = np.sqrt(variances)
            ft         = 1 / np.sqrt(pars[2])
            final_temp = pars[0] / ft  
            err_temp   = np.sqrt( (errors[0]/ft)**2 + (0.5 * pars[0] * errors[2] * ft)**2 )  
            final_rad  = pars[1]   / R_sun_mpc 
            err_rad    = errors[1] / R_sun_mpc 

            if np.any(np.isinf(cov)) == False:
                black_body_set['temperature_'+ str(epoch)] = final_temp
                black_body_set['radius_'     + str(epoch)] = final_rad
                black_body_set['zeta_'       + str(epoch)] = pars[2]
                black_body_set['err_temp_'   + str(epoch)] = err_temp
                black_body_set['err_rad_'    + str(epoch)] = err_rad
                black_body_set['err_zeta_'   + str(epoch)] = errors[2]

            if epoch == closest_time: 
                phot_at_max = ([final_temp, err_temp], [final_rad, err_rad], [pars[2], errors[2]])

            if epoch == 0.0: 
                prog_radius = ([final_rad, err_rad])

        return black_body_set, phot_at_max, prog_radius
    

    def progenitor_mass(self, mass_ejecta, err_mass_ej): 

        mass_ns     = 1.2 
        mass_bh     = 10 
        mass_pr     = (round(mass_ejecta + mass_ns, 2), round(mass_ejecta + mass_bh, 2)) 
        err_mass_pr = err_mass_ej 

        return mass_pr, err_mass_pr 
        

    def total_estimation(self, out_path): 

        tmax, err_tmax                              = self.time_at_maximum() 
        trise, err_trise                            = self.time_of_rising(tmax, err_tmax)
        Av, err_Av                                  = self.absorption(trise)

        save_files, doppler, err_dop, redshift, err_red, classe = line_fitting(self.t0, 
                                                                               self.pred_set, 
                                                                               self.time_series, 
                                                                               lines).doppler_and_redshift()
        velocity_lis, err_vel, vel_at_max, err_vel_max = self.expansion_velocity(doppler, err_dop, save_files, tmax)
        dist, err_dist                              = self.distance(redshift, err_red)
        luminosity, err_lum                         = self.bolometric_luminosity(trise, dist, err_dist)
        energy, err_energy                          = self.kinetic_energy(luminosity, err_lum, trise, err_trise)
        mass_ejecta, err_mass_ej                    = self.mass_of_the_ejecta(vel_at_max, err_vel_max, energy, err_energy)
        mass_nikel, err_mass_ni                     = self.mass_of_nichel(dist, trise)
        black_body_set, phot_at_max, prog_radius    = self.photosphere_fit(redshift, dist, trise)
        mass_pr, err_mass_pr                        = self.progenitor_mass(mass_ejecta, err_mass_ej)


        luminosity = luminosity / 10**(41) 
        err_lum    = err_lum / 10**(41) 
        energy, err_energy  = energy / 10**(51) , err_energy / 10**(51)
        data = [
            ['Hubble constant (assumed)', 'km/s/Mpc', Hubble_70],
            ['Class','-',  classe, '-'],
            ['Time of explosion',  'Mjd', round(self.t0, 2), round(self.err_t0, 2)],
            ['Time of maximum luminosity', 'Mjd', round(tmax, 3), round(err_tmax, 2)],
            ['Time of rising',  'day', round(trise, 2), round(err_trise, 2)],
            ['Absorption',  'mag', round(Av, 2), round(err_Av, 2)],
            ['Distance',  'Mpc', round(dist, 2), round(err_dist, 2)],
            ['Redshift', '-', round(redshift, 4), round(err_red, 4)],
            ['Shock velocity', 'km/s', round(float(velocity_lis[0]), 2), round(float(err_vel[0]), 2)], 
            ['Velocity at tmax', 'km/s', round(vel_at_max, 2), round(err_vel_max, 2)],
            ['Bolometric luminosity',  'L41', round(luminosity, 2), round(err_lum, 2)],
            ['Kinetic energy',  'E51', round(energy, 2), round(err_energy, 2)],
            ['Mass of the ejecta', 'Msun', round(mass_ejecta, 4), round(err_mass_ej, 4)],
            ["Mass of Nickel", 'Msun', round(mass_nikel, 3), round(err_mass_ni, 3)],
            ['Photospheric temperature at tmax',  'K', round(phot_at_max[0][0], 2), round(phot_at_max[0][1], 2)],
            ['Photospheric radius at tmax', 'Rsun', round(phot_at_max[1][0], 2), round(phot_at_max[1][1], 2)],
            ['Dilution factor at tmax', '-', round(phot_at_max[2][0], 2), round(phot_at_max[2][1], 2)],
            ["Progenitor's radius", 'Rsun', round(prog_radius[0], 2), round(prog_radius[1], 2)],
            ["Progenitor's mass",  "Msun", mass_pr, round(err_mass_pr, 3)]
        ]

        table = tabulate(data, headers=['Parameter', 'Units', 'Value', 'Error'], tablefmt='simple')
        print(table)            

        with open(out_path + 'results.txt', 'w') as f:
            f.write(table)    

        # Saving figures 

        fig_name = out_path + 'velocity.png' 
        plt.figure() 
        plt.scatter(save_files, velocity_lis, color = 'blue', label = 'Points') 
        plt.axvline(trise, 0, 1, ls = '--', color='red', label = 'tmax') 
        plt.xlabel('Days after the explosion')
        plt.ylabel('Velocity [km/s]') 
        plt.legend()
        plt.grid()
        plt.savefig(fig_name) 
        plt.close() 	


#
#
#
#
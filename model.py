import math
import numpy as np

c_ = 2.997e8          # [m/s]
h_ = 6.62606957E-34   # [m2 kg / s]
k_ = 1.3806488E-23    # [m2 kg s-2 K-1]

def wiens_temp(lam_peak):
    temperature = 2.898e-3/lam_peak # m*K/m
    
    return temperature   # [K]

def planck_function(wave, temperature):     

    Blam_A = 2.*h_*(c_**2)*(wave**-5 )
    
    Blam_B =  (math.e)**(h_*c_/(wave*k_*temperature))-1.0
    
    Blam = Blam_A / Blam_B
    

    return Blam      # [ W*sr-1*m-3]

## The funciton below allows beta to be fixed, or varied:
## Whether or not the photometric band coverage we have is sufficent enough to contrain beta...
## that may depend on who you ask.
## Another approach could be to generate a beta-temperature grid for each pixel

def modified_blackbody(wave, temperature, tau=1e-6, beta=1.50, wave0=550e-6):  
    ## Default beta set to 1.5, following Planck 2015 Results
   
    Blam = planck_function(wave, temperature) # # [ W*sr-1*m*3]
    
    Bnu = Blam * (wave / c_)

    unit_conversion_factor = 1e26 # [Jy] / [W/m^2/Hz]
    
    emissivity_factor = (wave0/wave)**beta # wave0 is the reference wavelength. 
 
    S = tau*Bnu*emissivity_factor*unit_conversion_factor             #tau is really "tau_wave0". The dust optical depth @ wave0
    
    return S   # [Jy]

def modBBfit(wavelengths, initial_guesses, region, color_corrected=False):
    
    if color_corrected == False:
        fluxes                = fd_all[region,7:] # The IR Photometry results (column 0 is the AME phot.)
        sigma                 = fd_err_all[region,7:]
    else:
        fluxes                = fd_all_cc[region,7:] # The IR Photometry results (column 0 is the AME phot.)
        sigma                 = fd_err_all[region,7:]
    
    
    popt, pcov = scipy.optimize.curve_fit( \
        modified_blackbody, # The modified blackbody function
        wavelengths, # The wavelength values for the IR photometry
        fluxes, # The IR photometry results
        p0=initial_guesses, # Initial guesses for Temperature, Beta (and optical depth?)
        sigma=sigma, # The IR photometry errors
        absolute_sigma=True, 
        check_finite=True)
    
    return popt, pcov
     
if __name__ == 'main':
    
    ## Read-in the data and set
    
    fd_all        = np.genfromtxt("photometry.dat", delimiter=',', comments='#')
    fd_all_err    = np.genfromtxt("errors.dat", delimiter=',', comments='#')
    band_centers  = [90e-6, 100e-6, 140e-6, 160e-6, 350e-6, 550e-6] # [m] Only bands > 70 microns for this modified blackbody fitting.
    nregions      = len(fd_all)
    nbands        = len(band_centers)
    
    

    tau_init = 0.000001           # Intial guess for dust optical depth at 550 microns. (the longest-wavelength band in our observational data set)
    
    T_all           = np.ones(nregions) # Intialize arrays to store the parameters and covariances
    tau_all         = np.ones(nregions)
    T_all_cov       = np.ones(nregions)
    tau_all_cov     = np.ones(nregions)

    for region in range(0,nregion):
        
        ## Setup the fitting for each region:
        
        fir_flux        = fd_all[7:,region]              # Take only the fluxes for the bands longer than 70 microns
        fir_max         = np.max(fir_flux)               # Find the highest flux among the far-IR bands
        fir_max_pos     = np.where(fir_flux == fir_max)  # Get the column number of the hightest-flux far-IR band
        lam_peak        = band_centers[fir_max_pos]      # Get the center wavelength corresponding the hightest-flux band
        T_init          = wiens_temp(lam_peak)           # Get the T of a blackbody peaking at the highest-flux band's center wavelength
        initial_guesses = (T_init, tau_init)             # Set the T from above as the intial guess for the modBB fitting. tau_init is given manually, above
        
        ####### The actual fitting  ###########
        #######################################
        
        par, cov = modBBfit(band_centers, initial_guesses, region)  ## Note: by default, beta is fixed at 1.5. tau is calculated at 550 microns.
    
        T_all[region]      = par[0]
        tau_all[region]    = par[1]
        T_all_cov[region]  = cov[0]
        tau_all_cov        = cov[1]

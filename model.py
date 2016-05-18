import math

from astropy.units import cds
cds.enable()  

c_ = 2.997e8          # [m/s]
h_ = 6.62606957E-34   # [m2 kg / s]
k_ = 1.3806488E-23    # [m2 kg s-2 K-1]


def planck_function(wave, temperature):     

    Blam_A = 2.*h_*(c_**2)*(wave**-5 )
    
    Blam_B =  (math.e)**(h_*c_/(wave*k_*temperature))-1.0
    
    Blam = Blam_A / Blam_B
    

    return Blam      # [ W·sr−1·m−3]

## The funciton below allows beta to be fixed, or varied:
## Whether or not the photometric band coverage we have is sufficent enough to contrain beta...
## that may depend on who you ask.
## Another approach could be to generate a beta-temperature grid for each pixel

def modified_blackbody(wave, temperature, tau=1e-6, beta=1.50, wave0=550e-6):  
    ## Default beta set to 1.5, following Planck 2015 Results
   
    Blam = planck_function(wave, temperature) # # [ W·sr−1·m−3]
    
    Bnu = Blam * (wave / c_)

    unit_conversion_factor = 1e26 # [Jy] / [W/m^2/Hz]
    
    emissivity_factor = (wave0/wave)**beta # wave0 is the reference wavelength. 
                                           #The wave0 you choose determines at which wavelength the optical depth is given   
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
     

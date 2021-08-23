library(cmdstanr)
library(tidyverse)

options(mc.cores = parallel::detectCores())

d = read.csv('../data/h3.csv')

d = d %>% 
    filter(
        FLAG == 0 & 
        SNR > 3 &
        Vrot < 5 & 
        Teff < 7000 &
        V_tan < 1000 &
        abs(V_gsr) < 400 &
        R_gal > 50 &
        !is.na(GAIAEDR3_RA) &
        !is.na(GAIAEDR3_PMRA) &
        !is.na(GAIAEDR3_PARALLAX) &
        !is.na(Vrad) &
        !is.na(X_gal) 
    )

nrow(d)

# are variables observed?
posobs = !is.na(d$GAIAEDR3_RA) %>% as.integer
pmobs = !is.na(d$GAIAEDR3_PMRA) %>% as.integer 
distobs = !is.na(d$dist_adpt) %>% as.integer
vlosobs = !is.na(d$Vrad) %>% as.integer 
d[is.na(d)] = 0

# Format the data to be passed in to Stan. 
# MULTIPLY kpc by 1000 to get pc
# DIVIDE mas by 1000 to get arcsec

mas_to_deg = 2.777777777777778e-07

standata = list(
    N = nrow(d),
    pos_obs = posobs,
    dist_obs = distobs,
    pm_obs = pmobs,
    vlos_obs = vlosobs,
    
    ra_measured = d$GAIAEDR3_RA, 
    dec_measured = d$GAIAEDR3_DEC,
    dist_measured = d$dist_adpt, # kpc
    
    ra_err = d$GAIAEDR3_RA_ERROR * mas_to_deg, 
    dec_err = d$GAIAEDR3_DEC_ERROR * mas_to_deg,
    dist_err = d$dist_adpt_err, # kpc
    
    Xgc = d$X_gal,
    Ygc = d$Y_gal,
    Zgc = d$Z_gal,
    
    pmra_measured = d$GAIAEDR3_PMRA, # mas/yr
    pmdec_measured = d$GAIAEDR3_PMDEC, # mas/yr
    vlos_measured = d$Vrad,
    
    pmra_err = d$GAIAEDR3_PMRA_ERROR, # mas/yr
    pmdec_err = d$GAIAEDR3_PMDEC_ERROR, # mas/yr
    vlos_err = d$Vrad_err,
    
    pos_corr = d$GAIAEDR3_RA_DEC_CORR,
    pm_corr = d$GAIAEDR3_PMRA_PMDEC_CORR,
    ra_pmra_corr = d$GAIAEDR3_RA_PMRA_CORR,
    dec_pmdec_corr = d$GAIAEDR3_DEC_PMDEC_CORR,
    ra_pmdec_corr = d$GAIAEDR3_RA_PMDEC_CORR,
    dec_pmra_corr = d$GAIAEDR3_DEC_PMRA_CORR,
    
    grainsize = 1,
    
    # try different priors!
    alpha_mean = 4.,
    alpha_sigma = 0.1,
    beta_mean = 0.3,
    beta_sigma = 0.05
)

# Define a function which generates a dictionary with the initial values. Initial values must be specified for either all parameters xor none of the parameters. The type for each parameter must match what is defined in the Stan model.

initfun = function() {
    list(
        
        p_gamma = 0.4,
        p_phi0 = 65,
        p_beta = 0.5,
        p_alpha = 4.,
        
        pmra = standata$pmra_measured / 5,
        pmdec = standata$pmdec_measured / 5,
        dist = standata$dist_measured,
        vlos = standata$vlos_measured / 5
    )
}

mod = cmdstan_model("../models/h3p.stan")

fit = mod$sample(
    data = standata,
    chains = 2,
    init = initfun,
    iter_warmup = 1000,
    iter_sampling = 1000,
    # output_dir = '../saved/'
    # max_treedepth = 13,
    # adapt_delta = 0.95,
)

fit = stan(file = "../models/gc.stan", # or use another model (.stan file) in the same directory
           data=stan_data,
           # warmup=4e3,
           # iter=4e3 + 2e3,
           chains = 2, # use this for diagnostics
           # chains = parallel::detectCores(),  # use all cores
           init = initfun,
           control=list(
               # max_treedepth=13, 
               # adapt_delta=0.95
           ),
           verbose = TRUE
)

# Run diagnostics.
fit$cmdstan_diagnose()

# Refer to `gme.r` for examples of what to do next.
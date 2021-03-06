library(cmdstanr)
library(tidyverse)
library(latex2exp)
library(posterior)
library(bayesplot)
library(reshape2)
library(ggridges)

# Detect the number of usable cores. Also make sure that Stan reuses models that have not changed so that we don't have to recompile every time.

options(mc.cores = parallel::detectCores())

# Define a function to read in the appropriate data set and do some preprocessing.

import_data = function(dataset) {
    # for globular cluster data set
    if (dataset == 'gc') {
        d = read.csv('../data/GC_Vasiliev.csv') %>%
            filter(!is.na(Vlos)) %>% # require Vlos to exist
            filter(Rgc > 15) # only use GCs outside 15kpc radius)
        # formatting GC names
        d$Fullname = d$Fullname %>% trimws
    }

    # for dwarf galaxy data set
    else if (dataset == 'dg') {
        # ignore these DGs, they are possibly associated with the LMC
        
        d = read.csv('../data/DG_Fritz.csv') %>% 
            filter(!is.na(Vlos)) %>% # require Vlos
            filter(Rgc > 15) # only use DGs outside 15kpc radius
        
        # formatting DG names
        d$Fullname = d$Fullname %>% trimws
        
        remove = c("Carina II", "Carina III", "Horologium I", "Hydrus I")
        d = d %>% filter(!Fullname %in% remove) # filter out DGs we don't want to include
    }
    
    return(d)
}

# read data
d = import_data('gc')

# Format the data to be passed in to Stan. 
stan_data = list(
    N = nrow(d),
 
    ra = d$RA,
    dec = d$DEC,
    plx = d$plx,
    
    Xgc = d$Xgc,
    Ygc = d$Ygc,
    Zgc = d$Zgc,
    
    pmra_measured = d$PMra,
    pmdec_measured = d$PMdec,
    vlos_measured = d$Vlos,
    r_measured = d$Rgc,
    
    pmra_err = d$ePMra,
    pmdec_err = d$ePMdec,
    vlos_err = d$eVlos,
    r_err = d$eRgc
)

# Define a function which generates a dictionary with the initial values. Initial values must be specified for either all parameters xor none of the parameters. The type for each parameter must match what is defined in the Stan model.

initfun = function() {
    list(
        
        # for the GC dataset, these values seem to work well
        p_gamma = runif(1, 0.4, 0.55),
        p_phi0 = runif(1, 60, 80),
        p_beta = runif(1, -0.1, 0.5),
        p_alpha = runif(1, 3.1, 3.5),
        
        # for the DG dataset, the (hyper)priors are different so 
        # use these (random) initial values instead
        # p_beta = runif(1, -2.5, 0.5),
        # p_alpha = runif(1, 3.3, 3.7),
        
        pmra = rnorm(stan_data$N, 0, 1e-6),
        pmdec = rnorm(stan_data$N, 0, 1e-6),
        vlos = rnorm(stan_data$N, 0, 1),
        r = runif(stan_data$N, 0, max(stan_data$r_measured))
    )
}

# Given a Stan file, compile the model and run it, passing in at the minimum the data to use. Here we also specify the number of chains (where each chain is run on a separate thread), the initial values, and some additional parameters controlling the algorithm. 

#`adapt_delta` effectively controls the step size of the numerical integrator. Run the code with lower values (default is 0.8) first, and if the number of divergences is significant enough to impact the validity of the samples then increase this value. This increases the model robustness but may lead to slower sampling. Do not increase this if you only have a few divergences out of a few thousand draws, or if your `n_eff` is already reasonably large.

#`max_treedepth` increases the depth of the binary tree created for each NUTS iteration when the sampler integrates forward and backward in time. If the code complains about saturated tree depths, this value may need to be increased. This will increase the runtime. Note that saturated tree depths are an efficiency problem, not a validity problem. 

#`chains` controls the number of Markov chains to use. The default is 4. Each chain runs on a separate core, and you can set the number of chains to as many cores as you have. However, since `n_eff` of a few hundred is usually enough for inference, more chains are not always necessary. Increasing the number of chains may also increase the runtime. Use `chains=1` for diagnostic purposes. 

#`iter_warmup` controls the number of iterations to perform in the warmup stage. This helps NUTS find the optimal number of steps and the step size, and the samples drawn here are not used for inference purposes. `iter_sampling` controls the total number of iterations to perform. The number of actual usable draws will be `iter - warmup`.  

mod = cmdstan_model("../models/gc.stan")

fit = mod$sample(
    data = stan_data,
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

# Draw samples from the posterior. 
samples = fit$draws() %>% as_draws_df
posterior = fit$draws() %>% as.array
np = nuts_params(fit)

# Here are the parameters of most interest. 
params = c("p_phi0", "p_gamma", "p_alpha", "p_beta")

# Make some plots of the posterior samples.
mcmc_dens_overlay(posterior, pars = params)

# Parallel plot. Check where divergences occur. 
mcmc_parcoord(posterior, pars = params, transformations = scale, np = np)

# Pairs plot
mcmc_pairs(posterior, pars = params, np = np, off_diag_fun = 'hex')

# Show a traceplot of the chains.  These are like usual traceplots, and should look like Gaussian noise with a stable means and (non-zero) variances.
mcmc_trace(posterior, pars = params)

# ridgeline plot of distances, proper motions, los velocities
plot_posterior_ridgelines = function(samples, varname) {
    select(samples, contains(paste(varname, "[", sep = ""))) %>% melt %>% ggplot() + geom_density_ridges(aes(x=value, y=variable), quantile_lines = TRUE, quantiles = 2) 
}
plot_posterior_ridgelines(samples, "pmra")
plot_posterior_ridgelines(samples, "pmdec")
plot_posterior_ridgelines(samples, "r")
plot_posterior_ridgelines(samples, "vlos")


# Define a function which, given the radius in kpc, returns either a large number of mass estimates at that radius, or some summary statistics for the mass at that radius. The resulting mass is in units of $10^{12}$ M$_\odot$.
mass_at_radius = function(r, data, full=FALSE) {
    m = data$p_gamma * data$p_phi0 * 2.325e-3 * (r)^(1 - data$p_gamma) * 1e12
    if (full) {
        return(m)
    } else {
        return (quantile(m, c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975)))
    }
}

# Plot the mass distribution out to some radius, along with various credible intervals.
plot_masses = function(massfunc, data, add=FALSE, color="blue", upto=200) {
    radii = seq(1, upto, length.out = 1000)
    masses = sapply(radii, massfunc, data = data) / 1e12
    par(mfrow=c(1,1))
    if (add) {
        lines(radii, masses['50%',], type='l', ylab=TeX("Mass ($10^{12} M_{sol}$)"), xlab="Radius (kpc)")
    } else {
        plot(radii, masses['50%',], type='l', ylab=TeX("Mass ($10^{12} M_{sol}$)$"), xlab="Radius (kpc)", ylim=c(0, 1.5))
    }
    polygon(c(radii, rev(radii)), c(masses['25%',], rev(masses['75%',])), col=adjustcolor(color,alpha.f=0.4) , fillOddEven = TRUE)
    polygon(c(radii, rev(radii)), c(masses['12.5%',], rev(masses['87.5%',])), col=adjustcolor(color,alpha.f=0.3) , fillOddEven = TRUE)
    polygon(c(radii, rev(radii)), c(masses['2.5%',], rev(masses['97.5%',])), col=adjustcolor(color,alpha.f=0.2) , fillOddEven = TRUE)
}

# Plot the mass distribution out to 300 kpc.
plot_masses(mass_at_radius, samples, color='#5E81AC', upto=300)


# Plot the density of the mass within 200 kpc.
mass_at_radius(200, samples, full=TRUE) %>% 
    data.frame %>% 
    rename(mass='.') %>% 
    ggplot(aes(x=mass)) + 
    geom_density(col='#5E81AC', size=1) +
    geom_vline(xintercept=mass_at_radius(200, samples)['50%'], col='#5E81AC', size=1) +
    xlim(5e11, 2e12)
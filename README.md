# gmestan

This repository provides examples for using Stan to create a hierarchical Bayesian model to estimate the mass of the Milky Way using kinematic data from tracers.

## GME model

This repository contains a Stan port of [GME](https://github.com/gweneadie/GME). There are two models, [`models/gc.stan`](models/gc.stan) and [`models/dg.stan`](models/dg.stan), which can be used to reproduce the results of [Eadie and Jurić 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...875..159E/abstract) and Slizewski et. al (submitted to ApJ) respectively. These models assume independence between the two proper motions (RA and DEC), do not incorporate covariances between phase-space parameters, and the results are produced with complete data only.

Two accompanying datasets are provided:

- [`GC_Vasiliev.csv`](data/GC_Vasiliev.csv): a globular cluster dataset from [Vasiliev 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2832V/abstract)
- [`DG_Fritz.csv`](data/DG_Fritz.csv): a dwarf galaxy dataset from [Fritz 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A.103F/abstract)

## Extended model

The extended model used for analysis of H3 halo stars is provided in [`models/h3.stan`](models/h3.stan). Bring your own data!

## Python

To install [CmdStanPy](https://github.com/stan-dev/cmdstanpy) (a lightweight Python interface to Stan), see the instructions [here](https://cmdstanpy.readthedocs.io/en/v0.9.76/installation.html).

The Python file, [`examples/gme_py.ipynb`](examples/gme_py.ipynb), is a Jupyter notebook. Apart from Pystan and Jupyter, the necessary dependencies to run all the code in the notebook can be installed through `pip` with the following line:

```
pip3 install numpy pandas matplotlib arviz seaborn
```

## R

To install [RStan](https://github.com/stan-dev/rstan) (the R interface to Stan), see the instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

The R file is [`examples/gme.r`](examples/gme.r), which I recommend opening in [RStudio](https://rstudio.com/products/rstudio/). Apart from RStan, the necessary depenencies to run all the code can be installed with the following line:

```
install.packages(c("ggplot2", "dplyr", "latex2exp", "bayesplot", "shinystan"), dependencies=TRUE)
```

## Useful Resources

- [Stan User Guide](https://mc-stan.org/docs/2_25/stan-users-guide/index.html)
- [Stan Reference Manual](https://mc-stan.org/docs/2_25/reference-manual/index.html)
- [Brief Guide to Stan's Warnings](https://mc-stan.org/misc/warnings.html)
- [NUTS Paper](https://arxiv.org/abs/1111.4246)
- [A Conceptual Introduction to Hamiltonian Monte Carlo](https://arxiv.org/abs/1701.02434)

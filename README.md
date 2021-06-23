# gmestan

This repository provides examples for using Stan to create a multilevel model to estimate the mass of the Milky Way using kinematic data from tracers.

Two datasets are provided:

- [`GC_Vasiliev.csv`](data/GC_Vasiliev.csv): a globular cluster dataset from [Vasiliev 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2832V/abstract)
- [`DG_Fritz.csv`](data/DG_Fritz.csv): a dwarf galaxy dataset from [Fritz 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A.103F/abstract)

The Stan code is in [`models`](models). There are several models available, with each variant have two separate files to deal with the GC dataset and the DG dataset. In order from the simplest to the most complex:

- Two models which do not incorporate measurement uncertainties are given in [`gc_nouncertainty.stan`](models/gc_nouncertainty.stan) and [`dg_nouncertainty.stan`](models/dg_nouncertainty.stan).
- The models matching [GME](https://github.com/gweneadie/GME), used to reproduce the results of [Eadie and Jurić 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...875..159E/abstract) and Slizewski et. al in prep, are in [`gc.stan`](models/gc.stan) and [`dg.stan`](models/dg.stan) respectively. These models assume independence between the two proper motions (RA and DEC), and the results are produced with complete data only.
- Two models which incorporate covariance matrices for the proper motions are given in [`gc_covmat.stan`](models/gc_covmat.stan) and [`dg_covmat.stan`](models/dg_covmat.stan). This reduces down to the previous model if the correlations given are zero.
- Two models which deal with missing data (where there are no proper motions and/or line-of-sight velocities available) and incorporate covariance matrices for the available proper motions are given in [`gc_covmat_incomplete.stan`](models/gc_covmat_incomplete.stan) and [`dg_covmat_incomplete.stan`](models/dg_covmat_incomplete.stan). This reduces down to the previous model if all objects have complete data.

## Python

To install [PyStan](https://github.com/stan-dev/pystan) (the Python interface to Stan), see the instructions [here](https://github.com/stan-dev/pystan/blob/develop/doc/installation_beginner.rst).

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

# gmestan

This repository provides examples for using Stan to create a hierarchical Bayesian model to estimate the mass of the Milky Way using kinematic data from tracers. See the associated paper describing the model and code here: [Shen et al. 2021](https://arxiv.org/abs/2111.09327)

## GME model

This repository contains a Stan port of [GME](https://github.com/gweneadie/GME). There are two models, [`models/gc.stan`](models/gc.stan) and [`models/dg.stan`](models/dg.stan), which can be used to reproduce the results of [Eadie and Jurić 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...875..159E/abstract) and Slizewski et. al (submitted to ApJ) respectively. These models assume independence between the two proper motions (RA and DEC), do not incorporate covariances between phase-space parameters, and the results are produced with complete data only.

Example files for how to use the models are in [`examples/gme_py.ipynb`](examples/gme_py.ipynb) for Python and [`examples/gme.r`](examples/gme.r) for R. 

Two accompanying datasets are provided:

- [`GC_Vasiliev.csv`](data/GC_Vasiliev.csv): a globular cluster dataset from [Vasiliev 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2832V/abstract)
- [`DG_Fritz.csv`](data/DG_Fritz.csv): a dwarf galaxy dataset from [Fritz 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A.103F/abstract)

## Extended model

The extended model used for analysis of H3 halo stars is provided in [`models/h3p.stan`](models/h3p.stan). The model can deal with missing positions and/or velocities, and incorporates within-source covariances.

There are again two examples, [`examples/h3.ipynb`](examples/h3.ipynb) and [`examples/h3.r`](examples/h3.r). Bring your own data!

Also provided is a version of the model where RA and Dec. are fixed. This model can be found in [`models/h3p-fixed.stan`](models/h3p-fixed.stan). This model is recommended when uncertainties in position are small (e.g., from Gaia). 

## Python

To install [CmdStanPy](https://github.com/stan-dev/cmdstanpy) (a lightweight Python interface to Stan), see the instructions [here](https://cmdstanpy.readthedocs.io/en/v0.9.76/installation.html).

The Python examples are Jupyter notebooks. Apart from Pystan and Jupyter, the necessary dependencies to run all the code in the notebook can be installed through `pip` with the following line:

```bash
pip3 install numpy pandas matplotlib arviz seaborn
```

## R

To install [CmdStanR](https://github.com/stan-dev/cmdstanr) (the R interface to CmdStan), see the instructions [here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

I recommend opening the provided R files in [RStudio](https://rstudio.com/products/rstudio/). Apart from CmdStanR, the necessary depenencies to run all the code can be installed with the following line:

```r
install.packages(c("tidyverse", "latex2exp", "bayesplot", "shinystan", "ggridges", "reshape2", "posterior"), dependencies=TRUE)
```

## Citation

If you make use of the code in this repository and find it useful, please consider citing [Shen et al. 2021](https://arxiv.org/abs/2111.09327). The BibTeX entry for this paper is: 
```bibtex
@article{Shen_2021,
  title={The Mass of the Milky Way from the H3 Survey},
  author={Shen, Jeff and Eadie, Gwendolyn M. and Murray, Norman and Zaritsky, Dennis and Speagle, Joshua S. and Ting, Yuan-Sen and Conroy, Charlie and Cargile, Phillip A. and Johnson, Benjamin D. and Naidu, Rohan P. and Han, Jiwon Jesse},
  year={2021},
  eprint={2111.09327v1},
  archivePrefix={arXiv},
  primaryClass={astro-ph.GA}
}
```

## Useful Resources

- [Stan User Guide](https://mc-stan.org/docs/2_25/stan-users-guide/index.html)
- [Stan Reference Manual](https://mc-stan.org/docs/2_25/reference-manual/index.html)
- [Brief Guide to Stan's Warnings](https://mc-stan.org/misc/warnings.html)
- [NUTS Paper](https://arxiv.org/abs/1111.4246)
- [A Conceptual Introduction to Hamiltonian Monte Carlo](https://arxiv.org/abs/1701.02434)

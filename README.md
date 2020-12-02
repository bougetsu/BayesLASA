# BayesLASA

We propsed a framework called Bayesian Landmark-based Shape Analysis (BayesLASA) for landmark detection in polygonal chain data.


#### Landmark identification

MCMC_shape

Usage:

```{r}
MCMC_shape(dat,  iter = n,  estK = est.K, gamma_i = gamma_i, alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)

#####
# argument
#
# dat: polygonal chain
# iter: interations for MCMC, default = 100*n, where n is total number of points in the chain.
# estK: estimated number of landmark points, default = 3
# gamma_i: initial gamma, if not speficifed, will generate a random one
# alpha_sigma: default = 3
# beta_sigma: default = 1/n for normalized chain
# ppm_store: if return the ppm matrix
#####
```

#### Case study

Distance- and model-based features were extracted and fit into Cox model.

# Bayesian Landmark-based Shape Analysis (BayesLASA)

## Introduction
Bayesian Landmark-based Shape Analysis (BayesLASA) is a framework for landmark detection in polygonal chain data. Given a polygonal chain, BayesLASA will identify and extract landmark points using Bayesian inference.

## Directory structure

* `code`:
  * `landmark_detection/MCMC_shape.cpp`: The BayesLASA function
  * `landmark_detection/sim_polygon_gaussian.R`: Functions for generating simulated polygons
  * `toolbox/functions.R`: Auxiliary functions for related analyses
* `demo`: Demo dataset for BayesLASA
* `data`: The four datasets to which BayesLASA was applied in the BayesLASA manuscript
* `manuscript_reproducibility`: Files to reproduce the analysis in the BayesLASA manuscript

## Usage
Below, we demonstrate the usage of BayesLASA for landmark identification on the MPEG-7 dataset. The core landmark detection functionality is performed by the `MCMC_shape` function in `code/landmark_detection/MCMC_shape.cpp`.

```{r}
####################
# Main function:
#
# MCMC_shape(dat, 
#           iter = 10000, 
#           estK = 4, 
#           gamma_i = c(1, 0, 0, 1, 0, 1, 0, 0, 1), 
#           updateGp = c(0.8, 0.1, 0.1), 
#           alpha_sigma = 3, 
#           beta_sigma = 0.01, 
#           ppm_store = TRUE,
#           open = FALSE)
####################


####################
# Required arguments:
#
# dat: matrix - a polygonal chain
# iter: numerical scalar - number of MCMC iterations (default = 100 * n, where n is total number of points in the chain)
# estK: integer - estimated number of landmark points (default = 3)
# gamma_i: integer vector - initial gamma (default = random)
# updateGp: numerical vector - probability vector for updating gamma (default = c(0.8, 0.1, 0.1), following the pattern c("add-or-delete", "swap",  "shift") where the vector sum must be 1; to fix K, set the first probability value to 0)
# alpha_sigma: numerical scalar (default = 3)
# beta_sigma: numerical scalar (default = 1/n for normalized chain)
# ppm_store: boolean - returns PPM matrix if TRUE
# open: boolean - polygonal chain is open if TRUE (default is FALSE)
####################

####################
# Function output:
#
# res: List - iter: total number of MCMC iterations,
#             burn: number of iterations for burn-in, 
#             gamma_map: MAP gamma value,
#             gamma_map_index: position index of MAP gamma, 
#             Llist: list of Ls in each iteration,
#             hastings: difference in posterior values in each iteration,
#             posteriors: posterior values in each iteration, 
#             Ks: K in each iteration,
#             accept_r_ad: acceptance rate of new proposed gamma by add-or-delete,
#             accept_r_swap: acceptance rate of new proposed gamma by swap,
#             accept_r_shift: acceptance rate of new proposed gamma by shift,
#             larger_lklh: if new MCMC iteration has larger posterior)
####################
```

#### Case study

The [MPEG-7 dataset](http://www.dabi.temple.edu/∼shape/MPEG7/dataset.html) is a well-known benchmark dataset used in the development of computer vision techniques. Here, we demonstrate BayesLASA's landmark detection using [a version of MPEG-7 that has been converted to polygonal chains](https://github.com/jd-strait/ALDUQ) by Strait et al. for their work ["Automatic Detection and Uncertainty Quantification of Landmarks on Elastic Curves"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6781625/). 

```{r}
set.seed(9080)

## Load packages
library(mcclust) 
library(Rcpp)
library(doParallel)
library(foreach)
registerDoParallel(4)
library(R.matlab)

## Load functions
sourceCpp("code/landmark_detection/MCMC_shape.cpp")
source("code/toolbox/functions.R")

set.seed(9080)

## Load packages
library(mcclust) 
library(Rcpp)
library(doParallel)
library(foreach)
registerDoParallel(4)
library(R.matlab)

## Load functions
sourceCpp("code/landmark_detection/MCMC_shape.cpp")
source("code/toolbox/functions.R")

## Read deer dataset
f = "demo/MPEG7closed.mat"
mdat <- readMat(f)$C.cl
k = 461 # Number of vertices in deer shape
pc = cbind(mdat[1,,k], mdat[2,,k])

## Preprocess data and scale
temp = pc_normalizor(pc)
dat = temp$pc;
dat = dat[-nrow(dat),]
n = nrow(dat)

## Set hyperparameters and algorithm arguments
fold = 100
est.K = round(n/100)
if(est.K <3){est.K = 3}
beta_sigma <- 0.001 

## 4 MCMC chains
res = foreach(i = 1:4) %dopar%{
  # Generate gamma_0
  gamma_i = generate_gamma(n, est.K)
  # Run MCMC in parallel
  MCMC_shape(dat,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
             alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
}

## Posterior inference
burnin = res[[1]]$burn
iter = res[[1]]$iter
ppm <- matrix(0, 100, 100)
current_post <- max(res[[1]]$posteriors)
for(i in 1:4){
  ppm = ppm+res[[i]]$ppm
  if(max(res[[i]]$posteriors) > current_post){
    current_post = max(res[[i]]$posteriors)
    L_map = which(res[[i]]$gamma_map > 0) # Select MAP based on posteriors
  }
}

ppm <- ppm / 4
z_ppm <- minbinder(ppm, method = "comp")$cl

## Get landmark positions
L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)


## Plot deer shape with landmarks
pc <- as.data.frame(dat)
colnames(pc) <- c("x", "y")
landmark_ppm <- pc[L_ppm, ] %>% mutate(method = "BasyesLASA (PPM)")
landmark_map <- pc[L_map, ] %>% mutate(method = "BasyesLASA (MAP)")
landmarks <- landmark_ppm %>%
  rbind(landmark_map)

p_deer <- ggscatter(landmarks, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(
    data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
    size = 0.5, color = "black"
  ) +
  geom_point(data = landmarks, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(data = landmarks, aes(x = x, y = y), colour = "#FC4E07", fill = NA, size = 0.8) +
  facet_wrap(~method) +
  theme(
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_rect(
      colour = "grey",
      fill = "grey"
    ), panel.border = element_rect(fill = NA)
  )
p_deer
ggsave("demo/deer_application.png")
```
![BayesLASA applied to a complex shape (deer) from MPEG-7](demo/deer_application.png)

## Collaborators
* Cong Zhang
* Guanghua Xiao
* Chul Moon
* Min Chen
* Qiwei Li
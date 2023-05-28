##############################################################
#* Simulation study
#* Generate simulated polygonal chains for the manuscript
##############################################################


# generate simulated polygonal chain --------------------------------------

#****************
#* file path
#****************
code_file = "code/landmark_detection"
input <- "data/simulated_data"
output = "data/simulated_data"

#****************
#* load lib and functions
#****************

library(mcclust) 
library(Rcpp)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
registerDoParallel(6)

source(file.path(code_file,"sim_polygon_gaussian.R"))
sourceCpp(file.path(code_file,"MCMC_shape.cpp"))
source(file.path("code/toolbox/functions.R"))

#****************
#*parameter
#****************

Kn = c(4, 5, 6)
sigma = c(0.5, 1, 2)
nn = c(100, 150, 200, 300, 500)

#****************
#*generate polygonal chains
#****************

for(K in Kn){
  for(pn in nn){
    for(equil in c(T, F)){
      
      for(seedd in c(1:25, 27:51)){
        fname = paste0("Normal_", K, "_equil_", equil, "_pn_", pn, "_seed_", seedd)
        if(file.exists(file.path(output, "sim_data", paste0(fname, ".Rdata")))){
          next
        }
        set.seed(2020+seedd)
        edgel = runif(1, 50, 100)
        dat = sim_randon_polygon_generator(K, equil = equil, sigma = sigma, pn = pn,
                                           edgel = edgel, kernel = "normal",
                                           l = 1,  trim = 0.05, rotation = T, seed = 2020+seedd)
        temp = pc_normalizor(dat[,c(1, 2)])
        pc_normalized <- temp$pc;
        center <- temp$center;
        scale <- temp$scale;
        dist <- dist_calculator(pc_normalized, cumsum = TRUE);
        sm001 <- pc_gsmoother(pc_normalized, dist, sigma = 0.01);
        sm001 = pc_normalizor(sm001)$pc
        polyg = list(original_dat = dat, normalized = temp, smooth001 = sm001)
        save(polyg, file = file.path(output, paste0(fname, ".Rdata")))
      }
    }
  }
}



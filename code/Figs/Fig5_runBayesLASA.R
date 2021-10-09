#############
##deer plot in MPGE
##############

#***************
#*file path
#***************
code_file = "code/landmark_detection/"
output = "application"
fig.output = "figs"
data.loc = "data"
fname = "MPEG7closed.mat"

#***************
#*load package
#***************
library(R.matlab)
library(mcclust) 
library(Rcpp)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
registerDoParallel(6)
source(file.path(code_file,"sim_polygon_guassian.R"))
sourceCpp(file.path(code_file,"MCMC_shape.cpp"))
source(file.path("code/toolbox/functions.R"))



f = file.path(data.loc,"demo" ,fname)
mdat <- readMat(f)
mdat = mdat$C.cl

dir.create(file.path(output, "Application_deer"))

#****************************************************************************
#*run 8 MCMC chain, check convergence, select the 4 chain with hightest cor
#*1. run MCMC chain
#****************************************************************************
k = 461
fold = 1000
bs = c(0.01, 0.001, 0.0001, 0.00001)
for(beta_sigma in bs){
  res = foreach(i = 1:8) %dopar%{
    
    pc = cbind(mdat[1,,k], mdat[2,,k])
    temp = pc_normalizor(pc)
    dat = temp$pc;
    dat = dat[-nrow(dat),]
    n = nrow(dat)
    
    est.K = round(n/100)
    
    if(est.K <3){
      est.K = 3
    }
    gamma_i = generate_gamma(n, est.K)
    
    
    MCMC_shape(dat,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
               alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
    
  }
  
  fname = paste0("MPEG7closed_", "deer_", k, "_beta_sigma_",beta_sigma, ".Rds")
  saveRDS(res, file = file.path(output, "Application_deer", fname))
}

#****************************************************************************
#*2. calculate correlation and convergence
#****************************************************************************

ff = dir(file.path(output, "Application_deer"), pattern = "deer_461.*Rds")
library(corrplot)
conv_cor = list()
for(i in 1:length(ff)){
  f = ff[i]
  
  res = readRDS(file.path(output, "Application_deer", f))
  bs = unlist(str_split(f, "_"))[6]
  bs = gsub(".Rds", "", bs)
  used_chain = sel_chain[bs][[1]]
  bs = as.numeric(bs)
  
  ##ppm
  
  PPI = matrix(NA, nrow = 100, ncol = 8)
  
  burnin = res[[1]]$burn
  iter = res[[1]]$iter
  
  for(j in 1:8){
    PPI[,j] = colSums(res[[j]]$PPI)/dim(res[[j]]$PPI)[1]
  }
  
  cor_mat = cor(PPI, method = "pearson")
  colnames(cor_mat) = 1:8
  rownames(cor_mat) = 1:8
  cor_vec = cor_mat[upper.tri(cor_mat)]
  conv_cor[[i]] = cbind.data.frame(bs = bs, cor_vec)
  #pdf(file = file.path(fig.output, paste0("Converge_application_deer.pdf")), width = 4.1, height = 4)
   p_converge = corrplot(cor_mat,method="circle", type = "upper", order = "hclust" ,
                         #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                         addCoef.col = "black", # Add coefficient of correlation,
                         addgrid.col = "grey90",
                         tl.col="black")
  #dev.off()
  
}

conv_cor = do.call("rbind.data.frame", conv_cor)
write.csv(conv_cor, file = file.path(output, "Application_deer", "deer_461_convergence.csv"))

#****************************************************************************
#*3. use select chain to get landmarks
#****************************************************************************

sel_chain = list(c(1, 4, 6, 7),
                 c(1, 5, 6, 8),
                 c(1, 4, 5, 8), 
                 c(2, 4, 7, 8))
names(sel_chain) = c(0.00001, 0.0001, 0.001, 0.01)

for(f in ff){
  res = readRDS(file.path(output, "Application_deer", f))
  bs = unlist(str_split(f, "_"))[6]
  bs = gsub(".Rds", "", bs)
  
  used_chain = sel_chain[bs][[1]]
  bs = as.numeric(bs)
  ##ppm
  ppm = 0
  PPI = 0
  Ks = c()
  burnin = res[[1]]$burn
  iter = res[[1]]$iter
  L_map = which(res[[1]]$gamma_map > 0)
  current_post = max(res[[1]]$posteriors)
  for(i in used_chain){
    ppm = ppm+res[[i]]$ppm
    PPI = PPI+res[[i]]$PPI
    Ks = c(Ks, res[[1]]$Ks[(burnin+1):iter])
    if(max(res[[i]]$posteriors) > current_post){
      current_post = max(res[[i]]$posteriors)
      L_map = which(res[[i]]$gamma_map > 0)
    }
  }
  
  ppm = ppm/4
  z_ppm <- minbinder(ppm, method = "comp")$cl
  L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
  
  len = dim(ppm)[1]
  
  ##ppi
  ppi = PPI/4
  
  data = list("L_ppm" = L_ppm,"L_map" = L_map, "PPI" = ppi, "Ks" = Ks, n = len)
  fname = paste0("MPEG7closed_", "deer_461_beta_sigma_",bs, "select_chain.Rdata")
  save(data, file = file.path(output, "Application_deer", fname))
  
}

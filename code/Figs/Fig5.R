#************************************************
#*Fig.5. application deer
#***********************************************

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(latex2exp)
library(ggpubr)

#***************
#*file path
#***************
code_file = "code/landmark_detection/"
output = "application"
fig.output = "figs"
data.loc = "data"
fname = "MPEG7closed.mat"

source(file.path(code_file,"sim_polygon_guassian.R"))
sourceCpp(file.path(code_file,"MCMC_shape.cpp"))
source(file.path("code/toolbox/functions.R"))


##input data

f = "~/Dropbox/shape_analysis_cong/sandbox/ALDUQ/DetailedCode/MPEG7closed.mat"
mdat <- readMat(f)
mdat = mdat$C.cl
pc = cbind(mdat[1,,k], mdat[2,,k])
#selected deep plot index = 461
k = 461


# BayesLASA ---------------------------------------------------------------


##############################
#### run BayesLASA
##############################


#iteration 100*N
fold = 100
bs = c(0.01, 0.001, 0.0001, 0.00001)
for(beta_sigma in bs){
  res = foreach(i = 1:4) %dopar%{
    
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
  
  fname = paste0("MPEG7closed_", "deer_", k, "_beta_sigma_",beta_sigma, "_MCMC.Rds")
  saveRDS(res, file = file.path(output, "Application_deer", fname))
}

################################################################
### read a data from bayesLASA and get Landmarks
################################################################
ff = dir(file.path(output, "Application_deer"), pattern = "MCMC.Rds")
for(f in ff){
  bs = unlist(str_split(f, "_"))[6]
  bs = as.numeric(gsub(".Rds", "", bs))
  
  #####
  ##read bayes lasa result
  #####
  ##ppm
  ppm = 0
  PPI = 0
  Ks = c()
  burnin = res[[1]]$burn
  iter = res[[1]]$iter
  L_map = which(res[[1]]$gamma_map > 0)
  current_post = max(res[[1]]$posteriors)
  for(i in 1:4){
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
  fname = gsub("MCMC.Rds", "MCMC_L.Rdata", f)
  save(data, file = file.path(output, "Application_deer",  fname))
}


################################################
## plot BayesLASA
################################################

landmarks = list()
Ks = list()
ff = dir(file.path(output, "Application_deer"), pattern = "MCMC_L.Rdata")
for(f in ff){
  load(file.path(output, "Application_deer", f),verbose = T)
  
  beta_sigma = as.numeric( unlist(str_split(f, "_"))[6])
  
  L_ppm = data$L_ppm
  L_ppm = sort(L_ppm)
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_ppm,], cmethod = "BayesLASA (PPM)",
                                                 beta_sigma = beta_sigma, Method = "BayesLASA")))
  ####################
  ###map
  ####################
  L_map = data$L_map
  L_map= sort(L_map)
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_map,], cmethod = "BayesLASA (MAP)",
                                                 beta_sigma = beta_sigma, Method = "BayesLASA")))
  
  ####################
  ###Ks
  ####################
  
  Ks = c(Ks, list(cbind.data.frame(K = data$Ks, beta_sigma = beta_sigma, Method = "BayesLASA")))
}

landmarks = do.call("rbind", landmarks)
Ks = do.call("rbind", Ks)
pc = as.data.frame(pc)
colnames(pc) = c("x", "y")
colnames(landmarks)[c(1, 2)] = c("x", "y")
landmarks$beta_sigma = factor(landmarks$beta_sigma, levels = rev(c(1e-05, 1e-04, 0.001, 0.01)))

landmark_ppm = subset(landmarks, cmethod == "BayesLASA (PPM)")


p_ppm = ggscatter(landmark_ppm, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "black") +
  geom_point(data = landmark_ppm, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = landmark_ppm, aes(x = x, y = y), colour = "#FC4E07",fill = NA, size = 0.8)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 

p_ppm_panel = facet(p_ppm, facet.by = c("beta_sigma"),nrow = 1)

Ks$beta_sigma = factor(Ks$beta_sigma, levels = rev(c(1e-05, 1e-04, 0.001, 0.01)))

p_Ks = gghistogram(Ks, x = "K", y = "..density..", fill = "lightgray", binwidth = 1)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 
p_Ks_panel = facet(p_Ks, facet.by = c("beta_sigma"),nrow = 1, scales = "free_x")

landmark_map = subset(landmarks, cmethod == "BayesLASA (MAP)")

p_map = ggscatter(landmark_map, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "black") +
  geom_point(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07",fill = NA, size = 0.8)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 

p_map_panel = facet(p_map, facet.by = c("beta_sigma"),nrow = 1)

p = ggarrange(p_Ks_panel, p_ppm_panel, p_map_panel,
              nrow = 3)

ggsave(p, file = file.path(fig.output, "Fig5_bayeslasa.pdf"), width = 12, height = 9)

# ALDUQ deer --------------------------------------------------------------


#################################################################
##run aldque use their codes on k = 461, 4 chains, iter = 100N
##alduq point_confidentce interval
##################################################################


#################
ff = dir(file.path(output, "Application_deer"), pattern = "ALDUQ.*mat")
landmarks_alduq = list()
Ks_alduq = list()
pn_t = 101
for(f in ff){
  mdat <- readMat(file.path(output, "Application_deer", f))
  
  k_occur = unlist(lapply(mdat$a, function(x) length(x[[1]])))
  
  k_freq = sort(table(k_occur),decreasing=TRUE)[1]
  k_freq = as.numeric(names(k_freq))
  k_index = which(k_occur == k_freq)
  
  k_post = NULL
  for(i in k_index){
    k_post = cbind(k_post, mdat$a[[i]][[1]])
  }
  Lp = numeric(nrow(k_post))
  L_post = apply(k_post, 2, function(x) close_point(x, pn_t-1))

  for(k in 2:ncol(L_post)){
    dis_1 = sum(round_distance(L_post[,k], L_post[,k-1], pn_t-1))
    L_shift = c(lead(L_post[,k]), L_post[1,k])
    L_shift = L_shift[!is.na(L_shift)]
    dis_2 = sum(round_distance(L_shift, L_post[,k-1], pn_t-1))
    if(dis_1 > dis_2){
      L_post[,k] = L_shift
    }
  }
  
  L_jasa_post = sort(apply(L_post, 1, function(x) {  as.numeric(names(sort(table(x),decreasing=TRUE)[1])) }))
  lambda = gsub(".mat", "",gsub(".*lambda_", "", f))
  landmarks_alduq = c(landmarks_alduq, list(cbind.data.frame(pc[L_jasa_post,],
                                                             Method = "ALDUQ", lambda = lambda, cmethod = "ALDUQ")))
  
  Ks_alduq = c(Ks_alduq , list(data.frame(K = k_occur, lambda = lambda, cmethod = "ALDUQ")))
  
}


landmarks_alduq = do.call("rbind", landmarks_alduq)
Ks_alduq = do.call("rbind",  Ks_alduq)
pc = as.data.frame(pc)
colnames(pc) = c("x", "y")
colnames(landmarks_alduq)[c(1, 2)] = c("x", "y")

Ks_alduq$lambda = factor(Ks_alduq$lambda, levels = c("0.00001", "0.0001", "0.001", "0.01"))

p_Ks = gghistogram(Ks_alduq, x = "K", y = "..density..", fill = "lightgray", binwidth = 1)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 
p_Ks_panel = facet(p_Ks, facet.by = c("lambda"),nrow = 1, scales = "free_x")
landmarks_alduq$lambda = factor(landmarks_alduq$lambda, levels = c("0.00001", "0.0001", "0.001", "0.01"))


p_alduq = ggscatter(landmarks_alduq, x = "x", y = "y", color = "#00AFBB",facet.by = "lambda",nrow = 1, size = 2, shape = 17)+
  geom_polygon(data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "black") +
  geom_point(data = landmarks_alduq, aes(x = x, y = y), colour = "#00AFBB", size = 2, shape = 17)+
  geom_polygon(data = landmarks_alduq, aes(x = x, y = y), colour = "#00AFBB",fill = NA, size = 0.8)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 
p = ggarrange(p_Ks_panel, p_alduq,nrow = 2)

ggsave(p, file = file.path(fig.output,"Fig5_alduq.pdf"), width = 12, height = 6)

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
library(R.matlab)

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

f = file.path(data.loc,"demo" ,fname)
mdat <- readMat(f)
mdat = mdat$C.cl
#selected deep plot index = 461
k = 461
pc = cbind(mdat[1,,k], mdat[2,,k])



# BayesLASA ---------------------------------------------------------------


##############################
#### run BayesLASA
##############################


#iteration 100*N
bs = c(0.01, 0.001, 0.0001, 0.00001)


################################################################
### read a data from bayesLASA and get Landmarks
################################################################
ff = dir(file.path(output,"Application_deer"), pattern = "*select_chain.Rdata")
for(f in ff){
  bs = unlist(stringr::str_split(f, "_"))[6]
  bs = as.numeric(gsub("select", "", bs))
  
}


################################################
## plot BayesLASA
################################################

landmarks = list()
Ks = list()
ff = dir(file.path(output,"Application_deer"), pattern = "*select_chain.Rdata")
for(f in ff){
  load(file.path(output, "Application_deer", f),verbose = T)
  
  beta_sigma = unlist(stringr::str_split(f, "_"))[6]
  beta_sigma = as.numeric(gsub("select", "", beta_sigma))
  
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

ff = dir(file.path(output, "Application_deer"), pattern = "ALDUQ.*csv")
landmarks_alduq = list()
Ks_alduq = list()
for(f in ff){
  df_sum = read.csv(file.path(output, "Application_deer", f))
  lambda = gsub("_sum.csv", "",gsub(".*lambda_", "", f))
  L_jasa_post = sort(close_point(df_sum$Median, 100))
  landmarks_alduq = c(landmarks_alduq,
                      list(cbind.data.frame(pc[L_jasa_post,], Method = "ALDUQ", lambda = lambda, cmethod = "ALDUQ")))
  
  f2 = gsub("_sum.csv", ".mat", f)
  mmdat <- readMat(file.path(output, "Application_deer",  f2))
  k_occur = unlist(lapply(mmdat$a, function(x) length(x[[1]])))
  
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


p_alduq = ggscatter(landmarks_alduq, x = "x", y = "y", color = "#00AFBB",facet.by = "lambda",nrow = 1)+
  geom_polygon(data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "black") +
  geom_point(data = landmarks_alduq, aes(x = x, y = y), colour = "#00AFBB", size = 2)+
  geom_polygon(data = landmarks_alduq, aes(x = x, y = y), colour = "#00AFBB",fill = NA, size = 0.8)+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 


p = ggarrange(p_Ks_panel, p_alduq,nrow = 2)


ggsave(p, file = file.path(fig.output,"Fig5_alduq.pdf"), width = 12, height = 6)

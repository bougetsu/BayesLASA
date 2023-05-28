#************************************************
#*Fig. 3. Simulated study, 
#*example polygonal chain with identified landmarks
#*landmarks were identified by BayesLASA, ALDUQ and convex hull
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
input = "manuscript_reproducibility/data/simulated_data/"
input_raw = "data/simulated_data/"
fig.output = "manuscript_reproducibility/figures_and_tables/"

source("code/toolbox/functions.R")

ffs = c("Normal_4_equil_FALSE_pn_150_seed_30",
        "Normal_5_equil_TRUE_pn_150_seed_16",
        "Normal_6_equil_TRUE_pn_150_seed_16")


# list to save landmark, CI and polygon
landmarks = list()
landmarks_cred_list = list()
pcs = list()
#* load data
for(fname in ffs){
  #get num of landmarkpoints
  s = unlist(strsplit(fname, "_"))
  num = s[2]
  ##load chull result
  load(file.path(input, "Chull", paste0(fname, "_hpts.Rdata")), verbose = T)
  ##load ADLUQ result
  load(file.path(input, "ALDUQ", paste0(fname, "peak_L_ci.rdata")),verbose = T)
  ##load BayesLASA result
  load(file = file.path(input, "BayesLASA", paste0(fname, "_MCMC_L.Rdata")), verbose = T)
  ##load original data
  load(file = file.path(input_raw, paste0(fname, ".Rdata")), verbose = T)
  
  pc = polyg$original_dat[,c(1, 2)]
  pcs = c(pcs, list(cbind.data.frame(pc, K = num, sigma2 = polyg$original_dat[,4])))
  ####################
  ##get true L and gamma
  ####################
  gamma_true = polyg$original_dat[,3]
  gamma_true = gamma_true[-length(gamma_true)]
  n = length(gamma_true)
  L_true = which(gamma_true != 0)
  
  landmarks = c(landmarks,
                list(cbind.data.frame(pc[L_true,], Method = "True landmarks", loc = "Lamdmarks",
                                      cloc = "Lamdmarks",K = num, cmethod = "True landmarks")))
  
  ####################
  ### PPM
  ####################
  L_ppm = data$L_ppm
  L_ppm = sort(L_ppm)
  landmarks = c(landmarks,
                list(cbind.data.frame(pc[L_ppm,],  Method = "BayesLASA (PPM)",  loc = "Lamdmarks",
                                      cloc = "Lamdmarks",K = num, cmethod = "BayesLASA")))
  
  df_cred_map = data$df_cred_ppm
  landmarks_cred = list()
  for(j in 1:nrow(df_cred_map)){
    
    if( (df_cred_map$lwr[j] <= df_cred_map$upr[j]) & ( df_cred_map$upr[j] - df_cred_map$lwr[j] < 90) ){
      landmarks_cred[[j]] = data.frame(x = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],1],
                                       y = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],2],
                                       Method = "BayesLASA (PPM)",  loc = paste0("Credible ",j), cloc = "Credible",K = num, cmethod = "BayesLASA")
    }else{
      start = max(df_cred_map$lwr[j], df_cred_map$upr[j])
      end = min(df_cred_map$lwr[j], df_cred_map$upr[j])
      landmarks_cred[[j]] = data.frame(x = pc[c(start:149, 1:end),1],
                                       y = pc[c(start:149, 1:end),2],
                                       Method = "BayesLASA (PPM)", loc = paste0("Credible ",j),  cloc = "Credible", K = num, cmethod = "BayesLASA")
    }
  }
  
  landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
  
  ####################
  ### MAP
  ####################
  L_ppi = data$df_cred_ppi$L
  L_ppi = sort(L_ppi)
  landmarks = c(landmarks,
                list(cbind.data.frame(pc[L_ppi,],  Method = "BayesLASA (MAP)",  loc = "Lamdmarks",
                                      cloc = "Lamdmarks",K = num, cmethod = "BayesLASA")))
  
  df_cred_map = data$df_cred_ppi
  landmarks_cred = list()
  for(j in 1:nrow(df_cred_map)){
    
    if( (df_cred_map$lwr[j] <= df_cred_map$upr[j]) & ( df_cred_map$upr[j] - df_cred_map$lwr[j] < 90) ){
      landmarks_cred[[j]] = data.frame(x = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],1],
                                       y = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],2],
                                       Method = "BayesLASA (MAP)",  loc = paste0("Credible ",j), cloc = "Credible",K = num, cmethod = "BayesLASA")
    }else{
      start = max(df_cred_map$lwr[j], df_cred_map$upr[j])
      end = min(df_cred_map$lwr[j], df_cred_map$upr[j])
      landmarks_cred[[j]] = data.frame(x = pc[c(start:149, 1:end),1],
                                       y = pc[c(start:149, 1:end),2],
                                       Method = "BayesLASA (MAP)", loc = paste0("Credible ",j),  cloc = "Credible", K = num, cmethod = "BayesLASA")
    }
  }
  
  landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
  
  ####################
  ####ALDUQ
  ####################
  df_L = as.data.frame(df_L)
  L_jasa_post = sort(df_L$medtheta)
  landmarks = c(landmarks,
                list(cbind.data.frame(pc[L_jasa_post,],  Method = "ALDUQ",  loc = "Lamdmarks",
                                      cloc = "Lamdmarks",K = num, cmethod = "ALDUQ")))
  landmarks_cred = list()
  for(j in 1:nrow(df_L)){
    
    if( (df_L$lwr[j] <= df_L$upr[j]) & ( df_L$upr[j] - df_L$lwr[j] < 90) ){
      landmarks_cred[[j]] = data.frame(x = pc[df_L$lwr[j]:df_L$upr[j],1],
                                       y = pc[df_L$lwr[j]:df_L$upr[j],2],
                                       Method = "ALDUQ",  loc = paste0("Credible ",j), cloc = "Credible",K = num, cmethod = "ALDUQ")
    }else{
      start = max(df_L$lwr[j], df_L$upr[j])
      end = min(df_L$lwr[j], df_L$upr[j])
      landmarks_cred[[j]] = data.frame(x = pc[c(start:149, 1:end),1],
                                       y = pc[c(start:149, 1:end),2],
                                       Method = "ALDUQ", loc = paste0("Credible ",j),  cloc = "Credible", K = num, cmethod = "ALDUQ")
    }
  }
  landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
  
  ###############
  ##convex hull
  #############
  
  L_chull = sort(unique(hpts))
  L_chull[L_chull == (n+1)] = 1
  L_chull = sort(L_chull)
  landmarks = c(landmarks,
                list(cbind.data.frame(pc[L_chull,],  Method = "Convex Hull",  loc = "Lamdmarks",
                                      cloc = "Lamdmarks",K = num, cmethod = "Convex Hull")))
}
#* combine result
landmarks = do.call("rbind", landmarks)
landmarks_cred_list = do.call("rbind",landmarks_cred_list)
pcs =  do.call("rbind", pcs)  %>% mutate(K = paste0("K = ", K))
#* save landmark points and intervals
landmarks2 = landmarks %>% rbind(landmarks_cred_list) %>%
  mutate(ccmethod = factor(Method,levels =  c("True landmarks", "BayesLASA (MAP)", "BayesLASA (PPM)", "ALDUQ", "Convex Hull")),
         Method = factor(cmethod, levels = c("True landmarks", "BayesLASA", "ALDUQ", "Convex Hull")),
         K = paste0("K = ", K))

#* plot
p_example = ggscatter(data = subset(landmarks2, cloc == "Lamdmarks"), x = "x", y = "y", color  = "Method", shape = "Method",
                      size = 0.01,
                      facet.by = c("K", "ccmethod"),
                      palette = c("#BB3099", "#FC4E07", "#00AFBB",  "#E7B800"))+
  geom_point(data = as.data.frame(pcs), aes(x = x, y = y), color = "grey70", size = 0.5, alpha = 0.6)+
  geom_polygon(data = as.data.frame(pcs), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "grey70") +
  geom_point(data = subset(landmarks2, cloc == "Lamdmarks"), aes(x = x, y = y, color = Method, shape = Method), size = 2)+
  geom_path(data = subset(landmarks2, cloc == "Credible"), aes(x = x, y = y, group = loc, color = Method),
            size = 2, alpha = 0.3)+
  scale_shape_manual(values = c(9, 15, 17, 19))+
  theme(strip.text=element_text(size=12, colour="grey50"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(size=0.5)) 

ggsave(p_example, file = file.path(fig.output, paste0("figure_3.pdf")), width = 12, height = 9)




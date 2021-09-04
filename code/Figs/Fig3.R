#************************************************
#*Fig. 3. Simulated study, 
#*example polygonal chain with identified landmarks
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
input = "simulation"
fig.output = "figs"

ffs = c("Normal_4_equil_FALSE_pn_150_seed_30_smoothed_001peak_L.rdata",
        "Normal_5_equil_TRUE_pn_150_seed_16_smoothed_001peak_L.rdata",
        "Normal_6_equil_TRUE_pn_150_seed_16_smoothed_001peak_L.rdata")

landmarks = list()
pcs = list()
for(f in ffs){
  
  fname = gsub("_smoothed.*", "", f)
  s = unlist(strsplit(f, "_"))
  dataset = s[1]
  num = s[2]
  ##load chull result
  load(file.path(input, "chull", paste0(fname, "_hpts.Rdata")), verbose = T)
  ##load ADLUQ result
  load(file.path(input, "ALDUQ", f),verbose = T)
  ##load BayesLASA result
  load(file = file.path(input, "MCMC", paste0(fname, "_MCMC_L.Rdata")), verbose = T)
  ##load original data
  load(file = file.path(input, "sim_data", paste0(fname, ".Rdata")), verbose = T)
  
  
  pc = polyg$original_dat[,c(1, 2)]
  pcs = c(pcs, list(cbind.data.frame(pc, K = num, sigma2 = polyg$original_dat[,4])))
  ####################
  ##get true L and gamma
  ####################
  gamma_true = polyg$original_dat[,3]
  gamma_true = gamma_true[-length(gamma_true)]
  n = length(gamma_true)
  L_true = which(gamma_true != 0)
  
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_true,], Method = "True landmarks", K = num, cmethod = "True landmarks")))
  
  L_ppm = data$L_ppm
  
  L_ppm_5 = L_window(L_true, L_ppm, n, window = 5)
  
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_ppm,], Method = "BayesLASA (PPM)", K = num, cmethod = "BayesLASA")))
  
  
  ####################
  ###ppi
  ####################
  L_ppi_ci = sort(data$L_ppi)
  
  ppi_df = NULL
  flag = F
  for( i in 1:length(L_ppi_ci)){
    if(flag == FALSE){
      flag = T
      low = L_ppi_ci[i]
    }
    if(flag == T & (i == length(L_ppi_ci) | L_ppi_ci[i] +1 != L_ppi_ci[min(i+1,length(L_ppi_ci)) ] )){
      flag = F
      high = L_ppi_ci[i]
      ppi_df = rbind(ppi_df, c(low, high))
    }
  }
  
  if(ppi_df[1,1] ==1 & ppi_df[nrow(ppi_df), 2] == n){
    ppi_df[1, 1] = ppi_df[nrow(ppi_df), 1]
    ppi_df = ppi_df[-nrow(ppi_df),]
  }
  
  L_ppi = Ci2Bin(ppi_df, L_true)
  L_ppi_5 = L_window(L_true,  L_ppi, n, window = 5)
  
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_ppi,], Method = "BayesLASA (MAP)", K = num, cmethod = "BayesLASA")))
  ####################
  ########JASA method
  ####################
  
  
  L_jasa_post = sort( L_jasa_post)
  L_jasa_5 = sort(L_window(L_true, L_jasa_post, n, window = 5))
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_jasa_post,], Method = "ALDUQ", K = num, cmethod = "ALDUQ")))
  ###############
  ##convex hull
  #############
  
  L_chull = sort(unique(hpts))
  L_chull[L_chull == (n+1)] = 1
  L_chull = sort(L_chull)
  L_chull_5 = L_window(L_true, L_chull, n, window = 5)
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_chull,], Method = "Convex Hull", K = num, cmethod = "Convex Hull")))
}
landmarks = do.call("rbind", landmarks)

landmarks2 = landmarks %>% mutate(ccmethod = factor(Method,levels =  c("True landmarks", "BayesLASA (MAP)", "BayesLASA (PPM)", "ALDUQ", "Convex Hull")),
                                  Method = factor(cmethod, levels = c("True landmarks", "BayesLASA", "ALDUQ", "Convex Hull")),
                                  K = paste0("K = ", K))
pcs =  do.call("rbind", pcs)

pcs = pcs %>% mutate(K = paste0("K = ", K))
pcs$sigma2 = factor(pcs$sigma2)

p_example = ggscatter(landmarks2, x = "x", y = "y", color  = "Method", shape = "Method", facet.by = c("K", "ccmethod"),
                      palette = c("#BB3099", "#FC4E07", "#00AFBB", "#708238"))+
  geom_point(data = as.data.frame(pcs), aes(x = x, y = y), color = "grey70", size = 0.5)+
  geom_polygon(data = as.data.frame(pcs), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "grey70") +
  geom_point(data = landmarks2, aes(x = x, y = y, color = Method, shape = Method), size = 2)+
  theme(strip.text=element_text(size=12, colour="grey50"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(size=0.5)) 


p_example = ggscatter(data = landmarks2, x = "x", y = "y", color  = "Method", shape = "Method", size = 0.01,facet.by = c("K", "ccmethod"),
                      palette = c("#BB3099", "#FC4E07", "#00AFBB",  "#E7B800"))+
  geom_point(data = as.data.frame(pcs), aes(x = x, y = y), color = "grey70", size = 0.5)+
  geom_polygon(data = as.data.frame(pcs), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "grey70") +
  geom_point(data = landmarks2, aes(x = x, y = y, color = Method, shape = Method), size = 2)+
  scale_shape_manual(values = c(9, 15, 17, 19))+
  theme(strip.text=element_text(size=12, colour="grey50"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(size=0.5)) 

ggsave(p_example, file = file.path(fig.output, paste0("Simulation_example_pn150_selected.pdf")), width = 12, height = 9)




#* plot deer add CI

#read data
library(ggpubr)
library(dplyr)
library(R.matlab)
f = "~/Dropbox/shape_analysis_cong/sandbox/ALDUQ/DetailedCode/MPEG7closed.mat"
mdat <- readMat(f)
mdat = mdat$C.cl
pc = cbind(mdat[1,,461], mdat[2,,461])
landmarks = list()


################
##BayesLASA
################
landmarks = list()
landmarks_cred_list = list()
ff = dir(file.path(output, "Application_deer_new"), pattern = "deer_461")
ff = grep("Rdata", ff, value = T)
for(f in ff){
  load(file.path(output, "Application_deer_new", f),verbose = T)
  
  beta_sigma = as.numeric(gsub("select", "", unlist(str_split(f, "_"))[6]))
  
  L_ppm = data$L_ppm
  L_ppm = sort(L_ppm)
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_ppm,], cmethod = "BayesLASA (PPM)",
                                                 beta_sigma = beta_sigma, Method = "BayesLASA")))
  
  
  df_cred_ppm = get_cred_interval(L_ppm, data$PPI, 100, pp = 0.05)
  landmarks_cred = list()
  for(j in 1:nrow(df_cred_ppm)){
    
    if( (df_cred_ppm$lwr[j] <= df_cred_ppm$upr[j]) & ( df_cred_ppm$upr[j] - df_cred_ppm$lwr[j] < 90) ){
      landmarks_cred[[j]] = data.frame(x = pc[df_cred_ppm$lwr[j]:df_cred_ppm$upr[j],1],
                                       y = pc[df_cred_ppm$lwr[j]:df_cred_ppm$upr[j],2],
                                       Method = "BayesLASA (PPM)",  loc = paste0("Credible ",j), cloc = "Credible",beta_sigma  = beta_sigma , cmethod = "BayesLASA")
    }else{
      start = max(df_cred_ppm$lwr[j], df_cred_ppm$upr[j])
      end = min(df_cred_ppm$lwr[j], df_cred_ppm$upr[j])
      landmarks_cred[[j]] = data.frame(x = pc[c(start:99, 1:end),1],
                                       y = pc[c(start:99, 1:end),2],
                                       Method = "BayesLASA (PPM)", loc = paste0("Credible ",j),  cloc = "Credible", beta_sigma  = beta_sigma,  cmethod = "BayesLASA")
    }
  }
  
  landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
  
  ####################
  ###map
  ####################
  L_map = data$L_map
  L_map= sort(L_map)
  landmarks = c(landmarks, list(cbind.data.frame(pc[L_map,], cmethod = "BayesLASA (MAP)",
                                                 beta_sigma = beta_sigma, Method = "BayesLASA")))
  
  df_cred_map = get_cred_interval(L_map, data$PPI, 100, pp = 0.05)
  landmarks_cred = list()
  for(j in 1:nrow(df_cred_map)){
    
    if( (df_cred_map$lwr[j] <= df_cred_map$upr[j]) & ( df_cred_map$upr[j] - df_cred_map$lwr[j] < 90) ){
      landmarks_cred[[j]] = data.frame(x = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],1],
                                       y = pc[df_cred_map$lwr[j]:df_cred_map$upr[j],2],
                                       Method = "BayesLASA (MAP)",  loc = paste0("Credible ",j), cloc = "Credible",beta_sigma  = beta_sigma , cmethod = "BayesLASA")
    }else{
      start = max(df_cred_map$lwr[j], df_cred_map$upr[j])
      end = min(df_cred_map$lwr[j], df_cred_map$upr[j])
      landmarks_cred[[j]] = data.frame(x = pc[c(start:99, 1:end),1],
                                       y = pc[c(start:99, 1:end),2],
                                       Method = "BayesLASA (MAP)", loc = paste0("Credible ",j),  cloc = "Credible", beta_sigma  = beta_sigma,  cmethod = "BayesLASA")
    }
  }
  
  landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
}

landmarks = do.call("rbind", landmarks)

pc = as.data.frame(pc)
colnames(pc) = c("x", "y")
colnames(landmarks)[c(1, 2)] = c("x", "y")
landmarks$beta_sigma = factor(landmarks$beta_sigma,
                              levels = rev(c(1e-05, 1e-04, 0.001, 0.01)),
                              labels = c(expression(beta[sigma] == "0.01"),
                                         expression(beta[sigma] == "0.001"),
                                         expression(beta[sigma] == "0.0001"),
                                         expression(beta[sigma] == "0.00001")))
landmarks_cred_list = do.call("rbind", landmarks_cred_list)
landmarks_cred_list$beta_sigma = factor(landmarks_cred_list$beta_sigma,
                       levels = rev(c(1e-05, 1e-04, 0.001, 0.01)),
                       labels = c(expression(beta[sigma] == "0.01"),
                                  expression(beta[sigma] == "0.001"),
                                  expression(beta[sigma] == "0.0001"),
                                  expression(beta[sigma] == "0.00001")))


landmark_map = subset(landmarks, cmethod == "BayesLASA (MAP)")
landmark_map = landmark_map %>% mutate(loc = "landmark", cloc = "landmark")
landmarks_cred_list_map = subset(landmarks_cred_list, Method == "BayesLASA (MAP)")
landmarks_cred_list_map = landmarks_cred_list_map %>% rename(cmethod2 = Method, Method = cmethod) %>% rename(cmethod = cmethod2)
p_map = ggscatter(landmark_map, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
               size = 0.5, color = "black") +
  geom_point(data = as.data.frame(pc), aes(x = x, y = y),
             colour = "grey20", size = 0.9)+
  geom_path(data = subset(landmarks_cred_list_map, cloc == "Credible"), 
            aes(x = x, y = y, group = loc, color = Method),
            size = 2, alpha = 0.5)+
  geom_point(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15)+
  geom_polygon(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07",fill = NA, size = 0.8)+
  theme_void()+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(fill=NA)) 
p_map_panel = p_map+ scale_x_continuous(breaks = integer_breaks(n = 4)) +
  facet_wrap( ~ beta_sigma,nrow = 1,
              labeller = labeller(beta_sigma = label_parsed))
p_map_panel


ggsave(p_map_panel, file = file.path("~/Dropbox/shape_analysis_cong/sandbox", "pics", "manuscript",
                           "deer_application_map_ci.pdf"), width = 12, height = 3)

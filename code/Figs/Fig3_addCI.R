library(Rcpp)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
registerDoParallel(6)

ffs = c("Normal_4_equil_FALSE_pn_150_seed_30_smoothed_001peak_L.rdata",
        "Normal_5_equil_TRUE_pn_150_seed_16_smoothed_001peak_L.rdata",
        "Normal_6_equil_TRUE_pn_150_seed_16_smoothed_001peak_L.rdata")



foreach(f = ffs) %dopar%{
  cat(f)
  fname = gsub("_smoothed.*", "", f)
  s = unlist(strsplit(f, "_"))
  dataset = s[1]
  num = s[2]
  load(file = file.path(output, "sim_data", paste0(fname, ".Rdata")), verbose = T)
  dat = polyg$normalized$pc
  dat = dat[-nrow(dat),]
  n = nrow(dat)
  beta_sigma = 1/n
  est.K = round(n/100)
  fold = 100
  if(est.K <3){
    est.K = 3
  }
  gamma_i = generate_gamma(n, est.K)
  a = system.time({r = MCMC_shape(dat,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
                                  alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)})
  ss = unlist(strsplit(f, "_"))
  num = ss[2]
  equil = ss[4]
  pn = ss[6]
  
  ##ppm
  ppms = r$ppm
  z_ppm <- minbinder(r$ppm, method = "comp")$cl
  L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
  ##ppi
  ppi = r$PPI
  ppis = colSums(ppi)/(r$iter-r$burn)
  #L_ppi = which(ppis >= 0.5)
  L_ppi = which(r$gamma_map > 0)
  L_ppi_cor = NULL
  for(i in 1:ncol(ppi)){
    if(i == ncol(ppi)){
      j = 1
    }else{
      j = i+1
    }
    cr = cor.test(ppi[,i], ppi[,j],  alternative = "less")
    if(!is.na(cr$p.value)){
      if(cr$p.value <= 0.05){
        L_ppi_cor = c(L_ppi_cor, i, j)
      }
    }
  }
  L_ppi_cor = unique(c(L_ppi_cor, L_ppi))
  df_cred_ppm = get_cred_interval(L_ppm, ppi, as.numeric(pn)-2, pp = 0.05)
  df_cred_ppi = get_cred_interval(L_ppi, ppi, as.numeric(pn)-2, pp = 0.05)
  data = list("L_ppm" = L_ppm,"L_ppi" = L_ppi_cor, "PPI" = ppis,
              "df_cred_ppm" = df_cred_ppm,
              "df_cred_ppi" = df_cred_ppi)
  save(data, file = file.path(output, "MCMC", paste0(fname,"_MCMC_L.Rdata")))
}
###plot figs and credible interval
landmarks = list()
landmarks_cred_list = list()
pcs = list()
for(f in ffs){
  
  fname = gsub("_smoothed.*", "", f)
  s = unlist(strsplit(f, "_"))
  dataset = s[1]
  num = s[2]
  ##load chull result
  load(file.path(output, "chull", paste0(fname, "_hpts.Rdata")), verbose = T)
  ##load ADLUQ result
  load(file.path(output, "ALDUQ", paste0(fname, "peak_L_ci.rdata")),verbose = T)
  ##load BayesLASA result
  load(file = file.path(output, "MCMC", paste0(fname, "_MCMC_L.Rdata")), verbose = T)
  ##load original data
  load(file = file.path(output, "sim_data", paste0(fname, ".Rdata")), verbose = T)
  
  
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
  ###ppm
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
  ###map
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
  ########JASA method
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
landmarks = do.call("rbind", landmarks)
landmarks_cred_list = do.call("rbind",landmarks_cred_list)
landmarks2 = landmarks %>% rbind(landmarks_cred_list) %>%
  mutate(ccmethod = factor(Method,levels =  c("True landmarks", "BayesLASA (MAP)", "BayesLASA (PPM)", "ALDUQ", "Convex Hull")),
                                  Method = factor(cmethod, levels = c("True landmarks", "BayesLASA", "ALDUQ", "Convex Hull")),
                                  K = paste0("K = ", K))
pcs =  do.call("rbind", pcs)

pcs = pcs %>% mutate(K = paste0("K = ", K))


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



ggsave(p_example, file = file.path("~/Dropbox/shape_analysis_cong/sandbox/pics/manuscript/",
                                   paste0("Simulation_example_pn150_selected_CI.pdf")), width = 12, height = 9)





##run MCMC on map data
library(mcclust) 

library(Rcpp)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doParallel)
library(foreach)
registerDoParallel(2)

code_file = "~/Dropbox/shape_analysis_cong/code/landmark_detection/"
output = "~/Dropbox/shape_analysis_cong/sandbox/output/1225_sim"
#code_file = "/home/cxz163430/Project/shape/code"
#output = "/home/cxz163430/Project/shape/data/1122_sim"
output = "~/Dropbox/shape_analysis_cong/sandbox/output/1225_sim"

#source(file.path(code_file,"sim_polygon_guassian.R"))
sourceCpp(file.path(code_file,"MCMC_shape.cpp"))
source(file.path(code_file,"functions.R"))



##
load("data/states_outline.Rdata", verbose = T)
state_list = names(outline_polygon)
state_pns = unlist(lapply(outline_polygon, nrow))
for(mm in c(27, 42)){
  state_name = state_list[mm]
  original_dat = outline_polygon[[mm]]
  normalized = pc_normalizor(original_dat)
  
  dat = normalized$pc
  dat = dat[-nrow(dat),]
  n = nrow(dat)
  beta_sigma = c(0.01, 0.001, 0.0001, 0.00001)
  names(beta_sigma) = c("0.01","0.001",  "0.0001","0.00001")
  est.K = 4
  fold = max(5*1e4, 100*n)
  
  hpts <- chull(dat)
  hpts = sort(hpts)
  gamma_i = rep(0, n)
  gamma_i[hpts] = 1
  #gamma_i = L2gamma(hpts,n)
  
  foreach(k = 1:4) %dopar% {
    r = MCMC_shape(dat,  iter = fold,  estK = est.K, gamma_i = gamma_i,
                   alpha_sigma = 3, beta_sigma = beta_sigma[k], ppm_store = T)
    
    # ##ppm
    # ppms = r$ppm
    # z_ppm <- minbinder(r$ppm, method = "comp")$cl
    # L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
    # ##ppi
    # ppi = r$PPI
    # ppis = colSums(ppi)/(r$iter-r$burn)
    # #L_ppi = which(ppis >= 0.5)
    # L_ppi = which(r$gamma_map > 0)
    # # L_ppi_cor = NULL
    # # for(i in 1:ncol(ppi)){
    # #   if(i == ncol(ppi)){
    # #     j = 1
    # #   }else{
    # #     j = i+1
    # #   }
    # #   cr = cor.test(ppi[,i], ppi[,j],  alternative = "less")
    # #   if(!is.na(cr$p.value)){
    # #     if(cr$p.value <= 0.05){
    # #       L_ppi_cor = c(L_ppi_cor, i, j)
    # #     }
    # #   }
    # # }
    # # L_ppi_cor = unique(c(L_ppi_cor, L_ppi))
    # df_cred_ppm = get_cred_interval(L_ppm, ppi, n, pp = 0.05)
    # df_cred_ppm = df_cred_ppm %>%
    #   mutate(Method = "PPM", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
    # df_cred_ppi = get_cred_interval(L_ppi, ppi, n, pp = 0.05)
    # df_cred_ppi = df_cred_ppi %>%
    #   mutate(Method = "MAP", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
    # df_cred_ppm_list = c(df_cred_ppm_list,
    #                      list(df_cred_ppm))
    # df_cred_ppi_list = c(df_cred_ppi_list,
    #                      list(df_cred_ppi))
    saveRDS(r, file = file.path("output/map", paste0("state_", mm,"_bs_", names(beta_sigma)[k],"_MCMCres.rds")))
    rm(r)
    }
}


df_cred_ppm_list = list()
df_cred_ppi_list = list()
for(mm in 1:length(state_list)){
  state_name = state_list[mm]
  original_dat = outline_polygon[[mm]]
  normalized = pc_normalizor(original_dat)
  
  dat = normalized$pc
  dat = dat[-nrow(dat),]
  n = nrow(dat)
  beta_sigma = c(0.01, 0.001, 0.0001, 0.00001)
  names(beta_sigma) = c("0.01","0.001",  "0.0001","0.00001")
  
  
  for(k in 1:4){
    r = readRDS(file = file.path("output/map", paste0("state_", mm,"_bs_", names(beta_sigma)[k],"_MCMCres.rds")))
    
    ##ppm
    ppms = r$ppm
    z_ppm <- minbinder(r$ppm, method = "comp")$cl
    L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
    ##ppi
    ppis = colSums(r$PPI)/(r$iter-r$burn)
    #L_ppi = which(ppis >= 0.5)
    L_ppi = which(r$gamma_map > 0)
    # L_ppi_cor = NULL
    # for(i in 1:ncol(ppi)){
    #   if(i == ncol(ppi)){
    #     j = 1
    #   }else{
    #     j = i+1
    #   }
    #   cr = cor.test(ppi[,i], ppi[,j],  alternative = "less")
    #   if(!is.na(cr$p.value)){
    #     if(cr$p.value <= 0.05){
    #       L_ppi_cor = c(L_ppi_cor, i, j)
    #     }
    #   }
    # }
    # L_ppi_cor = unique(c(L_ppi_cor, L_ppi))
    df_cred_ppm = data.frame(L_ppm)
    df_cred_ppm = df_cred_ppm %>%
      mutate(Method = "PPM", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
    df_cred_ppi = data.frame(L_ppi)
    df_cred_ppi = df_cred_ppi %>%
      mutate(Method = "MAP", bs = names(beta_sigma)[k], state = state_name, state_ind = mm)
    df_cred_ppm_list = c(df_cred_ppm_list,
                         list(df_cred_ppm))
    df_cred_ppi_list = c(df_cred_ppi_list,
                         list(df_cred_ppi))
    rm(r)
  }
}

do.call('rbind', df_cred_ppi_list) %>%
  write.csv(file.path("output", "map_df_cred_ppi.csv"))

do.call('rbind', df_cred_ppm_list) %>%
  write.csv(file.path("output", "map_df_cred_ppm.csv"))


# idx = lapply(df_cred_ppi_list, function(x) { (unique(x$state_ind) %in% c(19, 42)) & (unique(x$bs) == "1d10n") })
# df_cred_ppi_list[[166]] = NULL
# df_cred_ppi_list[[189]] = NULL
df_cred_ppi = do.call('rbind', df_cred_ppi_list)
df_cred_ppi = df_cred_ppi %>% rename(L = L_ppi)
df_cred_ppm = do.call('rbind', df_cred_ppm_list)
df_cred_ppm = df_cred_ppm %>% rename(L = L_ppm)
landmarks = rbind(df_cred_ppi, df_cred_ppm)

landmarks_list = list()
landmarks_cred_list = list()
pcs = list()


for(ind in unique(landmarks$state_ind)){
  landmarks2 = landmarks %>% filter(state_ind == ind)
  
  normalized = pc_normalizor(outline_polygon[[ind]])
  dat = normalized$pc
  centroid = st_centroid(st_polygon(list(dat)))
  
  dat_c = dat - rep(centroid, each = nrow(dat))
  
  pc = as.data.frame(dat_c) %>%
    rename(x = long, y = lat) %>%
    mutate(state = unique(landmarks2$state), state_ind = unique(landmarks2$state_ind))
  pns = nrow(outline_polygon[[ind]]) -1
  pcs = c(pcs, list(cbind(pc)))
  landmarks_list = c(landmarks_list,
                     list(cbind.data.frame(pc[landmarks2$L,c(1, 2)],cloc = "Landmarks", loc = "Landmarks",landmarks2)))
  
  
  # landmarks_cred = list()
  # for(j in 1:nrow(landmarks2)){
  #   
  #   if( (landmarks2$lwr[j] <= landmarks2$upr[j]) & ( landmarks2$upr[j] - landmarks2$lwr[j] < 0.4*pns) ){
  #     landmarks_cred[[j]] = data.frame(x = pc[landmarks2$lwr[j]:landmarks2$upr[j],1],
  #                                      y = pc[landmarks2$lwr[j]:landmarks2$upr[j],2],
  #                                      loc = paste0("Credible ",j), cloc = "Credible",landmarks2[j,])
  #   }else{
  #     start = max(landmarks2$lwr[j], landmarks2$upr[j])
  #     end = min(landmarks2$lwr[j], landmarks2$upr[j])
  #     landmarks_cred[[j]] = data.frame(x = pc[c(start:pns, 1:end),1],
  #                                      y = pc[c(start:pns, 1:end),2],
  #                                      loc = paste0("Credible ",j),  cloc = "Credible",landmarks2[j,])
  #   }
  # }
  # landmarks_cred_list = c(landmarks_cred_list, list(do.call("rbind", landmarks_cred)))
}


pcs =  do.call("rbind", pcs)
landmarks_list = do.call("rbind", landmarks_list)
rownames(landmarks_list) = NULL
landmarks_cred_list = do.call("rbind",landmarks_cred_list)
landmarks_plot = rbind(landmarks_list, landmarks_cred_list)

landmarks_plot = landmarks_list
(1/(state_pns-1) < 0.005)
rownames(landmarks_plot) = NULL


stat_ind_list = split(1:48, ceiling(1:48/8))
library(latex2exp)

landmarks_plot = landmarks_plot %>%
  mutate(bs = factor(bs, levels = c("0.01", "0.001", "0.0001", "0.00001"),
                     labels = c(expression(b[sigma] == "0.01"),
                                expression(b[sigma] == "0.001"),
                                expression(b[sigma] == "0.0001"),
                                expression(b[sigma] == "0.00001"))),
         state = factor(state))

# list of state names
landmarks_plot = landmarks_plot %>%
  mutate(state = stringr::str_to_title(state),
                          state_pn = paste0("(n = ",state_pns[state_ind], ")"),
         state_label = paste0(state, "\n", state_pn))%>%
  group_by(state, bs, Method)%>%
  mutate(landmark_n = n(), landmark_n = paste0("K = ", landmark_n)) #number of identified landmarks
pcs = pcs %>%
  mutate(state = stringr::str_to_title(state),
         state_pn = paste0("(n = ",state_pns[state_ind], ")"),
         state_label = paste0(state, "\n", state_pn))


llandmarks_plot %>%
  group_by(state, bs, Method) %>%
  mutate(legend_x = max(x), legend_y = max(y)) %>%
  slice(1)


pdf(file.path("output/", paste0("MAP_rev_foormat.pdf")), height = 14, width = 10)
for(i in 1:length(stat_ind_list)){
  ind = stat_ind_list[[i]]
  landmarks_plot2 = landmarks_plot %>% filter(state_ind %in% ind, Method == "MAP")
  landmark_legend = subset(pcs, state_ind %in% ind) %>%
    ungroup()%>%
    group_by(state) %>%
    mutate(legend_x = max(x), legend_y = max(y)) %>%
    slice(1)
  
  landmark_legend = landmarks_plot2 %>%
    group_by(state, bs) %>%
    slice(1) %>%
    dplyr::select(bs, state, landmark_n) %>%
    left_join(landmark_legend, by = c("state"))
  
  p_map = ggscatter(data = subset(landmarks_plot2, cloc == "Landmarks"),
                    x = "x", y = "y", 
                    group = "state",
                    color = "#FC4E07",
                    size = 1.5, shape = 15)+
    geom_point(data = subset(pcs, state_ind %in% ind),
               aes(x = x, y = y, group = state), color = "grey20",
               size = 0.6)+
    geom_polygon(data = subset(pcs, state_ind %in% ind),
                 aes(x = x, y = y, group = state), fill = NA, linetype = "solid",
                 size = 0.6, color = "grey20") +
    geom_point(data = subset(landmarks_plot2, cloc == "Landmarks"),
               aes(x = x, y = y),
               colour = "#FC4E07", size = 1.5, shape = 15)+
    geom_polygon(data = subset(landmarks_plot2, cloc == "Landmarks"),
                 aes(x = x, y = y),
                 colour = "#FC4E07",fill = NA,
                 size = 0.7, linetype = "longdash", alpha = 0.2)+
    geom_text(data = landmark_legend,
                    aes(x = 0, y = legend_y+0.02, label = landmark_n)) +
    theme_void()+
    theme(strip.text=element_text(size=12, colour = "grey30"),
          strip.background=element_rect(colour="grey", 
                                        fill="grey"),
          strip.text.y = element_text(angle = 90, colour = "grey30"),
          panel.border =element_rect(fill=NA))+
    #theme(aspect.ratio=1)
    coord_quickmap()
  p_maps = p_map + 
    facet_grid(state_label ~ bs, scales = "free",labeller = labeller(bs = label_parsed))
  print(p_maps)
}
dev.off()


texas_pc = pcs %>% filter(state=="Texas") 
texas_pc %>%
  ggplot(aes(x, y)) +geom_polygon()+theme(aspect.ratio=1)

texas_pc %>%
  ggplot(aes(x, y)) +geom_polygon()+coord_quickmap()

for(ind in 1: unique(landmarks_plot$state_ind)){
  landmarks_plot2 = landmarks_plot %>% filter(state_ind %in% ind, Method == "MAP")
  p_example = ggscatter(data = subset(landmarks_plot2, cloc == "Landmarks"), 
                        x = "x", y = "y", color  = "Method", shape = "Method",
                        size = 0.01,
                        facet.by = c("bs"),
                        palette = c("#BB3099", "#FC4E07", "#00AFBB",  "#E7B800"))+
    geom_point(data = subset(pcs, state_ind == ind),
               aes(x = x, y = y), color = "grey70", size = 0.5)+
    geom_polygon(data =  subset(pcs, state_ind == ind),
                 aes(x = x, y = y), fill = NA, linetype = "solid",
                 size = 0.5, color = "grey70") +
    geom_point(data = subset(landmarks_plot2, cloc == "Landmarks"),
               aes(x = x, y = y, color = Method, shape = Method), size = 2)+
    # geom_path(data = subset(landmarks_plot2, cloc == "Credible"), aes(x = x, y = y, group = loc, color = Method),
    #           size = 2, alpha = 0.3)+
    #scale_shape_manual(values = c(9, 15, 17, 19))+
    theme(strip.text=element_text(size=12, colour="grey50"),
          strip.background=element_rect(colour="grey", 
                                        fill="grey"), panel.border =element_rect(size=0.5))+
    ggtitle(paste0(unique(landmarks_plot2$state), " 1/n = ", round(1/(state_pns[ind]-1), 5)))+
    theme(aspect.ratio=1)
  
  p_map = ggscatter(data = subset(landmarks_plot2, cloc == "Landmarks"),
                    x = "x", y = "y", 
                    group = "state",
                    color = "#FC4E07",
                    size = 2, shape = 15,
                    facet.by = c("bs", "state"), scales = "free")+
    geom_point(data = subset(pcs, state_ind %in% ind),
               aes(x = x, y = y, group = state_ind), color = "grey20", size = 0.8)+
    geom_polygon(data = subset(pcs, state_ind %in% ind),
                 aes(x = x, y = y), fill = NA, linetype = "solid",
                 size = 0.8, color = "grey20") +
    geom_point(data = subset(landmarks_plot2, cloc == "Landmarks"),
               aes(x = x, y = y),
               colour = "#FC4E07", size = 2, shape = 15)+
    geom_polygon(data = subset(landmarks_plot2, cloc == "Landmarks"),
                 aes(x = x, y = y),
                 colour = "#FC4E07",fill = NA, size = 0.8, linetype = "dotted")+
    theme_void()+
    theme(strip.text=element_text(size=12, colour="black"),
          strip.background=element_rect(colour="grey", 
                                        fill="grey"), panel.border =element_rect(fill=NA))+
    theme(aspect.ratio=1)
  
  
  
  print(p_example)
}



# Create a first facet variable with examples of math formulas
iris2 <- iris %>%
  mutate(species_math = factor(Species,
                               levels = c("setosa", "versicolor", "virginica"),
                               labels = c("m^2",
                                          expression(beta[sigma] == 0.01),
                                          bquote(pi == .(pi)))))




landmarks_plot = landmarks_plot %>%
  mutate(bs2 = factor(bs, levels = c("0.01", "0.001", "0.0001", "0.00001"),
                     labels = c(expression(beta[sigma] == 0.01),
                                expression(beta[sigma] == "0.001"),
                                expression(beta[sigma] == "0.0001"),
                                expression(beta[sigma] == "0.00001"))))

# Create a second facet variable with mean lengths
# This illustrates how to pass a numeric vector inside a formula
iris_mean <- iris2 %>%
  group_by(Species) %>%
  summarise(across(ends_with("Length"), mean), .groups="drop")

iris2$mean_length <- factor(iris2$Species,
                            levels =  c("setosa", "versicolor", "virginica"),
                            labels = mapply(function(p, s) bquote(bar(p) == .(p) ~ bar(s) ==.(s)),
                                            round(iris_mean$Petal.Length,3), round(iris_mean$Sepal.Length,3)))


iris2 %>%
  ggplot(aes(x = Petal.Length, y = Petal.Width)) +
  geom_point() +
  facet_wrap(species_math ~ mean_length + Species, labeller = labeller(species_math = label_parsed, mean_length = label_parsed))

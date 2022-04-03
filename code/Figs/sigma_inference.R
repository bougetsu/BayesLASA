#* sigma inference
library(ggpubr)
library(tidyr)
library(Rcpp)
source("code/toolbox/functions.R")
sourceCpp("code/landmark_detection/MCMC_shape.cpp")
# simulation example -----------------------------------------------------------------
pic_out = file.path("figs/sigma_infer/")
dir.create(pic_out)
ffs = c("Normal_4_equil_FALSE_pn_150_seed_30",
        "Normal_5_equil_TRUE_pn_150_seed_16",
        "Normal_6_equil_TRUE_pn_150_seed_16")
for(fname in ffs) {
  load(file = file.path("data/manuscript/", paste0("sim_",fname, ".Rdata")))
  L_list = apply(data$ppi[-nrow(data$ppi),], 1, function(x) {which(x != 0)})
  # pc = data$raw_data[,c(1, 2)]
  normalized = pc_normalizor(data$raw_data[,c(1, 2)])
  pc = normalized$pc
  perim = normalized$length
  n = nrow(pc)
  alpha_sigma = 3
  beta_sigma = 1/n
  #split di by segments defined by landmarks
  sigma_true = data$raw_data[,4]
  sigma_mat = lapply(L_list, function(xx){
    P.reduce = PointOnReducedP(xx-1, pc)
    di = get_di(pc, P.reduce, open = T)
    gamma = rep(0, length(di))
    gamma[xx] = 1
    clusters = cumsum(gamma)
    if(clusters[1] == 0){
      clusters[which(clusters == 0)] = clusters[length(clusters)]
    }
    segments = split(di, clusters)
    sigma_c = clusters
    #estimate sigma2
    for(k in sort(unique(clusters))){
      di_k = segments[[k]][-1]
      kn = length(di_k)
      #get sigma from IG dist
      a = alpha_sigma+kn/2
      b = beta_sigma+sum(di_k^2)/2
      sigma2_k = 1/rgamma(1,shape = a,rate = b)
  
      sigma_c[clusters == k] = sqrt(sigma2_k)
      sigma_c[which(clusters == k)[1]] = 0
    }
    return(sigma_c)
  })
  
  sigma_mat = do.call("rbind", sigma_mat)
  sigma_sum = apply(sigma_mat, 2, function(x) return(data.frame(mean = mean(x),
                                                                lwr = quantile(x, 0.025),
                                                                upr = quantile(x, 0.975))))
  sigma_sum = do.call("rbind", sigma_sum)
  sigma_sum =  sigma_sum %>% mutate(index = 1:n, true_val = sigma_true/perim)
  p = ggplot(data =  sigma_sum, aes(x = index, y = mean))+
    geom_line()+
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", 
                            alpha=0.2, 
                            linetype="dashed",
                            color="grey")+
    geom_line(aes(x = index, y = true_val),linetype = "longdash", col = "red")+
    ylab("sigma")+theme_pubr()
  pdf(file= file.path(pic_out, paste0(fname, "_di_dist.pdf")), width = 4, height = 3.5)
  print(p)
  dev.off()
  # di_mat = lapply(L_list, function(xx){
  #   P.reduce = PointOnReducedP(xx-1, pc)
  #   di = get_di(pc, P.reduce, open = T)
  #   return(di)
  #   })
  # di_mat = do.call("cbind", di_mat)
  # sd_td = data.frame(index = 1:nrow(di_mat), di_var = apply(di_mat, 1, var), true_val = sigma2_true)
  # 
  # jpeg(filename = file.path(pic_out, paste0(fname, "di_var.jpg")))
  # p2 = ggplot(data =  sd_td, aes(x = index, y = di_var))+
  #   geom_line()+
  #   geom_line(aes(x = index, y = true_val),linetype = "longdash", col = "red")+
  #   ylab("di_var")+theme_bw()
  # print(p2)
  # dev.off()
  
}

# deer application example -----------------------------------------------------------------
#* to be edited
#* 
#* 
##sigma infer for deer application
ffs = dir("data/manuscript/", pattern = "MPEG7closed_deer")
for(fname in ffs) {
  load(file = file.path("data/manuscript/", fname))
  L_list = apply(data$ppi[-nrow(data$ppi),], 1, function(x) {which(x != 0)})
  pc = data$raw_data[,c(1, 2)]
  n = nrow(pc)
  alpha_sigma = 3
  beta_sigma = 1/(n)
  #split di by segments defined by landmarks
  sigma_mat = lapply(L_list, function(xx){
    P.reduce = PointOnReducedP(xx-1, pc)
    di = get_di(pc, P.reduce, open = T)
    gamma = rep(0, length(di))
    gamma[xx] = 1
    clusters = cumsum(gamma)
    if(clusters[1] == 0){
      clusters[which(clusters == 0)] = clusters[length(clusters)]
    }
    segments = split(di, clusters)
    sigma2_c = clusters
    #estimate sigma2
    for(k in sort(unique(clusters))){
      di_k = segments[[k]][-1]
      kn = length(di_k)
      
      a = alpha_sigma+kn/2
      b = beta_sigma+sum(di_k^2)/2
      sigma2_k = 1/rgamma(1,a,b)
      
      sigma2_c[clusters == k] = sigma2_k
      sigma2_c[which(clusters == k)[1]] = 0
    }
    return(sigma2_c)
  })
  
  sigma_mat = do.call("rbind", sigma_mat)
  sigma_sum = apply(sigma_mat, 2, function(x) return(data.frame(mean = mean(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))))
  sigma_sum = do.call("rbind", sigma_sum)
  sigma_sum =  sigma_sum %>% mutate(index = 1:nrow( sigma_sum))
  p = ggplot(data =  sigma_sum, aes(x = index, y = mean))+
    geom_line()+
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "lightblue", 
                alpha=0.2, 
                linetype="dashed",
                color="grey")+
    ylab("sigma2")+theme_pubr()
  pdf(file= file.path(pic_out, gsub(".Rdata","_di_dist.pdf",fname)), width = 4, height = 3.5)
  print(p)
  dev.off()
}


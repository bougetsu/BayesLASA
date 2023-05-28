#************************************************
#*Fig. S3. Sigma inference of simulated polygon
#* read the simulated polygonal chain and identified landmark poitns
#* infer the sigma
#* plot the estimated sigma interval against true value
#***********************************************
library(dplyr)
library(tidyr)
library(Rcpp)
library(latex2exp)
library(ggpubr)
source("code/toolbox/functions.R")
sourceCpp("code/landmark_detection/MCMC_shape.cpp")

set.seed(9080)
#***************
#*file path
#***************
code_file = "code/landmark_detection/"
input_raw = "data/simulated_data/"
fig.output = "manuscript_reproducibility/figures_and_tables/"

# simulation example -----------------------------------------------------------------

ffs = c("Normal_4_equil_FALSE_pn_150_seed_30",
        "Normal_5_equil_TRUE_pn_150_seed_16",
        "Normal_6_equil_TRUE_pn_150_seed_16")

sigma_list = list()
for(fname in ffs) {
  
  K = gsub("Normal_(\\d)_.*","\\1" , fname)
  
  load(file = file.path(input_raw, paste0("sim_",fname, ".Rdata")))
  # idenfitied landmarks by BayesLASA
  L_list = apply(data$ppi[-nrow(data$ppi),], 1, function(x) {which(x != 0)})
  # normalized polygonal chain
  normalized = pc_normalizor(data$raw_data[,c(1, 2)])
  pc = normalized$pc
  perim = normalized$length
  pc = data$raw_data[,c(1, 2)]
  n = nrow(pc)
  alpha_sigma = 3
  beta_sigma = 1/n
  #true siam
  sigma_true = data$raw_data[,4]
  # infer sigma
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
      b = beta_sigma/(perim^2)+sum(di_k^2)/2

      sigma2_k = 1/rgamma(1,shape = a,rate = b)
  
      sigma_c[clusters == k] = sqrt(sigma2_k)
      sigma_c[which(clusters == k)[1]] = 0
    }
    return(sigma_c)
  })
  
  sigma_mat = do.call("rbind", sigma_mat)
  sigma_sum = apply(sigma_mat, 2,
                    function(x) return(data.frame(mean = mean(x),
                                                  lwr = quantile(x, 0.025),
                                                  upr = quantile(x, 0.975))))
  sigma_sum = do.call("rbind", sigma_sum)
  #* make into data frame
  sigma_sum =  sigma_sum %>%
    mutate(index = 1:n, true_val = sigma_true, K = paste0("K = ", K), perim = perim)
  sigma_list = c(sigma_list, list(sigma_sum))
  
}
sigma_list = do.call("rbind", sigma_list)

#* remove invalid for plot
sigma_plot = sigma_list %>% 
  filter((index < 151) & !(K == "K = 4" & index == 149)) 

#segments for true sigma
seg_sigma = sigma_list %>%
  mutate(pre_sigma = lag(true_val),
         seg_sigma = lead(true_val)) %>%
  filter(true_val == 0) %>%
  mutate(start = index, end = lead(index)) %>%
  filter(start < end) %>%
  mutate(end = ifelse(end >= 149, end-1, end)) %>%
  group_by(K, start) %>%
  mutate(vertical_sigma = max(pre_sigma, seg_sigma, na.rm = T))

p = ggplot(data = sigma_plot, aes(x = index, y = mean))+
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey", 
              alpha=0.2, 
              linetype="solid",
              color="grey20")+
  geom_segment(data = seg_sigma,
               aes(x = start, xend = end, y = seg_sigma, yend = seg_sigma),
               linetype = "solid", col = "red") +
  geom_vline(data = subset(sigma_plot, true_val == 0),
             aes(xintercept = index), linetype = "dashed",color = "red")+
  ylab(TeX("$\\sigma$"))+xlab("Index")+theme_pubr() +facet_wrap(~K)

p = p+scale_x_continuous(name="Index", breaks = c(1, 50, 100, 150)) +
  scale_y_continuous(name=TeX("$\\sigma$"), breaks = seq(0, 5, 0.5))+
  theme(strip.text=element_text(size=12, colour="grey30"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"),
        panel.border =element_rect(size=1, fill=NA)) 
pdf(file= file.path(fig.output, "figure_s3.pdf"), width = 10, height = 3.5)
print(p)
dev.off()

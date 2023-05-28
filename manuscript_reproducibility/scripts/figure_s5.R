#************************************************
#* Fig.5. Compare performance of different method on simulated study 
#* read the raw data and output of BayesLASA, ALDUQ and convex hull on simulated data
#* calculated ARI
#* do violin and boxplot of ARI
#***********************************************

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(latex2exp)
source("code/toolbox/functions.R")

#***************
#*file path
#***************
code_file = "code/landmark_detection/"
input = "manuscript_reproducibility/data/simulated_data/"
input_raw = "data/simulated_data/"
fig.output = "manuscript_reproducibility/figures_and_tables/"
output <- file.path(input, "Summary_statistics")
dir.create(output)


# Summary results ---------------------------------------------------------
#* to get MCC and ARI for each method
#* 
#* NOTE: TO GET THE MCC and ARI may take a while
#* PLEASE BE PATIENT or skip this section to use the ready-to-plot data in the `Summary_statistics` folder
#* 
#* 
# 
# alduq_loc = file.path(input, "ALDUQ")
# ff = dir(alduq_loc, pattern = "peak_L.rdata" )
# for(f in ff){
#   
#   fname = gsub("_smoothed.*", "", f)
#   ##load chull result
#   load(file.path(input, "chull", paste0(fname, "_hpts.Rdata")), verbose = T)
#   ##load ADLUQ result
#   load(file.path(input, "ALDUQ", f),verbose = T)
#   ##load BayesLASA result
#   load(file = file.path(input, "BayesLASA", paste0(fname, "_MCMC_L.Rdata")), verbose = T)
#   ##load original data
#   load(file = file.path(input_raw, paste0(fname, ".Rdata")), verbose = T)
#   
#   res_tb = data.frame(sample = fname, method = c("PPM", "PPM_window5", "PPM_window10", 
#                                                  "PPI", "PPI_window5", "PPI_window10",
#                                                  "ALDUQ",  "ALDUQ_window5", "ALDUQ_window10",
#                                                  "chull", "chull_window5", "chull_window10"),
#                       MCC = numeric(12), ARI = numeric(12))
#   ####################
#   ##get true L and gamma
#   ####################
#   gamma_true = polyg$original_dat[,3]
#   gamma_true = gamma_true[-length(gamma_true)]
#   n = length(gamma_true)
#   L_true = which(gamma_true != 0)
#   
#   
#   #############
#   ###PPM
#   #############
#   L_ppm = data$L_ppm
#   gamma_ppm = rep(0, n)
#   gamma_ppm[L_ppm] = 1
#   MCC_ppm = L2mcc(gamma_true, L_ppm)
#   ARI_ppm = L2adj(gamma_true, L_ppm)
#   
#   res_tb[res_tb$method == "PPM","MCC"] = MCC_ppm
#   res_tb[res_tb$method == "PPM","ARI"] = ARI_ppm
#   
#   
#   L_ppm_5 = L_window(L_true, L_ppm, n, window = 5)
#   gamma_ppm_5 = rep(0, n)
#   gamma_ppm_5[L_ppm_5] = 1
#   MCC_ppm_5 = L2mcc(gamma_true, L_ppm_5)
#   ARI_ppm_5 = L2adj(gamma_true, L_ppm_5)
#   
#   res_tb[res_tb$method == "PPM_window5","MCC"] = MCC_ppm_5
#   res_tb[res_tb$method == "PPM_window5","ARI"] = ARI_ppm_5
#   
#   
#   L_ppm_10 = L_window(L_true, L_ppm, n, window = 10)
#   gamma_ppm_10 = rep(0, n)
#   gamma_ppm_10[L_ppm_10] = 1
#   MCC_ppm_10 = L2mcc(gamma_true, L_ppm_10)
#   ARI_ppm_10 = L2adj(gamma_true, L_ppm_10)
#   
#   res_tb[res_tb$method == "PPM_window10","MCC"] = MCC_ppm_10
#   res_tb[res_tb$method == "PPM_window10","ARI"] = ARI_ppm_10
#   
#   ####################
#   ###ppi
#   ####################
#   L_ppi_ci = sort(data$L_ppi)
#   
#   ppi_df = NULL
#   flag = F
#   for( i in 1:length(L_ppi_ci)){
#     if(flag == FALSE){
#       flag = T
#       low = L_ppi_ci[i]
#     }
#     if(flag == T & (i == length(L_ppi_ci) | L_ppi_ci[i] +1 != L_ppi_ci[min(i+1,length(L_ppi_ci)) ] )){
#       flag = F
#       high = L_ppi_ci[i]
#       ppi_df = rbind(ppi_df, c(low, high))
#     }
#   }
#   
#   if(ppi_df[1,1] ==1 & ppi_df[nrow(ppi_df), 2] == n){
#     ppi_df[1, 1] = ppi_df[nrow(ppi_df), 1]
#     ppi_df = ppi_df[-nrow(ppi_df),]
#   }
#   
#   L_ppi = Ci2Bin(ppi_df, L_true)
#   MCC_ppi = L2mcc(gamma_true, L_ppi)
#   ARI_ppi = L2adj(gamma_true, L_ppi)
#   
#   
#   res_tb[res_tb$method == "PPI","MCC"] = MCC_ppi
#   res_tb[res_tb$method == "PPI","ARI"] = ARI_ppi
#   
#   ######
#   
#   L_ppi_5 = L_window(L_true,  L_ppi, n, window = 5)
#   gamma_ppi_5 = rep(0, n)
#   gamma_ppi_5[L_ppi_5] = 1
#   MCC_ppi_5 = L2mcc(gamma_true, L_ppi_5)
#   ARI_ppi_5 = L2adj(gamma_true, L_ppi_5)
#   
#   res_tb[res_tb$method == "PPI_window5","MCC"] = MCC_ppi_5
#   res_tb[res_tb$method == "PPI_window5","ARI"] = ARI_ppi_5
#   
#   
#   L_ppi_10 = L_window(L_true,  L_ppi, n, window = 10)
#   gamma_ppi_10 = rep(0, n)
#   gamma_ppi_10[L_ppi_10] = 1
#   MCC_ppi_10 = L2mcc(gamma_true, L_ppi_10)
#   ARI_ppi_10 = L2adj(gamma_true, L_ppi_10)
#   
#   res_tb[res_tb$method == "PPI_window10","MCC"] = MCC_ppi_10
#   res_tb[res_tb$method == "PPI_window10","ARI"] = ARI_ppi_10
#   
#   
#   ####################
#   ####ALDUQ
#   ####################
#   
#   
#   L_jasa_post = sort( L_jasa_post)
#   MCC_jasa = L2mcc(gamma_true, L_jasa_post)
#   ARI_jasa = L2adj(gamma_true, L_jasa_post)
#   
#   res_tb[res_tb$method == "ALDUQ","MCC"] = MCC_jasa
#   res_tb[res_tb$method == "ALDUQ","ARI"] = ARI_jasa
#   
#   L_jasa_5 = sort(L_window(L_true, L_jasa_post, n, window = 5))
#   gamma_jasa_5 = rep(0, n)
#   gamma_jasa_5[L_jasa_5] = 1
#   MCC_jasa_5 = L2mcc(gamma_true, L_jasa_5)
#   ARI_jasa_5 = L2adj(gamma_true, L_jasa_5)
#   
#   res_tb[res_tb$method == "ALDUQ_window5","MCC"] = MCC_jasa_5
#   res_tb[res_tb$method == "ALDUQ_window5","ARI"] = ARI_jasa_5
#   
#   
#   L_jasa_10 = L_window(L_true, L_jasa_post, n, window = 10)
#   gamma_jasa_10 = rep(0, n)
#   gamma_jasa_10[L_jasa_10] = 1
#   MCC_jasa_10 = L2mcc(gamma_true, L_jasa_10)
#   ARI_jasa_10 = L2adj(gamma_true, L_jasa_10)
#   
#   res_tb[res_tb$method == "ALDUQ_window10","MCC"] = MCC_jasa_10
#   res_tb[res_tb$method == "ALDUQ_window10","ARI"] = ARI_jasa_10
#   
#   
#   
#   ###############
#   ##convex hull
#   #############
#   
#   L_chull = sort(unique(hpts))
#   L_chull[L_chull == (n+1)] = 1
#   L_chull = sort(L_chull)
#   gamma_chull = rep(0, n)
#   
#   gamma_chull[L_chull] = 1
#   MCC_chull = L2mcc(gamma_true, L_chull)
#   ARI_chull = L2adj(gamma_true, L_chull)
#   
#   res_tb[res_tb$method == "chull","MCC"] = MCC_chull
#   res_tb[res_tb$method == "chull","ARI"] = ARI_chull
#   
#   
#   L_chull_5 = L_window(L_true, L_chull, n, window = 5)
#   gamma_chull_5 = rep(0, n)
#   gamma_chull_5[L_chull_5] = 1
#   MCC_chull_5 = L2mcc(gamma_true, L_chull_5)
#   ARI_chull_5 = L2adj(gamma_true, L_chull_5)
#   
#   res_tb[res_tb$method == "chull_window5","MCC"] = MCC_chull_5
#   res_tb[res_tb$method == "chull_window5","ARI"] = ARI_chull_5
#   
#   
#   L_chull_10 = L_window(L_true, L_chull, n, window = 10)
#   gamma_chull_10 = rep(0, n)
#   gamma_chull_10[L_chull_10] = 1
#   MCC_chull_10 = L2mcc(gamma_true, L_chull_10)
#   ARI_chull_10 = L2adj(gamma_true, L_chull_10)
#   
#   res_tb[res_tb$method == "chull_window10","MCC"] = MCC_chull_10
#   res_tb[res_tb$method == "chull_window10","ARI"] = ARI_chull_10
#   write.csv(res_tb, 
#             file = file.path(output, paste0(fname, "_summary_coef.csv")),
#             quote = F, row.names = F)
# }


# Read and plot data ---------------------------------------------------------
#* Violin compare our, JASA, Convex hull

#*read data
ff = dir(output, pattern = "coef")
df = NULL
for(f in ff){
  s = unlist(strsplit(f, "_"))
  dataset = s[1]
  num = s[2]
  equil = s[4]
  pn = s[6]
  tmp = read.csv(file.path(output, f), stringsAsFactors = F)
  tmp$kernel = dataset
  tmp$num = num
  tmp$equil = equil
  tmp$pn = pn
  if(is.null(df)){
    df = tmp
  }else{
    df = rbind(df, tmp)
  }
}

#* combine data
df_plot = df %>%
  gather(metric, value, 3:4) %>%
  filter(str_detect(method, "window5"))


#* ARI plot
p_ARI = df_plot %>% 
  filter(metric == "ARI") %>%
  mutate(K = paste0("K = ", num), pn = paste0("n = ", pn))%>%
  mutate(cmethod = ifelse(str_detect(method, "chull"), "Convex Hull", 
                          ifelse(str_detect(method, "ALDUQ"), "ALDUQ", 
                                 ifelse(str_detect(method, "PPM"), "BayesLASA (PPM)", "BayesLASA (MAP)"))),
         Method = factor(gsub(" \\(.*\\)", "", cmethod),levels = c("BayesLASA", "ALDUQ", "Convex Hull")),
         cmethod = factor(cmethod, levels = c("BayesLASA (MAP)", "BayesLASA (PPM)", "ALDUQ", "Convex Hull"))) %>%
  ggviolin(x = "cmethod", y = "value", color = "Method",
           facet.by = c("K", "pn"),
           trim = T, add = "boxplot", x.text.angle = 45, 
           palette = c("#FC4E07", "#00AFBB", "#E7B800")) +
  ylab("MCC") +xlab(NULL) +
  rremove("legend")+
  theme(strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="grey", 
                                      fill="grey"), panel.border =element_rect(size=0.5)) 

ggsave(p_ARI, file = file.path(fig.output, paste0("figure_s5.pdf")), width = 12, height = 9)

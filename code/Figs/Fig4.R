#************************************************
#*Fig.4. Simulated study, 
#*violin and boxplot of MCCs and ARIs
#***********************************************

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(latex2exp)


#***************
#*file path
#***************
output = "simulation"
fig.output = "figs"

##Violin compare our, JASA, Convex hull

##read data
ff = dir(file.path(output, "res_sum"), pattern = "coef")
df = NULL
for(f in ff){
  s = unlist(strsplit(f, "_"))
  dataset = s[1]
  num = s[2]
  equil = s[4]
  pn = s[6]
  tmp = read.csv(file.path(output, "res_sum", f), stringsAsFactors = F)
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


############

df_plot = df %>%
  gather(metric, value, 3:4) 
  filter(str_detect(method, "window5"))

p_MCC = df_plot %>% 
  filter(metric == "MCC") %>%
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

ggsave(p_MCC, file = file.path(fig.output, paste0("Simulation_violin_MCC.pdf")), width = 12, height = 9)
ggsave(p_ARI, file = file.path(fig.output, paste0("Simulation_violin_ARI.pdf")), width = 12, height = 9)

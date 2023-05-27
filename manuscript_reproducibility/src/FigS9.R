###################################################################################
## Fig.S9 Comparison of runtime of ALDUQ and BayesLASA on simutalted dataset
#################################################################################
library(dplyr)
library(ggpubr)
#***************
#*file path
#***************
input = "manuscript_reproducibility/data/simulated_data/"
fig.output = "manuscript_reproducibility/figures_and_tables/"

#***************
#* read data
#***************
time_alduq = read.csv(file = file.path(input, "ALDUQ_time_table.csv")) %>%
  dplyr::select(-method) %>%
  mutate(Method = "ALDUQ") 

time_mcmc <- read.csv(file = file.path(input, "MCMC_time_table.csv")) %>%
  mutate(Method = "BayesLASA") %>%
  dplyr::select(-X)

tb_time = time_alduq %>% rbind(time_mcmc) %>%
  mutate(K = factor(pn), Method = factor(Method, levels = c("BayesLASA", "ALDUQ")))
#***************
#* plot
#***************
p_time = ggline(tb_time, x = "K", y = "elapse", color = "Method",
              add =c("mean_se", "boxplot"), error.plot = "linerange",
              ylab = "Runtime (s)", 
              palette = c("#FC4E07", "#00AFBB")) 
ggsave(p_time, file = file.path(fig.output, paste0("Fig_S9.pdf")), width = 10, height = 5)

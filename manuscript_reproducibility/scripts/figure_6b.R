#***********************************************
#* Survival analysis of lung cancer study
#* read the calculated roughness measurement based on BayesLASA output
#* Figure 6b KM plot of high vs low risk of lung cancer patients
#* predicted using model-based roughness measurement 
#* HMM transit prob, Negative to positive
#***********************************************

library(ggplot2)
library(tidyr)
library(dplyr)
library(moments) 
library(survival)
library(survminer)
library(knitr)
library(tab)

#***************
#*file path
#***************
code_file = "code/landmark_detection/"
input = "manuscript_reproducibility/data/real_data_pathology_images/"
rdatloc = file.path(input, "BayesLASA")
fig.output = "manuscript_reproducibility/figures_and_tables"
# Read data ---------------------------------------------------------------

##* load roughness
Roughness <- read.csv(file.path(input, "summary_statistics","Roughness_summary.csv"), row.names = 1)
rough = Roughness %>% dplyr::filter(shape == 1) %>%
  gather(roughness, value, 4:14) %>%
  group_by(sample, BetaSigma, roughness) %>%
  summarise(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), 
            kurtosis = kurtosis(value, na.rm = T), skewness = skewness(value, na.rm = T),
            q95 = quantile(value, 0.95, na.rm = T),q90 = quantile(value, 0.90, na.rm = T),
            q85 = quantile(value, 0.85, na.rm = T), q80 = quantile(value, 0.80, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T), q70 = quantile(value, 0.70, na.rm = T),
            q50 = quantile(value, 0.50, na.rm = T))


##* patient info
##* get sample from clinical data, use sample name to extract K and di
load(file.path(input, "clinical_info.Rdata"))
pat.dat = data %>% dplyr::select(patient_id, slide_id, dead, stage, female, tobacco, survival_time_new) %>%
  dplyr::filter(slide_id %in% rough$sample) %>%
  distinct() %>%
  mutate(slide_id = as.numeric(slide_id))


#* area, parameter and Ks of tumor region
perim = read.csv(file =  file.path(input, "perimeter.csv"))
Ks = read.csv(file.path(input, "LargestK_bs500.csv"), row.names = 1)
areas = read.csv(file = file.path(input, "areas.csv"), row.names = 1)
r_tb = areas %>% group_by(sample) %>%
  filter(area == max(area), sample %in% rough$sample) %>%
  dplyr::select(sample, shape, area) %>%
  left_join(rough, by = "sample") %>%
  left_join(Ks, by = "sample")

ssamples = r_tb %>% filter(sample %in% rough$sample) %>%
  inner_join(pat.dat, by = c("sample" = "slide_id")) %>%
  distinct(patient_id, sample) %>% pull(sample)

tt = Ks %>% left_join(perim, by = c("sample", "shape")) %>%
  filter(sample %in% ssamples) %>%
  mutate(KperLen = K/perim*100) %>%
  gather(Measure, value, c(2, 5, 6))

# Survival analysis with HMM (model_based) -------------------------------------

###############################
##* read the data and reformat
###############################

hmm_loc = file.path(input, "summary_statistics")

Transition <- read.csv(file = file.path(hmm_loc, "HMM_transition.csv"), row.names = 1, stringsAsFactors = F)
Coef <- read.csv(file = file.path(hmm_loc, "HMM_Coef.csv"), row.names = 1, stringsAsFactors = F)
neg_coef = Coef %>% filter(status == "-1") %>%
  dplyr::select(neg_Intercept = Re1..Intercept., neg_sd = Re1.sd, sample, shape, sigma = BetaSigma, segs)
pos_coef = Coef %>% filter(status == "1") %>%
  dplyr::select(pos_Intercept = Re1..Intercept., pos_sd = Re1.sd, sample, shape, sigma= BetaSigma, segs)

Transition2 = Transition %>% 
  dplyr::select(negneg = X.1, neg2pos = X.1to1, pos2neg = X1to.1, pospos = X1, sample, shape, sigma= BetaSigma, model,segs) %>%
  left_join(neg_coef, by = c("sample", "shape", "sigma", "segs"))%>%
  left_join(pos_coef, by = c("sample", "shape", "sigma", "segs"))


Transition2 = Transition2 %>% 
  dplyr::select(negneg, neg2pos, pos2neg, pospos, 
                neg_Intercept, neg_sd, pos_Intercept, pos_sd,
                sample, shape, sigma, model) %>%
  gather(HMM, value, 1:8) %>%
  group_by(sample, sigma, model,HMM) %>%
  summarise(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), 
            kurtosis = kurtosis(value, na.rm = T), skewness = skewness(value, na.rm = T))

Transition_max = areas %>% group_by(sample) %>%
  filter(area == max(area), sample %in% ssamples) %>%
  dplyr::select(sample, shape, area) %>%
  inner_join(Transition2, by = c("sample"))

Transition_tb = Transition_max %>% 
  filter( !(model == "sign" & HMM %in% c("neg_Intercept", "neg_sd", "pos_Intercept", "pos_sd"))) %>%
  left_join(Ks, by = "sample")

# Figure 6b -----------------------------------------------------------------
#######################################################
##* KM plot for LOOV predicted high-low rish using model-based roughness 
##* measurement, HMM transit prob, Negative to postive
#######################################################

HMM_n2p_tb = Transition_tb %>% filter(HMM =="neg2pos", model == "gaussian") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))
##predict LOOV
n = nrow(HMM_n2p_tb)
risk = numeric(n)
for (i in 1:n) {
  ## survival_time_new is defined as the time between biopsy and death or the end of study, which comes first
  ##survival_time_new is defined as the time between the beginning of enroll and death or the end of study, which comes first
  fit <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                 cluster(patient_id)+area+stage+tobacco+female,
               data = HMM_n2p_tb[-i,])
  risk[i] <- predict(fit, HMM_n2p_tb[i, ], type = "risk")
}
HMM_n2p_tb$risk_group = ifelse(risk >= median(risk, na.rm = T), "high", "low")

fit_hmm_risk = survfit(Surv(survival_time_new, dead) ~ risk_group, data = HMM_n2p_tb);
logrank <- survdiff(Surv(survival_time_new, dead) ~ risk_group, data = HMM_n2p_tb);
p_hmm_risk = ggsurvplot(fit_hmm_risk, pval = T, conf.int = T,
                        legend.title = "Predicted risk",
                        legend.labs = c("High", "Low"))
p_hmm_risk
ggsave(file.path(fig.output, "figure_6b.pdf"), width = 4, height = 4)

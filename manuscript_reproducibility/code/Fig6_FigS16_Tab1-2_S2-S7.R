#***********************************************
#* Survival analysis of lung cancer study
#* Fig. 6
#* Table.1, 2 and S2-S7, S16
#***********************************************
#library("depmixS4") #the HMM library we’ll use
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
data.loc =  file.path(input, "processed")
rdatloc = file.path(input, "BayesLASA")

fig.output = "manuscript_reproducibility/figs_tabs"
# Read data ---------------------------------------------------------------

##* load roughness
Roughness <- read.csv(file.path(input, "Roughness_summary.csv"), row.names = 1)
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



# Table 2 -----------------------------------------------------------------
###############################
##* survival analysis with Ra
###############################
##* get data
Ra_tb = r_tb %>% filter(sample %in% ssamples, roughness == "Ra") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))
##* cox model
fit_ra <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
               cluster(patient_id)+area+stage+tobacco+female, data = Ra_tb)
##* print model
fit_ra  %>% tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p"))
# ##*  latex format
# fit_ra  %>% tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p"), latex = T) %>%
#   kable(format = "latex")
fit_ra  %>%
  tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p")) %>%
  write.csv(file.path(fig.output, paste0("Table_1.csv")))

# Fig 6a -----------------------------------------------------------------
################################################
##Fig 6a, KM plot for high-risk group patients
##############################################

##predict LOOV
n = nrow(Ra_tb)
risk = numeric(n)
for (i in 1:n) {
  ## survival_time_new is defined as the time between biopsy and death or the end of study, which comes first
  ##survival_time_new is defined as the time between the beginning of enroll and death or the end of study, which comes first
  fit <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                 cluster(patient_id)+area+stage+tobacco+female,
               data = Ra_tb[-i,])
  risk[i] <- predict(fit, Ra_tb[i, ], type = "risk")
}
Ra_tb$risk_group = ifelse(risk >= median(risk, na.rm = T), "high", "low")

fit_ra_risk = survfit(Surv(survival_time_new, dead) ~ risk_group, data = Ra_tb);
logrank <- survdiff(Surv(survival_time_new, dead) ~ risk_group, data = Ra_tb);
psurve_ra = ggsurvplot(fit_ra_risk, pval = T, conf.int = T,
                    legend.title = "Predicted risk",
                    legend.labs = c("High", "Low"))
psurve_ra
ggsave(file.path(fig.output, "Fig_6a_Ra_predictSurv.pdf"),width = 4, height = 4)

# Table S2-S6 -----------------------------------------------------------------
###############################
##* survival analysis with surface roughness
###############################
measurement <- c("Rq", "Rp", "Rv", "Rx", "RzJIS")
tab_name <- paste0("S", 2:6)
names(tab_name) <- measurement
for(r in measurement){
  cat("\n#############\n", r, "\n#############\n")
  tmp = r_tb %>% filter(sample %in% ssamples,roughness == r) %>%
    inner_join(pat.dat, by = c("sample" = "slide_id")) 
  formula = paste0("Surv(time = survival_time_new, event = dead) ~ mean +sd +kurtosis+skewness+
                   cluster(patient_id)+area+stage+tobacco+female+K")
  fit.coxph1 <- coxph(as.formula(formula),data = tmp)
  fit.coxph1  %>%
    tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p")) %>%
    print()
  fit.coxph1  %>%
    tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p")) %>%
    write.csv(file.path(fig.output, paste0("Table_", tab_name[[r]], ".csv")))
}


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


# Table 2 -----------------------------------------------------------------
#######################################################
##* survival analysis with model-based roughness 
##* measurement, HMM transit prob, Negative to postive
#######################################################

HMM_n2p_tb = Transition_tb %>% filter(HMM =="neg2pos", model == "gaussian") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))

fit_hmm_n2p <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                   cluster(patient_id)+area+stage+tobacco+female, data = HMM_n2p_tb)
fit_hmm_n2p %>%
  tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p")) %>%
  write.csv(file.path(fig.output, paste0("Table_2.csv")))


# Figure 6b -----------------------------------------------------------------
#######################################################
##* KM plot for LOOV predicted high-low rish using model-based roughness 
##* measurement, HMM transit prob, Negative to postive
#######################################################

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
ggsave(file.path(fig.output, "Fig_6a_HMM_predictSurv.pdf"), width = 4, height = 4)

# Table S7 -----------------------------------------------------------------
#######################################################
##* survival analysis with model-based roughness 
##* measurement, HMM transit prob, pos to neg
#######################################################

HMM_p2n_tb = Transition_tb %>% filter(HMM =="pos2neg", model == "gaussian") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))

fit_hmm_p2n <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                   cluster(patient_id)+area+stage+tobacco+female, data = HMM_p2n_tb)
fit_hmm_p2n %>%
  tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p")) %>%
  write.csv(file.path(fig.output, paste0("Table_S7.csv")))


# Fig. S16 -----------------------------------------------------------------
###########################################
##* survival analysis using tumor
##*  boundary roughness as a covariate
###########################################

ff = dir(file.path(input, "Radius_based"), pattern = "TumorBoundaryRoughness_whole_radius.csv")
tb_tumor_rough <- 
  do.call(rbind,
          lapply(file.path(input, "Radius_based", ff),
                 function(x) read.csv(x,  stringsAsFactors = F, row.names = 1))) %>%
  dplyr::filter(shape == 1) 

tbr_tb = areas %>% group_by(sample) %>%
  filter(area == max(area), sample%in% ssamples) %>%
  dplyr::select(sample, shape, area) %>%
  left_join(tb_tumor_rough, by = "sample") 

param = unique(tbr_tb$L)
df_tbr = list()
for(l in param){
  tmp = tbr_tb %>% filter(sample %in% ssamples,L == l) %>%
    inner_join(pat.dat, by = c("sample" = "slide_id")) 
  formula = paste0("Surv(time = survival_time_new, event = dead) ~ tb_roughness+cluster(patient_id)+area+stage+tobacco+female")
  fit.coxph1 <- coxph(as.formula(formula),data = tmp)
  a = summary(fit.coxph1)$coef
  if(ncol(a) >5){
    df_tbr = c(df_tbr, list(cbind.data.frame(l = l, covariate = rownames(a), a)))
  }
}
df_tbr = do.call(rbind,df_tbr)
colnames(df_tbr)[8] = "pval"

p_TBR = df_tbr %>% filter(covariate == "tb_roughness", l <= 200) %>%
  mutate(Logp = log10(pval), L = factor(l, levels = sort(l))) %>%
  ggline(x = "L", y = "pval", xlab = "L", ylab = "p-value", ylim = c(0, 1)) +
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "red") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1))

ggsave(file.path(fig.output, "Fig_S16.pdf"), p_TBR, width = 5, height = 3)




sk = Roughness %>% filter(sample %in% c(sg, sp)) %>%
  dplyr::select(sample, Ra) %>%
  group_by(sample)%>%
  summarise(kurtosis = e1071::kurtosis(Ra, na.rm = T), skewness = e1071::skewness(Ra, na.rm = T))

sk %>% write.csv(file.path(fig.output, "Fig6_kurtosis_skewness.csv"))


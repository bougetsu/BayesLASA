#***********************************************
#*Fig.7, Table 2 & 3 survival analysis of lung cancer study
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
fig.output = "figs"
input = "data/nlst"
raw_dat = file.path(input, "raw")
data.loc =  file.path(input, "processed")
rdatloc = file.path(input, "BayesLASA")


# Read data ---------------------------------------------------------------

##using bs = 500
##select on, to produce table and predict
##add K into model


bs = 500
load(file.path(data.loc, "nlst", "clinical_info.Rdata"))
areas = read.csv(file = file.path(data.loc, "nlst", "areas.csv"), row.names = 1)
perim = read.csv(file =  file.path(data.loc, "nlst", "perimeter.csv"))
Ks = read.csv(file.path(data.loc, "nlst", paste0("LargestK_bs",bs,".csv")), row.names = 1)

ff = dir(file.path(data.loc, "nlst", "rdata_nr"), pattern = "500_Roughness.csv")
##get sample from clinical data, use sample name to extract K and di
ssample = unique(gsub("_bs_.*","", ff))

Roughness <- 
  do.call(rbind,
          lapply(file.path(outloc, "rdata_nr", ff), function(x) read.csv(x, row.names = 1, stringsAsFactors = F)))

pat.dat = data %>% dplyr::select(patient_id, slide_id, dead, stage, female, tobacco, survival_time_new) %>%
  dplyr::filter(slide_id %in% ssample) %>%
  distinct() %>%
  mutate(slide_id = as.numeric(slide_id))
rough = Roughness %>% dplyr::filter(shape == 1) %>%
  gather(roughness, value, 4:14) %>%
  group_by(sample, BetaSigma, roughness) %>%
  summarise(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T), 
            kurtosis = kurtosis(value, na.rm = T), skewness = skewness(value, na.rm = T),
            q95 = quantile(value, 0.95, na.rm = T),q90 = quantile(value, 0.90, na.rm = T),
            q85 = quantile(value, 0.85, na.rm = T), q80 = quantile(value, 0.80, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T), q70 = quantile(value, 0.70, na.rm = T),
            q50 = quantile(value, 0.50, na.rm = T))
r_tb = areas %>% group_by(sample) %>%
  filter(area == max(area), sample%in% ssample) %>% dplyr::select(sample, shape, area) %>%
  left_join(rough, by = "sample") %>%
  left_join(Ks, by = "sample")
r_tb %>% filter(sample %in% ssample) %>%
  inner_join(pat.dat, by = c("sample" = "slide_id")) %>%
  distinct(patient_id, sample) %>%
  ungroup() %>%
  summarise(n_patient  = length(unique(patient_id)), n_slides = length(unique(sample)))

ssamples = r_tb %>% filter(sample %in% ssample) %>%
  inner_join(pat.dat, by = c("sample" = "slide_id")) %>%
  distinct(patient_id, sample) %>% pull(sample)

tt = Ks %>% left_join(perim, by = c("sample", "shape")) %>%
  filter(sample %in% ssamples) %>%
  mutate(KperLen = K/perim*100) %>%
  gather(Measure, value, c(2, 5, 6))


pk = ggplot(tt, aes(x = value, color = Measure, fill = Measure)) +
  geom_histogram(bins = 30, alpha = 0.5) +
  facet_wrap(~Measure, scales = "free_x") +
  theme_bw()

##############
##Fit Sueface roughness Cox
##########

measurement = unique(r_tb$roughness)
df = list()

for(r in measurement){
  tmp = r_tb %>% filter(sample %in% ssample,roughness == r) %>%
    inner_join(pat.dat, by = c("sample" = "slide_id")) 
  formula = paste0("Surv(time = survival_time_new, event = dead) ~ mean +sd +kurtosis+skewness+
                   cluster(patient_id)+area+stage+tobacco+female+K")
  fit.coxph1 <- coxph(as.formula(formula),data = tmp)
  a = summary(fit.coxph1)$coef
  if(ncol(a) >5){
    df = c(df, list(cbind.data.frame(l = l, roughness = r,covariate = rownames(a), a)))
  }
}

df = do.call(rbind,df)
colnames(df)[9] = "pval"


##########################
##Take Ra for example
###########################

Ra_tb = r_tb %>% filter(sample %in% ssample, roughness == "Ra") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))

fit <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
               cluster(patient_id)+area+stage+tobacco+female, data = Ra_tb)


##########################
##Table 2
###########################

##print model
fit %>% tabcoxph(factor.compression = 1, columns = c( "beta","se", "hr","hr.ci", "p"), latex = T) %>%
  kable(format = "latex")



##########################
##Fig 7a
##########################

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

fit2 = survfit(Surv(survival_time_new, dead) ~ risk_group, data = Ra_tb);
logrank <- survdiff(Surv(survival_time_new, dead) ~ risk_group, data = Ra_tb);
psurve = ggsurvplot(fit2, pval = T, conf.int = T,
                    legend.title = "Predicted risk",
                    legend.labs = c("High", "Low"))
ggsave(file.path(fig.output, "Fig7_Ra_predictSurv.pdf"), width = 4, height = 4)


##############
##HMM
##############

csv.loc = file.path(data.loc, "nlst","rdata_nr")
ff = dir(csv.loc, pattern = "*500_HMM_Transition*")
ff_coef = dir(csv.loc, pattern = "500_HMM_Coeff.csv")
##get sample from clinical data, use sample name to extract K and di
ssample = unique(gsub("_bs_.*", "", ff))

Transition <- 
  do.call(rbind,
          lapply(file.path(csv.loc, ff), function(x) read.csv(x, row.names = 1, stringsAsFactors = F)))
Coef <- 
  do.call(rbind,
          lapply(file.path(csv.loc, ff_coef), function(x) read.csv(x, row.names = 1, stringsAsFactors = F)))

neg_coef = Coef %>% filter(status == "-1") %>%
  dplyr::select(neg_Intercept = Re1..Intercept., neg_sd = Re1.sd, sample, shape, sigma = BetaSigma, segs)
pos_coef = Coef %>% filter(status == "1") %>%
  dplyr::select(pos_Intercept = Re1..Intercept., pos_sd = Re1.sd, sample, shape, sigma= BetaSigma, segs)


Transition2 = Transition %>% dplyr::select(negneg = X.1, neg2pos = X.1to1, pos2neg = X1to.1, pospos = X1, sample, shape, sigma= BetaSigma, model,segs) %>%
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
  filter(area == max(area), sample%in% ssample) %>% dplyr::select(sample, shape, area) %>%
  inner_join(Transition2, by = c("sample"))



Transition_tb = Transition_max %>% 
  filter( !(model == "sign" & HMM %in% c("neg_Intercept", "neg_sd", "pos_Intercept", "pos_sd"))) %>%
  left_join(Ks, by = "sample")


########
##Example Negative 2 postive
############

HMM_n2p_tb = Transition_tb %>% filter(sample %in% ssample, HMM =="neg2pos", model == "gaussian") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))

fit_hmm <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                   cluster(patient_id)+area+stage+tobacco+female, data = HMM_n2p_tb)


HMM_p2n_tb = Transition_tb %>% filter(sample %in% ssample, HMM =="pos2neg", model == "gaussian") %>%
  inner_join(pat.dat, by = c("sample" = "slide_id"))

fit_hmm <- coxph(Surv(time = survival_time_new, event = dead) ~mean +sd +kurtosis+skewness+K+
                   cluster(patient_id)+area+stage+tobacco+female, data = HMM_p2n_tb)
fit_hmm


##########################
##Table 3
##########################

##print model
fit_hmm %>% tabcoxph(factor.compression = 1, decimals = 2,
                     columns = c( "beta","se", "hr","hr.ci", "p"), latex = T) %>%
  kable(format = "latex")


##########################
##Fig 7b
##########################

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

fit2 = survfit(Surv(survival_time_new, dead) ~ risk_group, data = HMM_n2p_tb);
logrank <- survdiff(Surv(survival_time_new, dead) ~ risk_group, data = HMM_n2p_tb);
psurve = ggsurvplot(fit2, pval = T, conf.int = T,
                    legend.title = "Predicted risk",
                    legend.labs = c("High", "Low"))
ggsave(file.path(fig.output, "Fig7_RealData_HMM_predictSurv.pdf"), width = 4, height = 4)




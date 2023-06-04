#***********************************************
#* Survival analysis of lung cancer study
#* Fig. S16 cox model using radius-based tumor boundary roughness measurement as predictors
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
#* file path
#***************
code_file <- "code/landmark_detection/"
input <- "manuscript_reproducibility/results/real_data_pathology_images/"
rdatloc <- file.path(input, "BayesLASA")

fig.output <- "manuscript_reproducibility/figures_and_tables"
# Read data ---------------------------------------------------------------

##* patient info
##* get sample from clinical data, use sample name to extract K and di
load(file.path(input, "nlst_clinical_info.Rdata"))
pat.dat <- data %>%
  dplyr::select(patient_id, slide_id, dead, stage, female, tobacco, survival_time) %>%
  distinct() %>%
  mutate(slide_id = as.numeric(slide_id))

#* area, parameter and Ks of tumor region
perim <- read.csv(file = file.path(input, "perimeter.csv"))
Ks <- read.csv(file.path(input, "LargestK_bs500.csv"), row.names = 1)
areas <- read.csv(file = file.path(input, "areas.csv"), row.names = 1)


# Fig. S16 -----------------------------------------------------------------
###########################################
##* survival analysis using tumor
##*  boundary roughness as a covariate
###########################################

ff <- dir(file.path(input, "Radius_based"), pattern = "TumorBoundaryRoughness_whole_radius.csv")
tb_tumor_rough <-
  do.call(
    rbind,
    lapply(
      file.path(input, "Radius_based", ff),
      function(x) read.csv(x, stringsAsFactors = F, row.names = 1)
    )
  ) %>%
  dplyr::filter(shape == 1)

tbr_tb <- areas %>%
  group_by(sample) %>%
  filter(area == max(area)) %>%
  dplyr::select(sample, shape, area) %>%
  left_join(tb_tumor_rough, by = "sample")

param <- unique(tbr_tb$L)
df_tbr <- list()
for (l in param) {
  tmp <- tbr_tb %>%
    filter(L == l) %>%
    inner_join(pat.dat, by = c("sample" = "slide_id"))
  formula <- paste0("Surv(time = survival_time, event = dead) ~ tb_roughness+cluster(patient_id)+area+stage+tobacco+female")
  fit.coxph1 <- coxph(as.formula(formula), data = tmp)
  a <- summary(fit.coxph1)$coef
  if (ncol(a) > 5) {
    df_tbr <- c(df_tbr, list(cbind.data.frame(l = l, covariate = rownames(a), a)))
  }
}
df_tbr <- do.call(rbind, df_tbr)
colnames(df_tbr)[8] <- "pval"

p_TBR <- df_tbr %>%
  filter(covariate == "tb_roughness", l <= 200) %>%
  mutate(Logp = log10(pval), L = factor(l, levels = sort(l))) %>%
  ggline(x = "L", y = "pval", xlab = "L", ylab = "p-value", ylim = c(0, 1)) +
  geom_hline(
    yintercept = 0.05, linetype = "dashed",
    color = "red"
  ) +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1))

ggsave(file.path(fig.output, "figure_s16.pdf"), p_TBR, width = 5, height = 3)

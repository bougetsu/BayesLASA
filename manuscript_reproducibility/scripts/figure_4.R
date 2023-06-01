#************************************************
#* Fig.4. application deer
#* read the raw data and BayesLASA output
#* plot the identified landmarks in deer application case
#* the histooram of landmark numbers K during MCMC
#***********************************************

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(latex2exp)
library(ggpubr)
library(R.matlab)
source("code/toolbox/functions.R")

#***************
#* file path
#***************
code_file <- "code/landmark_detection/"
input <- "manuscript_reproducibility/results/demo_data_mpeg7_deer/"
input_raw <- "data/demo_data_mpeg7_deer/"
fig.output <- "manuscript_reproducibility/figures_and_tables/"

#* load input data of deer shape
fname <- "MPEG7closed.mat"
f <- file.path(input_raw, fname)
mdat <- readMat(f)
mdat <- mdat$C.cl
# selected deep plot index = 461
k <- 461
pc <- cbind(mdat[1, , k], mdat[2, , k])

# BayesLASA ---------------------------------------------------------------
#* BayesLASA has been run with the beta_sigma and other parameters indicated as in the manuscript
#* with iteration 100*N
# bs = c(0.01, 0.001, 0.0001, 0.00001)

################################################
## plot BayesLASA
################################################

landmarks <- list()
Ks <- list()
ff <- dir(file.path(input, "BayesLASA"), pattern = "*select_chain.Rdata")
for (f in ff) {
  load(file.path(input, "BayesLASA", f), verbose = T)

  beta_sigma <- unlist(stringr::str_split(f, "_"))[6]
  beta_sigma <- as.numeric(gsub("select", "", beta_sigma))

  L_ppm <- data$L_ppm
  L_ppm <- sort(L_ppm)
  landmarks <- c(landmarks, list(cbind.data.frame(pc[L_ppm, ],
    cmethod = "BayesLASA (PPM)",
    beta_sigma = beta_sigma, Method = "BayesLASA"
  )))
  ####################
  ### map
  ####################
  L_map <- data$L_map
  L_map <- sort(L_map)
  landmarks <- c(landmarks, list(cbind.data.frame(pc[L_map, ],
    cmethod = "BayesLASA (MAP)",
    beta_sigma = beta_sigma, Method = "BayesLASA"
  )))

  ####################
  ### Ks
  ####################

  Ks <- c(Ks, list(cbind.data.frame(K = data$Ks, beta_sigma = beta_sigma, Method = "BayesLASA")))
}

landmarks <- do.call("rbind", landmarks)
Ks <- do.call("rbind", Ks)
pc <- as.data.frame(pc)
colnames(pc) <- c("x", "y")
colnames(landmarks)[c(1, 2)] <- c("x", "y")
landmarks$beta_sigma <- factor(landmarks$beta_sigma, levels = rev(c(1e-05, 1e-04, 0.001, 0.01)))


#* hist of Ks
Ks$beta_sigma <- factor(Ks$beta_sigma, levels = rev(c(1e-05, 1e-04, 0.001, 0.01)))
p_Ks <- gghistogram(Ks, x = "K", y = "..density..", fill = "lightgray", binwidth = 1) +
  theme(
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_rect(
      colour = "grey",
      fill = "grey"
    ), panel.border = element_rect(fill = NA)
  )
p_Ks_panel <- facet(p_Ks, facet.by = c("beta_sigma"), nrow = 1, scales = "free_x")

#* plot PPM landmarks
landmark_ppm <- subset(landmarks, cmethod == "BayesLASA (PPM)")
p_ppm <- ggscatter(landmark_ppm, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(
    data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
    size = 0.5, color = "black"
  ) +
  geom_point(data = landmark_ppm, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(data = landmark_ppm, aes(x = x, y = y), colour = "#FC4E07", fill = NA, size = 0.8) +
  theme(
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_rect(
      colour = "grey",
      fill = "grey"
    ), panel.border = element_rect(fill = NA)
  )
p_ppm_panel <- facet(p_ppm, facet.by = c("beta_sigma"), nrow = 1)


#* plot MAP landmarks
landmark_map <- subset(landmarks, cmethod == "BayesLASA (MAP)")
p_map <- ggscatter(landmark_map, x = "x", y = "y", color = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(
    data = as.data.frame(pc), aes(x = x, y = y), fill = NA, linetype = "solid",
    size = 0.5, color = "black"
  ) +
  geom_point(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07", size = 2, shape = 15) +
  geom_polygon(data = landmark_map, aes(x = x, y = y), colour = "#FC4E07", fill = NA, size = 0.8) +
  theme(
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_rect(
      colour = "grey",
      fill = "grey"
    ), panel.border = element_rect(fill = NA)
  )
p_map_panel <- facet(p_map, facet.by = c("beta_sigma"), nrow = 1)

#* put together
p <- ggarrange(p_Ks_panel, p_ppm_panel, p_map_panel,
  nrow = 3
)

ggsave(p, file = file.path(fig.output, "figure_4a.pdf"), width = 12, height = 9)

# ALDUQ deer --------------------------------------------------------------

#################################################################
## run ALDUQ use the matlab on k = 461, 4 chains, iter = 100N
##################################################################

#* K histogram
#*
#*
Ks_alduq <- data.frame(
  lambda = c(rep("lambda == 0.00001", 4), rep("lambda == 0.0001", 5), rep("lambda == 0.001", 9), rep("lambda == 0.01", 7)),
  K = c(3:6, 3:7, 13:21, 18:24),
  Probability = c(
    0.75, 0.01, 0.25, 0.01,
    0.005, 0.01, 0.55, 0.425, 0.025,
    0.005, 0.005, 0.04, 0.23, 0.17, 0.48, 0.06, 0.01, 0.005,
    0.17, 0.34, 0.27, 0.16, 0.06, 0.02, 0.005
  )
)

p_Ks_panel <- Ks_alduq %>%
  ggplot(aes(x = K, y = Probability)) +
  geom_bar(stat = "identity", width = 1, fill = "grey90", col = "black") +
  facet_grid(~lambda, scales = "free_x", labeller = label_parsed) +
  theme_bw() +
  labs(x = "") +
  scale_x_continuous(breaks = scales::extended_breaks(n = 7)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15), axis.title = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

p_4b <- p_Ks_panel

ggsave(p_4b, file = file.path(fig.output, "figure_4b.pdf"), width = 12, height = 3)

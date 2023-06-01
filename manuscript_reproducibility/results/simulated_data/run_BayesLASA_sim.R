##############################################################
#* Run BayesLASA on Simulation study
#* 1. read simulated polygonal chains
#* 2. run BayesLASA, convex hull
#* 3. summarize the results from BayesLASA, Convex Hull and ALDUQ
##############################################################


#****************
#* file path
#****************
code_file <- "code/landmark_detection"
input_raw <- "data/simulated_data"
output <- "manuscript_reproducibility/results/simulated_data"

#****************
#* load lib and functions
#****************

library(mcclust)
library(Rcpp)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
registerDoParallel(6)

source(file.path(code_file, "sim_polygon_gaussian.R"))
sourceCpp(file.path(code_file, "MCMC_shape.cpp"))
source(file.path("code/toolbox/functions.R"))

# run BayesLASA -----------------------------------------------------------
ff <- dir(file.path(input_raw), pattern = "Rdata")
tm_df <- NULL

tm_df <- foreach(f = ff) %dopar% {
  cat(f)
  load(file.path(input_raw, f))
  dat <- polyg$normalized$pc
  dat <- dat[-nrow(dat), ]
  n <- nrow(dat)
  beta_sigma <- 1 / n
  est.K <- round(n / 100)
  fold <- 100
  if (est.K < 3) {
    est.K <- 3
  }
  gamma_i <- generate_gamma(n, est.K)
  a <- system.time({
    r <- MCMC_shape(dat,
      iter = n * fold, estK = est.K, gamma_i = gamma_i,
      alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T
    )
  })
  ss <- unlist(strsplit(f, "_"))
  num <- ss[2]
  equil <- ss[4]
  pn <- ss[6]

  ## ppm
  ppms <- r$ppm
  z_ppm <- minbinder(r$ppm, method = "comp")$cl
  L_ppm <- which(diff(c(z_ppm[n], z_ppm)) != 0)
  ## ppi
  ppi <- r$PPI
  ppis <- colSums(ppi) / (r$iter - r$burn)
  # L_ppi = which(ppis >= 0.5)
  L_ppi <- which(r$gamma_map > 0)
  L_ppi_cor <- NULL
  for (i in 1:ncol(ppi)) {
    if (i == ncol(ppi)) {
      j <- 1
    } else {
      j <- i + 1
    }
    cr <- cor.test(ppi[, i], ppi[, j], alternative = "less")
    if (!is.na(cr$p.value)) {
      if (cr$p.value <= 0.05) {
        L_ppi_cor <- c(L_ppi_cor, i, j)
      }
    }
  }
  L_ppi_cor <- unique(c(L_ppi_cor, L_ppi))
  data <- list("L_ppm" = L_ppm, "L_ppi" = L_ppi_cor, "PPI" = ppis)
  fname <- gsub(".Rdata", "_MCMC_L.Rdata", f)
  save(data, file = file.path(output, "BayesLASA", fname))
  #* record time
  data.frame(num, equil, pn, elapse = a[3])
}
tm_df <- do.call("rbind.data.frame", tm_df)
write.csv(tm_df, file = file.path(output, "MCMC_time_table.csv"))


# ALDUQ -------------------------------------------------------------------

#**************************************************
#* ALDUQ was run using matlab codes from the repo https://github.com/jd-strait/ALDUQ of
#* the manuscript "Automatic Detection and Uncertainty Quantification of Landmarks on Elastic Curves"
#*  (https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1527224)
#* and convert into landmark position
#**************************************************

# Convex hull -------------------------------------------------------------

tm_df_curv <- list()
ff <- dir(file.path(input_raw))

for (f in ff) {
  load(file.path(input_raw, f))
  fname <- gsub(".Rdata", "", f)
  pc <- polyg$original_dat[, c(1, 2)]
  s <- unlist(strsplit(fname, "_"))
  num <- s[2]
  equil <- s[4]
  pn <- s[6]
  polyg <- list()

  a <- system.time({
    hpts <- chull(pc)
    hpts <- c(hpts, hpts[1])
  })

  tm_df_curv[[i]] <- data.frame(num, equil, pn = pn, elapse = a[3], method = "Chull")
  save(hpts, file = file.path(output, "Chull", paste0(fname, "_hpts.Rdata")))
  i <- i + 1
}

tm_df_curv_1 <- do.call("rbind", tm_df_curv)
write.csv(tm_df_curv_1, file = file.path(output, "Chull_time_table.csv"))

##############################################################
#* Run BayesLASA on pathology images data
#* 1. read processed pathology image data
#* 2. run BayesLASA **TAKE LONG TIME, SUGGEST TO RUN ON A SERVER WITH PARALLEL
#* 3. summarize the results from BayesLASA
#* 4. calculate distance-based and model-based roughness measurement
##############################################################


#****************
#* file path
#****************
code_file <- "code/landmark_detection"
input_raw <- "data/real_data_pathology_images/"
input <- "manuscript_reproducibility/results/real_data_pathology_images"
output <- "manuscript_reproducibility/results/real_data_pathology_images/BayesLASA/"

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
##************************************
##*
##* TAKE LONG TIME
##* SUGGEST TO RUN ON A SERVER WITH PARALLEL
##*
##************************************
alpha_sigma <- 3
beta_sigma <- 500
est.K <- 10
updateG.p <- c(0.8, 0.1, 0.1)
fold <- 100

ff <- dir(input_raw, "outline.Rdata")
foreach(f = ff) %dopar% {
  load(file.path(input_raw, f))
  dat <- outline_polygon_nr[[1]]
  dat <- dat[-nrow(dat), ]
  n <- nrow(dat)
  res <- foreach(j = 1:4) %dopar% {
    gamma_i <- generate_gamma(n, 10)
    MCMC_shape(dat,
      iter = n * fold, estK = est.K, gamma_i = gamma_i, updateGp = updateG.p,
      alpha_sigma = alpha_sigma, beta_sigma = beta_sigma
    )
  }

  fname <- paste0(gsub("_outline.Rdata", "", f), "Sigma_bs500_normal_4chain_largest.rdata")

  saveRDS(res, file = file.path(output, fname))
}



# calculate roughness -----------------------------------------------------

#############
## perimeter
#############

ff <- dir(file.path(input_raw), pattern = "_outline.Rdata")
p <- list()
j <- 1
for (f in ff) {
  load(file.path(input_raw, f))
  sample <- gsub("_outline.Rdata", "", f)
  for (k in 1:length(outline_polygon)) {
    peri <- perimeter(outline_polygon[[k]])
    p[[j]] <- data.frame(sample = sample, shape = k, perim = peri)
    j <- j + 1
  }
}
do.call("rbind", p) %>%
  write.csv(file = file.path(input, "perimeter.csv"), quote = F, row.names = F)
############
## records K
#############

bs <- 500
ff <- dir(output, pattern = "Sigma_bs500_normal_4chain_largest.rdata")
Ks <- list()
j <- 1
for (f in ff) {
  sample <- gsub("Sigma_bs.*", "", f)
  load(file.path(outloc, f), verbose = T)
  Ks[[j]] <- data.frame(sample = sample, K = length(L_ppm[[1]]), shape = 1, BetaSigma = bs)
  j <- j + 1
}
do.call("rbind", Ks) %>%
  write.csv(file = file.path(input, paste0("LargestK_bs", bs, ".csv")))


# calculate roughness and hmm ---------------------------------------------
SurfaceRoughness <- function(y) {
  Ra <- mean(abs(y))
  Rq <- sqrt(mean(y^2))
  Rv <- abs(min(y))
  Rp <- max(y)
  Rx <- Rp + Rv
  Rsk <- mean((y^3) / Rq^3)
  Rku <- mean((y^4) / Rq^4)
  RzJIS <- mean(sort(y, decreasing = T)[1:5] - sort(y)[1:5])
  Rtm <- get_Rtm(y)
  diff_abs_mean <- mean(abs(diff(y)))
  diff_sd <- sd(diff(y))
  di_mean <- mean(di)
  # y = c(y, 0)
  # n = length(dis)
  # Ra = auc(dis, y, absolutearea= T)/dis[n]
  # y2 = ifelse(y >=0, y^2, -y^2)
  # Rq = sqrt(auc(p.dist, y2, absolutearea= T)/dis[n])
  return(data.frame(di_mean, diff_abs_mean, diff_sd, Ra, Rq, Rv, Rp, Rx, Rsk, Rku, RzJIS))
}
library("depmixS4") # the HMM library we’ll use
library(moments)
pc.id.conversion <- function(pc1, pc2) {
  idx.conversion <- numeric(nrow(pc1))
  i <- 1
  j <- 1
  while (i <= nrow(pc1)) {
    if (all(pc1[i, ] == pc2[j, ])) {
      idx.conversion[i] <- j
      i <- i + 1
      j <- j + 1
      if (j > nrow(pc2)) {
        j <- j - nrow(pc2)
      }
    } else {
      j <- j + 1
      if (j > nrow(pc2)) {
        j <- j - nrow(pc2)
      }
    }
  }
  return(idx.conversion)
}

bs <- 500
Roughness_sum <- list()
ff <- dir(output, pattern = "Sigma_bs500_normal_4chain_largest.rdata")
for (f in ff) {
  sample <- gsub("Sigma_bs.*", "", f)
  load(file.path(output, f), verbose = T)
  ##### sample points along polygon?
  load(file.path(input_raw, paste0(sample, "_outline.Rdata")))

  L_nr <- list()
  di <- list()
  roughness <- list()
  transition_mat <- list()
  coeffs <- list()
  for (k in 1:length(L_ppm)) {
    pc1 <- outline_polygon[[k]]
    pc2 <- outline_polygon_nr[[k]]
    L <- L_ppm[[k]]

    idx.conversion <- pc.id.conversion(pc1, pc2)
    L_nr[[k]] <- idx.conversion[L]

    P.reduce <- PointOnReducedP(L_nr[[k]] - 1, pc2)

    di[[k]] <- get_di(pc2, P.reduce)
    yi <- di[[k]] - mean(di[[k]])
    n <- nrow(pc2) - 1
    gamma <- L2gamma(L_nr[[k]] - 1, n)
    segs <- gamma2seg(gamma) + 1
    for (ss in 1:nrow(segs)) {
      if (segs[ss, 1] < segs[ss, 2]) {
        ddi <- di[[k]][segs[ss, 1]:(segs[ss, 2] - 1)]
        # p.project = P.reduce[c(segs[ss,1]:segs[ss,2]),]
        # p.dist = dist_calculator(p.project, T)
      } else {
        ddi <- di[[k]][c(segs[ss, 1]:n, 1:(segs[ss, 2] - 1))]
      }

      if (length(roughness) < k) {
        roughness[[k]] <- data.frame(sample = sample, shape = k, BetaSigma = bs, SurfaceRoughness(ddi), di_mean = mean(ddi))
      } else {
        roughness[[k]] <- rbind(
          roughness[[k]],
          data.frame(sample = sample, shape = k, BetaSigma = bs, SurfaceRoughness(ddi), di_mean = mean(ddi))
        )
      }

      data_hmm <- data.frame(obs = ddi, status = ifelse(ddi >= 0, "1", "-1"))
      set.seed(1)
      ##### fit hmm model with guassian dist family
      mod1 <- depmix(response = obs ~ 1, data = data_hmm, nstates = 2, ntimes = nrow(data_hmm)) # use gaussian() for normally distributed data

      t <- try(
        {
          fit.mod1 <- fit(mod1)
        },
        silent = T
      )
      if (!inherits(t, "try-error")) {
        # cat("a")
        coeff <- summary(fit.mod1)
        ### set S1 < 0
        transit_gaussian <- c(fit.mod1@transition[[1]]@parameters$coefficients, fit.mod1@transition[[2]]@parameters$coefficients)
        if (coeff[1] > coeff[2]) {
          coeff <- coeff[2:1, ]
          transit_gaussian <- transit_gaussian[4:1]
        }

        ## fit hmm model with just status sequence
        mcFit <- markovchainFit(data_hmm$status)
        if (length(mcFit$estimate@transitionMatrix) > 1) {
          transit_sign <- as.vector(c(mcFit$estimate@transitionMatrix[1, ], mcFit$estimate@transitionMatrix[2, ]))
          if (mcFit$estimate@states[1] == "1") {
            transit_sign <- transit_sign[4:1]
          }
        } else {
          transit_sign <- c(1, 0, 0, 1)
        }

        transition <- rbind.data.frame(transit_gaussian, transit_sign)
        colnames(transition) <- c("-1", "-1to1", "1to-1", "1")

        if (length(transition_mat) < k) {
          transition_mat[[k]] <- transition %>% mutate(sample = sample, shape = k, BetaSigma = bs, model = c("sign", "gaussian"), segs = ss)
          coeffs[[k]] <- as.data.frame(coeff) %>% mutate(sample = sample, shape = k, BetaSigma = bs, status = c("-1", "1"), segs = ss)
        } else {
          transition_mat[[k]] <- rbind(
            transition_mat[[k]],
            transition %>% mutate(sample = sample, shape = k, BetaSigma = bs, model = c("sign", "gaussian"), segs = ss)
          )
          coeffs[[k]] <- as.data.frame(coeff) %>%
            mutate(sample = sample, shape = k, BetaSigma = bs, status = c("-1", "1"), segs = ss) %>%
            rbind(coeffs[[k]])
        }
      }
    }
  }
  Roughness_sum <- c(Roughness_sum, list(Roughness))
  Transition_mat <- do.call("rbind", transition_mat)
  Coeff <- do.call("rbind", coeffs)

  write.csv(Transition_mat, file = file.path(input, "summary_statistics", paste0(sample, "_bs_", bs, "_HMM_Transition_mat.csv")))
  write.csv(Coeff, file = file.path(input, "summary_statistics", paste0(sample, "_bs_", bs, "_HMM_Coeff.csv")))
}
Roughness_sum <- do.call("rbind", Roughness_sum)
write.csv(Roughness, file = file.path(input, "summary_statistics", "Roughness_summary.csv"))

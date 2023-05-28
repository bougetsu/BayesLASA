###################################################################################
##* Fig.S6
##* read the simulated data
##* run BayesLASA
##* plot the estimated # of landmkars K_hat vs different
##*  BayesLASA prior of beta_sigma or ALDUQ prior lambda
#################################################################################

##load pkgs

library(doParallel)
library(foreach)
library(R.matlab)
library(ggplot2)
library(latex2exp)

registerDoParallel(6)
source("code/toolbox/functions.R")

#***************
#*file path
#***************
input <- "manuscript_reproducibility/data/simulated_data/"
input_raw <- "data/simulated_data/"
fig.output = "manuscript_reproducibility/figures_and_tables/"


# Run BayesLASA on simulated data -----------------------------------------
####################################################################################
#* The MCMC step will take a while to finish
#* Ready-to-use results are already saved in "manuscript_reproducibility/data/simulated_data/"
#* Uncomment the below code section to re-run the MCMC if needed
####################################################################################
library(mcclust) 
library(Rcpp)
sourceCpp("code/landmark_detection/MCMC_shape.cpp")
# #### Reading simulated data
# # K = 4
# load(file.path(input_raw, "Normal_4_equil_FALSE_pn_100_seed_1.Rdata"))
# sim1 <- as.data.frame(polyg$original_dat)
# 
# # K = 5
# load(file.path(input_raw, "Normal_5_equil_TRUE_pn_100_seed_1.Rdata"))
# sim2 <- as.data.frame(polyg$original_dat)
# 
# # K = 6
# load(file.path(input_raw, "Normal_6_equil_TRUE_pn_100_seed_1.Rdata"))
# sim3 <- as.data.frame(polyg$original_dat)
# 
# # K = 4, 5, 6 (keep changing simulated data for different K)
# # K = 4
# ##pre-process data, scale
# temp1 = pc_normalizor(as.matrix(sim1[,1:2]))
# dat1 = temp1$pc;
# dat1 = dat1[-nrow(dat1),]
# n = nrow(dat1)
# 
# ###set hyper parameter and algorithm setting
# fold = 100
# est.K = round(n/100)
# if(est.K <3){est.K = 3}
# 
# beta_sigma_vec = c(0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
# K_mat <- NULL
# 
# ## 4 MCMC chains
# for(beta_sigma in beta_sigma_vec){
#   #browser()
#   res = foreach(i = 1:4) %dopar%{
#     #generate gamma_0
#     gamma_i = generate_gamma(n, est.K)
#     #run MCMC
#     MCMC_shape(dat1,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
#                alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
#   }
#   
#   ##post infer
#   
#   ##ppm
#   burnin = res[[1]]$burn
#   iter = res[[1]]$iter
#   L_map = which(res[[1]]$gamma_map > 0)
#   current_post = max(res[[1]]$posteriors)
#   
#   ppm = PPI = 0
#   for(i in 1:4){
#     ppm = ppm+res[[i]]$ppm
#     PPI = PPI+res[[i]]$PPI
#     if(max(res[[i]]$posteriors, na.rm = T) > current_post){
#       current_post = max(res[[i]]$posteriors)
#       L_map = which(res[[i]]$gamma_map > 0)
#     }
#   }
#   
#   ppm = ppm/4
#   z_ppm <- minbinder(ppm, method = "comp", max.k = n)$cl
#   L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
#   
#   K_mat <- rbind(K_mat, cbind(length(L_map), length(L_ppm)))
# }
# 
# K_mat.dfr.6 <- data.frame(beta_sigma = beta_sigma_vec, K_map = K_mat[,1], K_ppm = K_mat[,2])
# 
# 
# # K = 5
# ##pre-process data, scale
# temp2 = pc_normalizor(as.matrix(sim2[,1:2]))
# dat2 = temp2$pc;
# dat2 = dat2[-nrow(dat2),]
# n = nrow(dat2)
# 
# ###set hyper parameter and algorithm setting
# fold = 100
# est.K = round(n/100)
# if(est.K <3){est.K = 3}
# 
# beta_sigma_vec = c(0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
# K_mat <- NULL
# 
# ## 4 MCMC chains
# for(beta_sigma in beta_sigma_vec){
#   #browser()
#   res = foreach(i = 1:4) %dopar%{
#     #generate gamma_0
#     gamma_i = generate_gamma(n, est.K)
#     #run MCMC
#     MCMC_shape(dat2,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
#                alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
#   }
#   
#   ##post infer
#   
#   ##ppm
#   burnin = res[[1]]$burn
#   iter = res[[1]]$iter
#   L_map = which(res[[1]]$gamma_map > 0)
#   current_post = max(res[[1]]$posteriors)
#   
#   ppm = PPI = 0
#   for(i in 1:4){
#     ppm = ppm+res[[i]]$ppm
#     PPI = PPI+res[[i]]$PPI
#     if(max(res[[i]]$posteriors, na.rm = T) > current_post){
#       current_post = max(res[[i]]$posteriors)
#       L_map = which(res[[i]]$gamma_map > 0)
#     }
#   }
#   
#   ppm = ppm/4
#   z_ppm <- minbinder(ppm, method = "comp", max.k = n)$cl
#   L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
#   
#   K_mat <- rbind(K_mat, cbind(length(L_map), length(L_ppm)))
# }
# 
# K_mat.dfr.7 <- data.frame(beta_sigma = beta_sigma_vec, K_map = K_mat[,1], K_ppm = K_mat[,2])
# 
# 
# # K = 6
# ##pre-process data, scale
# temp3 = pc_normalizor(as.matrix(sim3[,1:2]))
# dat3 = temp3$pc;
# dat3 = dat3[-nrow(dat3),]
# n = nrow(dat3)
# 
# ###set hyper parameter and algorithm setting
# fold = 100
# est.K = round(n/100)
# if(est.K <3){est.K = 3}
# 
# beta_sigma_vec = c(0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
# K_mat <- NULL
# 
# ## 4 MCMC chains
# for(beta_sigma in beta_sigma_vec){
#   res = foreach(i = 1:4) %dopar%{
#     #generate gamma_0
#     gamma_i = generate_gamma(n, est.K)
#     #run MCMC
#     MCMC_shape(dat3,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
#                alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
#   }
#   
#   ##post infer
#   
#   ##ppm
#   burnin = res[[1]]$burn
#   iter = res[[1]]$iter
#   L_map = which(res[[1]]$gamma_map > 0)
#   current_post = max(res[[1]]$posteriors)
#   
#   ppm = PPI = 0
#   for(i in 1:4){
#     ppm = ppm+res[[i]]$ppm
#     PPI = PPI+res[[i]]$PPI
#     if(max(res[[i]]$posteriors, na.rm = T) > current_post){
#       current_post = max(res[[i]]$posteriors)
#       L_map = which(res[[i]]$gamma_map > 0)
#     }
#   }
#   
#   ppm = ppm/4
#   z_ppm <- minbinder(ppm, method = "comp", max.k = n)$cl
#   L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
#   
#   K_mat <- rbind(K_mat, cbind(length(L_map), length(L_ppm)))
# }
# 
# K_mat.dfr.8 <- data.frame(beta_sigma = beta_sigma_vec, K_map = K_mat[,1], K_ppm = K_mat[,2])
# 
# save(file = file.path(input, "K_prior_convergence.RData"), K_mat.dfr.6, K_mat.dfr.7, K_mat.dfr.8)


# Plot Fig S6 -------------------------------------------------------------
#* This is to show the estimated # of landmkars K_hat vs different
#*  BayesLASA prior of beta_sigma or ALDUQ prior lambda
## load("K_mat.RData") to skip the above codes

load(file.path(input, "K_prior_convergence.RData"))
# Simulation K = 4 Convergence plot

pdf(file.path(fig.output, "figure_s6.pdf"), width = 11, height = 5)

par(mfrow = c(1, 3))
par(mar = c(5, 5, 5, 3))
# Simulation K = 4 Convergence plot

lambda.vec <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 25, 50, 75, 90, 99)
K.vec <- c(rep(4, 3), 11, 12, 14, 15, 42, 55, 63, 67, 67)
ald.dfr <- cbind.data.frame(lambda = lambda.vec, K = K.vec)


matplot(log(K_mat.dfr.6$beta_sigma), log2(K_mat.dfr.6[,2:3]), col = c("orange", "#FC4E07"),
        type = "l", lty = 2, xlim = rev(range(log(K_mat.dfr.6$beta_sigma))), xlab = TeX(r'($\log(b_\sigma)$)'), ylab = "K",
        ylim = c(min(log2(ald.dfr$K)), max(log2(ald.dfr$K))),
        yaxt = "n", cex.main=1.5, cex.lab=1.5, cex.axis=1.5
        #ylim = c(min(ald.dfr$K), max(ald.dfr$K))
)
matpoints(log(K_mat.dfr.6$beta_sigma),log2(K_mat.dfr.6[,2:3]), col = c("orange", "#FC4E07"), pch = 16) # Create first plot
abline(h = log2(4), lty = 2, col = "grey")
axis(2, at = seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5), cex.axis =1.5, labels = round(2^seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5)))

par(new = TRUE)                             # Add new plot
plot(log(ald.dfr$lambda), log2(ald.dfr$K), type = "l", col = "#00AFBB", lty = 2,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "", ylim = c(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)))
     , cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(log(ald.dfr$lambda), log2(ald.dfr$K), pch = 16, col = "#00AFBB")
axis(side = 3, at = pretty(range(log(ald.dfr$lambda))), cex.axis =1.5)      # Add second axis
mtext(TeX(r'($\log(\lambda)$)'), side = 3, line = 3, cex = 1.25)

legend(x = "topleft",          # Position
       legend = c("BayesLASA (MAP)", "BayesLASA (PPM)", "ALDUQ", "True K"),  # Legend texts
       lty = 2,           # Line types
       col = c("orange", "#FC4E07", "#00AFBB", "grey"),           # Line colors
       lwd = 2, title = "Methods", cex = 1, bty = "n")



# Simulation K = 5 Convergence plot

lambda.vec <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 25, 50, 75, 90, 99)
K.vec <- c(5, 5, 6, 7, 8, 18, 20, 43, 56, 61, 65, 69)
ald.dfr <- cbind.data.frame(lambda = lambda.vec, K = K.vec)


matplot(log(K_mat.dfr.7$beta_sigma), log2(K_mat.dfr.7[,2:3]), col = c("orange", "#FC4E07"),
        type = "l", lty = 2, xlim = rev(range(log(K_mat.dfr.7$beta_sigma))), xlab = TeX(r'($\log(\b_\sigma)$)'), ylab = "K",
        yaxt = "n", cex.main=1.5, cex.lab=1.5, cex.axis=1.5
)
matpoints(log(K_mat.dfr.7$beta_sigma),log2(K_mat.dfr.7[,2:3]), col = c("orange", "#FC4E07"), pch = 16) # Create first plot
abline(h = log2(5), lty = 2, col = "grey")
axis(2, at = seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5), cex.axis =1.5, labels = round(2^seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5)))

par(new = TRUE)                             # Add new plot
plot(log(ald.dfr$lambda), log2(ald.dfr$K), type = "l", col = "#00AFBB", lty = 2,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",
     ylim = c(min(log2(K_mat.dfr.7[,2:3])), max(log2(K_mat.dfr.7[,2:3]))),
     cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(log(ald.dfr$lambda), log2(ald.dfr$K), pch = 16, col = "#00AFBB")
axis(side = 3, at = pretty(range(log(ald.dfr$lambda))), cex.axis =1.5)      # Add second axis
mtext(TeX(r'($\log(\lambda)$)'), side = 3, line = 3, cex = 1.25)



# Simulation K = 6 Convergence plot

K.vec <- c(rep(6, 4), 7, 9, 11, 44, 54, 61, 65, 68)
ald.dfr <- cbind.data.frame(lambda = lambda.vec, K = K.vec)

matplot(log(K_mat.dfr.8$beta_sigma), log2(K_mat.dfr.8[,2:3]), col = c("orange", "#FC4E07"),
        type = "l", lty = 2, xlim = rev(range(log(K_mat.dfr.8$beta_sigma))), xlab = TeX(r'($\log(\b_\sigma)$)'), ylab = "K",
        yaxt = "n", cex.main=1.5, cex.lab=1.5, cex.axis=1.5
)
matpoints(log(K_mat.dfr.8$beta_sigma),log2(K_mat.dfr.8[,2:3]), col = c("orange", "#FC4E07"), pch = 16) # Create first plot
abline(h = log2(6), lty = 2, col = "grey")
axis(2, at = seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5), cex.axis =1.5, labels = round(2^seq(min(log2(ald.dfr$K)), max(log2(ald.dfr$K)), length.out = 5)))

par(new = TRUE)                             # Add new plot
plot(log(ald.dfr$lambda), log2(ald.dfr$K), type = "l", col = "#00AFBB", lty = 2,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",
     ylim = c(min(log2(K_mat.dfr.8[,2:3])), max(log2(K_mat.dfr.8[,2:3]))),
     cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(log(ald.dfr$lambda), log2(ald.dfr$K), pch = 16, col = "#00AFBB")
axis(side = 3, at = pretty(range(log(ald.dfr$lambda))), cex.axis =1.5)      # Add second axis
mtext(TeX(r'($\log(\lambda)$)'), side = 3, line = 3, cex = 1.25)
dev.off()

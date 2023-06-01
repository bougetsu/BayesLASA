#######################################################################
##* Fig.S8 sensitivity analysis of simulated polygons
##* 1. read the simulated data
##* 2. perform BayesLASA
##* 3. calculate MCC and ARI
##* 4. plot MCC, ARI, TPR and FPR
##* The final FPR plot in the manuscript was overlaid using Adobe illustrator for better visualization
#######################################################################
library(RColorBrewer)
library(corrplot)

#***************
#* file path
#***************
input <- "manuscript_reproducibility/results/simulated_data/"
input_raw <- "data/simulated_data/"
code_path <- "code/"
fig.output <- "manuscript_reproducibility/figures_and_tables/"

source(file.path(code_path, "toolbox/functions.R"))

# Identify landmark points using different beta_sigma and initial K -------
####################################################################################
##** NOTE: It takes a while to run the MCMC
##** Results are already saved in the 'manuscript_reproducibility/data/simulated_data//sensitivity'
##** uncomment the below code section to re-run in needed.
####################################################################################
# library(mcclust)
# library(Rcpp)
# library(doParallel)
# library(foreach)
# registerDoParallel(4)
# source(file.path(code_path,"landmark_detection/sim_polygon_gaussian.R"))
# sourceCpp(file.path(code_path,"landmark_detection/MCMC_shape.cpp"))

# tm_df = NULL
# ff = dir(file.path(input_raw), pattern = "Normal_5.*pn_150.*Rdata")
# tm_df = foreach(f = ff) %dopar%{
#   cat(f)
#   load(file.path(output,"sim_data", f))
#   dat = polyg$normalized$pc
#   dat = dat[-nrow(dat),]
#   n = nrow(dat)
#   fold = 100
#
#   bss = c(1/(4*n), 1/(4*n), 1/n, 2/n, 4/n)
#   names(bss) = c("1/4n", "1/2n", "1/n", "2/n", "4/n")
#   Ks = c(3, 5, 7, 9, 11)
#
#   ss = unlist(strsplit(f, "_"))
#   num = ss[2]
#   equil = ss[4]
#   pn = ss[6]
#   sensitivity_res = list()
#   for(beta_sigma in bss){
#     for(est.K in Ks){
#       gamma_i = generate_gamma(n, est.K)
#       r = MCMC_shape(dat,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
#                                       alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)
#
#       ##ppm
#       ppms = r$ppm
#       z_ppm <- minbinder(r$ppm, method = "comp")$cl
#       L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
#       ##ppi
#       ppi = r$PPI
#       ppis = colSums(ppi)/(r$iter-r$burn)
#       L_ppi = which(r$gamma_map >= 0.5)
#       L_ppi_cor = NULL
#       for(i in 1:ncol(ppi)){
#         if(i == ncol(ppi)){
#           j = 1
#         }else{
#           j = i+1
#         }
#         cr = cor.test(ppi[,i], ppi[,j],  alternative = "less")
#         if(!is.na(cr$p.value)){
#           if(cr$p.value <= 0.05){
#             L_ppi_cor = c(L_ppi_cor, i, j)
#           }
#         }
#       }
#       L_ppi_cor = unique(c(L_ppi_cor, L_ppi))
#       data = list("beta_sigma" = names(beta_sigma), "estK" = est.K,"L_ppm" = L_ppm,"L_ppi" = L_ppi_cor, "PPI" = ppis)
#       sensitivity_res = c(sensitivity_res , list(data))
#     }
#   }
#   fname = gsub(".Rdata", "_sensitivity_L.Rdata", f)
#   save(sensitivity_res, file = file.path(input, "sensitivity", fname))
#   cat("finished ", fname, "\n")
# }



# Calculate MCC, ARI, TPR and FPR  -------
####################################################################################
##** NOTE: It takes a while to calculated the MCC, ARI, TPR and FPR
##** Results are already saved in the 'manuscript_reproducibility/data/simulated_data//sensitivity'
####################################################################################

ff <- dir(file.path(input, "sensitivity"), pattern = "sensitivity_L.Rdata")
for (f in ff) {
  fname <- gsub("_sensitivity.*", "", f)
  ## load BayesLASA result
  load(file = file.path(input, "sensitivity", f), verbose = T)
  ## load original data
  load(file = file.path(input_raw, paste0(fname, ".Rdata")), verbose = T)

  ####################
  ## get true L and gamma
  ####################
  gamma_true <- polyg$original_dat[, 3]
  gamma_true <- gamma_true[-length(gamma_true)]
  n <- length(gamma_true)
  L_true <- which(gamma_true != 0)

  #* used beta_sigma
  bss <- c("1/4n", "1/2n", "1/n", "2/n", "4/n")

  sens_tb <- list()
  for (i in 1:length(sensitivity_res)) {
    data <- sensitivity_res[[i]]
    bsn <- ceiling(i / 5)
    beta_sigma <- bss[bsn]
    estK <- data$estK
    res_tb <- data.frame(
      sample = fname, beta_sigma = beta_sigma, estK = estK,
      method = c(
        "PPM", "PPM_window5",
        "PPI", "PPI_window5"
      ),
      MCC = numeric(4), ARI = numeric(4),
      TPR = numeric(4), FPR = numeric(4)
    )

    #############
    ### PPM
    #############
    L_ppm <- data$L_ppm
    gamma_ppm <- rep(0, n)
    gamma_ppm[L_ppm] <- 1
    MCC_ppm <- L2mcc(gamma_true, L_ppm)
    ARI_ppm <- L2adj(gamma_true, L_ppm)


    res_tb[res_tb$method == "PPM", "MCC"] <- MCC_ppm
    res_tb[res_tb$method == "PPM", "ARI"] <- ARI_ppm


    TPR_ppm <- L2TPR(gamma_true, L_ppm)
    FPR_ppm <- L2FPR(gamma_true, L_ppm)
    res_tb[res_tb$method == "PPM", "TPR"] <- TPR_ppm
    res_tb[res_tb$method == "PPM", "FPR"] <- FPR_ppm


    L_ppm_5 <- L_window(L_true, L_ppm, n, window = 5)
    gamma_ppm_5 <- rep(0, n)
    gamma_ppm_5[L_ppm_5] <- 1
    MCC_ppm_5 <- L2mcc(gamma_true, L_ppm_5)
    ARI_ppm_5 <- L2adj(gamma_true, L_ppm_5)

    res_tb[res_tb$method == "PPM_window5", "MCC"] <- MCC_ppm_5
    res_tb[res_tb$method == "PPM_window5", "ARI"] <- ARI_ppm_5

    TPR_ppm_5 <- L2TPR(gamma_true, L_ppm_5)
    FPR_ppm_5 <- L2FPR(gamma_true, L_ppm_5)
    res_tb[res_tb$method == "PPM_window5", "TPR"] <- TPR_ppm_5
    res_tb[res_tb$method == "PPM_window5", "FPR"] <- FPR_ppm_5


    ####################
    ### ppi
    ####################
    L_ppi_ci <- sort(data$L_ppi)

    ppi_df <- NULL
    flag <- F
    for (i in 1:length(L_ppi_ci)) {
      if (flag == FALSE) {
        flag <- T
        low <- L_ppi_ci[i]
      }
      if (flag == T & (i == length(L_ppi_ci) | L_ppi_ci[i] + 1 != L_ppi_ci[min(i + 1, length(L_ppi_ci))])) {
        flag <- F
        high <- L_ppi_ci[i]
        ppi_df <- rbind(ppi_df, c(low, high))
      }
    }

    if (ppi_df[1, 1] == 1 & ppi_df[nrow(ppi_df), 2] == n) {
      ppi_df[1, 1] <- ppi_df[nrow(ppi_df), 1]
      ppi_df <- ppi_df[-nrow(ppi_df), ]
    }

    L_ppi <- Ci2Bin(ppi_df, L_true)
    MCC_ppi <- L2mcc(gamma_true, L_ppi)
    ARI_ppi <- L2adj(gamma_true, L_ppi)


    res_tb[res_tb$method == "PPI", "MCC"] <- MCC_ppi
    res_tb[res_tb$method == "PPI", "ARI"] <- ARI_ppi

    TPR_ppi <- L2TPR(gamma_true, L_ppi)
    FPR_ppi <- L2FPR(gamma_true, L_ppi)
    res_tb[res_tb$method == "PPI", "TPR"] <- TPR_ppi
    res_tb[res_tb$method == "PPI", "FPR"] <- FPR_ppi


    ######

    L_ppi_5 <- L_window(L_true, L_ppi, n, window = 5)
    gamma_ppi_5 <- rep(0, n)
    gamma_ppi_5[L_ppi_5] <- 1
    MCC_ppi_5 <- L2mcc(gamma_true, L_ppi_5)
    ARI_ppi_5 <- L2adj(gamma_true, L_ppi_5)

    res_tb[res_tb$method == "PPI_window5", "MCC"] <- MCC_ppi_5
    res_tb[res_tb$method == "PPI_window5", "ARI"] <- ARI_ppi_5

    TPR_ppi_5 <- L2TPR(gamma_true, L_ppi_5)
    FPR_ppi_5 <- L2FPR(gamma_true, L_ppi_5)
    res_tb[res_tb$method == "PPI_window5", "TPR"] <- TPR_ppi_5
    res_tb[res_tb$method == "PPI_window5", "FPR"] <- FPR_ppi_5

    sens_tb[[i]] <- res_tb
  }

  sens_tb <- do.call("rbind.data.frame", sens_tb)

  write.csv(sens_tb,
    file = file.path(input, "sensitivity", paste0(fname, "_summary_coef.csv")),
    quote = F, row.names = F
  )
}

# Load TPR FPR and plot -------

##* load data
ff <- dir(file.path(input, "sensitivity"), pattern = "_summary_coef.csv")

sensitivity_tb <- do.call(
  rbind,
  lapply(file.path(input, "sensitivity", ff), function(x) read.csv(x, stringsAsFactors = F))
)

sensitivity_tb2 <- sensitivity_tb %>%
  group_by(beta_sigma, estK, method) %>%
  summarise(MCC = mean(MCC), ARI = mean(ARI), TPR = mean(TPR), FPR = mean(FPR))

mat <- sensitivity_tb2 %>% filter(method == "PPM_window5")

mat_tpr <- xtabs(TPR ~ beta_sigma + estK, mat)
mat_tpr <- mat_tpr[c(2, 1, 3, 4, 5), ]
mat_fpr <- xtabs(FPR ~ beta_sigma + estK, mat)
mat_fpr <- mat_fpr[c(2, 1, 3, 4, 5), ]


##* make plot

pdf(file = file.path(fig.output, paste0("figure_s8a.pdf")), width = 4.1, height = 4)
corrplot(mat_tpr,
  method = "circle", col = c(
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ),
  is.corr = FALSE, cl.lim = c(0, 1),
  addCoef.col = "black", # Add coefficient of correlation,
  addgrid.col = "grey90",
  tl.col = "black"
)
dev.off()


#* the two FPR plots were overlayed to get the final one
pdf(file = file.path(fig.output, paste0("figure_s8b.pdf")), width = 4.1, height = 4)
corrplot(1 - mat_fpr,
  method = "circle", col = rev(c(
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  )),
  is.corr = FALSE, cl.lim = c(0, 1),
  addgrid.col = "grey90"
)
corrplot(mat_fpr,
  method = "number",
  is.corr = FALSE, cl.lim = c(0, 1),
  col = "black",
  addgrid.col = "grey90",
  number.cex = 0.66,
  number.digits = 3,
  number.font = 1,
  tl.col = "black",
  cl.pos = "n"
)
dev.off()

#****************
#*Simulation study
#*1. generate silulated polygonal chains
#*2. run BayesLASA, convex hull
#*3. summarize the results from BayesLASA, Convex Hull and ALDUQ
#***************



# generate simulated polygonal chain --------------------------------------

###################
##used functions
##################
cal_mcc = function(TP, TN, FP, FN){
  return((TP*TN -FP*FN)/sqrt(  (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ))
}

L2mcc = function(gamma_true, L_pred){
  n = length(gamma_true)
  gamma_pred = rep(0, n)
  gamma_pred[L_pred] = 1
  TP = sum(gamma_pred == 1 & gamma_true == 1)
  TN = sum(gamma_pred == 0 & gamma_true == 0)
  FP = sum(gamma_pred == 1 & gamma_true == 0)
  FN = sum(gamma_pred == 0 & gamma_true == 1)
  return((TP*TN -FP*FN)/sqrt(  (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ))
}

L2adj = function(gamma_true, L_pred){
  n = length(gamma_true)
  gamma_pred = rep(0, n)
  gamma_pred[L_pred] = 1
  tab <- table(gamma_true, gamma_pred)
  f2 <- function(n) n*(n-1)/2
  sum.f2 <- function(v) sum(f2(v))
  marg.1 <- apply(tab,1,sum)
  marg.2 <- apply(tab,2,sum)  
  n	<- sum(tab)
  prod <- sum.f2(marg.1)*sum.f2(marg.2)/f2(n)
  num <- (sum.f2(as.vector(tab))- prod)
  den <- 0.5*(sum.f2(marg.1)+sum.f2(marg.2))-prod
  return(num/den)
}

rround = Vectorize(function(x){
  dig = x %% 1
  if(dig < 0.5){
    return(floor(x))
  }else{
    return(ceiling(x))
  }
})

Ci2Bin = function(df, L_true){
  L_ci = NULL
  for(i in 1:nrow(df)){
    interv1 = c( max(df[i,]):n, 1:min(df[i,]))
    interv2 =min(df[i,]):max(df[i,])
    if(df[i,1] > df[i,2] ){
      interv = interv1
    }else{
      interv = interv2
    }
    if(length( intersect(L_true, interv)) > 0){
      L_ci = c(L_ci, intersect(L_true, interv))
    }else{
      L_ci = c(L_ci, interv[rround(length(interv)/2)] )
    }
    
  }
  return(sort(unique(L_ci)))
}

L_window = function(L_true, L_pred, n, window = 5){
  L_pred = sort(L_pred)
  L_window = NULL
  for(l in L_pred){
    ls = l+-window:window
    ls[ls <1] = ls[ls <1]+n
    ls[ls >n] = ls[ls >n]-n
    
    if(length( intersect(L_true, ls)) == 1){
      L_window = c(L_window, intersect(L_true, ls))
      L_true = setdiff(L_true, intersect(L_true, ls))
    }else if(length( intersect(L_true, ls)) == 0){
      L_window = c(L_window, l )
    }else{
      ll = intersect(L_true, ls)
      ll = ll[which.min(abs(ll -l))]
      L_window = c(L_window, ll)
      L_true = setdiff(L_true, ll)
    }
    
  }
  
  return(L_window)
}


#########################
##prepare simulation data set
##########################



#****************
#* file path
#****************
code_file = "code/landmark_detection/"
output = "simulation"
dir.create(output)
dir.create(file.path(output,"csv"), recursive = T)
dir.create(file.path(output,"sim_data"), recursive = T)
dir.create(file.path(output,"smoothed"), recursive = T)

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

source(file.path(code_file,"sim_polygon_guassian.R"))
sourceCpp(file.path(code_file,"MCMC_shape.cpp"))
source(file.path("code/toolbox/functions.R"))

#****************
#*parameter
#****************

Kn = c(4, 5, 6)
sigma2 = c(0.5, 1, 2)
nn = c(100, 150, 200, 300, 500)

#****************
#*generate polygonal chains
#****************

for(K in Kn){
  for(pn in nn){
    for(equil in c(T, F)){
      
      for(seedd in c(1:25, 27:51)){
        fname = paste0("Normal_", K, "_equil_", equil, "_pn_", pn, "_seed_", seedd)
        if(file.exists(file.path(output, "sim_data", paste0(fname, ".Rdata")))){
          next
        }
        set.seed(2020+seedd)
        edgel = runif(1, 50, 100)
        dat = sim_randon_polygon_generator(K, equil = equil, sigma = sigma2, pn = pn,
                                           edgel = edgel, kernel = "normal",
                                           l = 1,  trim = 0.05, rotation = T, seed = 2020+seedd)
        temp = pc_normalizor(dat[,c(1, 2)])
        pc_normalized <- temp$pc;
        center <- temp$center;
        scale <- temp$scale;
        dist <- dist_calculator(pc_normalized, cumsum = TRUE);
        sm001 <- pc_gsmoother(pc_normalized, dist, sigma = 0.01);
        sm001 = pc_normalizor(sm001)$pc
        polyg = list(original_dat = dat, normalized = temp, smooth001 = sm001)
        
        save(polyg, file = file.path(output, "sim_data", paste0(fname, ".Rdata")))
        write.csv(sm001, file = file.path(output, "smoothed",paste0(fname, "_smoothed_001.csv")), row.names = F, quote = F)
      }
    }
  }
}



# run BayesLASA -----------------------------------------------------------

dir.create(file.path(output, "MCMC"))
ff = dir(file.path(output,"sim_data"), pattern = "Rdata")
tm_df = NULL

generate_gamma <- function(n, k){
  gamma = rep(0, n)
  idx = sample(1:n, k)
  gamma[idx] = 1
  return(gamma)
}


tm_df = foreach(f = ff) %dopar%{
  cat(f)
  load(file.path(output,"sim_data", f))
  dat = polyg$normalized$pc
  dat = dat[-nrow(dat),]
  n = nrow(dat)
  beta_sigma = 1/n
  est.K = round(n/100)
  fold = 100
  if(est.K <3){
    est.K = 3
  }
  gamma_i = generate_gamma(n, est.K)
  a = system.time({r = MCMC_shape(dat,  iter = n*fold,  estK = est.K, gamma_i = gamma_i,
                                  alpha_sigma = 3, beta_sigma = beta_sigma, ppm_store = T)})
  ss = unlist(strsplit(f, "_"))
  num = ss[2]
  equil = ss[4]
  pn = ss[6]
  
  ##ppm
  ppms = r$ppm
  z_ppm <- minbinder(r$ppm, method = "comp")$cl
  L_ppm = which(diff(c(z_ppm[n], z_ppm)) != 0)
  ##ppi
  ppi = r$PPI
  ppis = colSums(ppi)/(r$iter-r$burn)
  #L_ppi = which(ppis >= 0.5)
  L_ppi = which(r$gamma_map > 0)
  L_ppi_cor = NULL
  for(i in 1:ncol(ppi)){
    if(i == ncol(ppi)){
      j = 1
    }else{
      j = i+1
    }
    cr = cor.test(ppi[,i], ppi[,j],  alternative = "less")
    if(!is.na(cr$p.value)){
      if(cr$p.value <= 0.05){
        L_ppi_cor = c(L_ppi_cor, i, j)
      }
    }
  }
  L_ppi_cor = unique(c(L_ppi_cor, L_ppi))
  data = list("L_ppm" = L_ppm,"L_ppi" = L_ppi_cor, "PPI" = ppis)
  fname = gsub(".Rdata", "_MCMC_L.Rdata", f)
  save(data, file = file.path(output, "MCMC", fname))
  
  data.frame(num, equil, pn, elapse = a[3])
  
}
tm_df = do.call("rbind.data.frame", tm_df)
write.csv(tm_df, file = file.path(output,  "MCMC_time_table.csv"))


# ALDUQ -------------------------------------------------------------------

###############
library(R.matlab)
close_point = Vectorize(function(pos, n){
  cp = pos*n+1
  if(cp%%1 >= 0.5){
    cp = ceiling(cp)
  }else{
    cp =  floor(cp)
  }
  if(cp >n) cp = cp-n
  return(cp)
})

dir.create(file.path(output, "JASA_post"))
#**************************************************
#*ALDUQ was run using matlab codes provided by the paper
#**************************************************
matloc="~/Dropbox/shape_analysis_cong/sandbox/output/1225_sim/ALDUQ/"
ff = dir(matloc, pattern = "_JASA_post.mat")
tm_df_jasa = NULL
for(f in ff){
  fname = gsub("_time_.*", "", f)
  ss = unlist(strsplit(f, "_"))
  num = as.numeric(ss[2])
  equil = ss[4]
  pn = as.numeric(ss[6])
  elapse = ss[12]
  if( (!file.exists(file = file.path(output, "JASA_post", file = paste0(fname, "peak_L.rdata"))))){
    
    smdata = gsub("_smoothed.*", "_smoothed_001.csv", fname)
    sm001 = read.csv(file.path(output, "smoothed",smdata))
    pn_t = nrow(sm001)
    
    if(is.null(tm_df_jasa)){
      tm_df_jasa = data.frame(num,  equil, pn, elapse, method = "ALDUQ")
    }else{
      tm_df_jasa = rbind(tm_df_jasa,data.frame(num, equil, pn, elapse, method = "ALDUQ"))
    }
    
    mdat <- readMat(file.path(matloc, f))
    
    k_occur = unlist(lapply(mdat$a, function(x) length(x[[1]])))
    
    k_freq = sort(table(k_occur),decreasing=TRUE)[1]
    k_freq = as.numeric(names(k_freq))
    k_index = which(k_occur == k_freq)
    
    k_post = NULL
    for(i in k_index){
      k_post = cbind(k_post, mdat$a[[i]][[1]])
    }
    Lp = numeric(nrow(k_post))
    L_post = apply(k_post, 2, function(x) close_point(x, pn_t-1))
    #for(k in 1:nrow(k_post)){
    #  ymax = which.max(density(k_post[k,])$y)
    #  Lp[k] = density(k_post[k,])$x[ymax]
    #}
    #L_jasa_post = close_point(Lp, pn_t-1)
    for(k in 2:ncol(L_post)){
      dis_1 = sum(round_distance(L_post[,k], L_post[,k-1], pn_t-1))
      L_shift = c(lead(L_post[,k]), L_post[1,k])
      L_shift = L_shift[!is.na(L_shift)]
      dis_2 = sum(round_distance(L_shift, L_post[,k-1], pn_t-1))
      if(dis_1 > dis_2){
        L_post[,k] = L_shift
      }
    }
    
    L_jasa_post = sort(apply(L_post, 1, function(x) {  as.numeric(names(sort(table(x),decreasing=TRUE)[1])) }))
    save(L_jasa_post, file = file.path(output, "ALDUQ", file = paste0(fname, "peak_L.rdata")))
  }
  
}


for(f in ff){
  fname = gsub("_time_.*", "", f)
  ss = unlist(strsplit(f, "_"))
  num = as.numeric(ss[2])
  equil = ss[4]
  pn = as.numeric(ss[6])
  elapse = ss[12]
  
  if(is.null(tm_df_jasa)){
    tm_df_jasa = data.frame(num,  equil, pn, elapse, method = "ALDUQ")
  }else{
    tm_df_jasa = rbind(tm_df_jasa,data.frame(num, equil, pn, elapse, method = "ALDUQ"))
  }
  
  
}


write.csv(tm_df_jasa, file = file.path(output, "ALDUQ_time_table.csv"), row.names = F, quote = F)



# Convex hull -------------------------------------------------------------

dir.create(file.path(output, "chull"))
tm_df_curv = list()
ff = dir(file.path(output, "sim_data"))

for(f in ff){
  load(file.path(output, "sim_data", f))
  fname = gsub(".Rdata", "", f)
  pc = polyg$original_dat[,c(1, 2)]
  s = unlist(strsplit(fname, "_"))
  num = s[2]
  equil = s[4]
  pn = s[6]
  polyg = list()
  
  a = system.time({hpts <- chull(pc);
  hpts <- c(hpts, hpts[1])})
  
  tm_df_curv[[i]] = data.frame(num,  equil, pn = pn, elapse = a[3], method = "Chull")
  save(hpts, file = file.path(output, "chull", paste0(fname, "_hpts.Rdata")))
  i = i+1
}

tm_df_curv_1 = do.call("rbind", tm_df_curv)
write.csv(tm_df_curv_1, file = file.path(output,  "Chull_time_table.csv"))



# Summary results ---------------------------------------------------------

library(dplyr)
library(stringr)
library(ggpubr)
dir.create(file.path(output, "res_sum"))


jasa_post_loc = file.path(output, "ALDUQ")
ff = dir(jasa_post_loc, pattern = "peak_L.rdata" )

for(f in ff){
  
  fname = gsub("_smoothed.*", "", f)
  ##load chull result
  load(file.path(output, "chull", paste0(fname, "_hpts.Rdata")), verbose = T)
  ##load ADLUQ result
  load(file.path(output, "ALDUQ", f),verbose = T)
  ##load BayesLASA result
  load(file = file.path(output, "MCMC", paste0(fname, "_MCMC_L.Rdata")), verbose = T)
  ##load original data
  load(file = file.path(output, "sim_data", paste0(fname, ".Rdata")), verbose = T)
  
  #load(file.path(output2, "curvature", paste0(fname, "_curvature.Rdata")))
  #res_tb = data.frame(sample = fname, method = c("PPM", "PPI", "JASA001", "JASA005", "JASA01"),
  #                    MCC = numeric(5), ARI = numeric(5))
  res_tb = data.frame(sample = fname, method = c("PPM", "PPM_window5", "PPM_window10", 
                                                 "PPI", "PPI_window5", "PPI_window10",
                                                 "ALDUQ",  "ALDUQ_window5", "ALDUQ_window10",
                                                 "chull", "chull_window5", "chull_window10"),
                      MCC = numeric(12), ARI = numeric(12))
  ####################
  ##get true L and gamma
  ####################
  gamma_true = polyg$original_dat[,3]
  gamma_true = gamma_true[-length(gamma_true)]
  n = length(gamma_true)
  L_true = which(gamma_true != 0)
  
  
  #############
  ###PPM
  #############
  L_ppm = data$L_ppm
  gamma_ppm = rep(0, n)
  gamma_ppm[L_ppm] = 1
  MCC_ppm = L2mcc(gamma_true, L_ppm)
  ARI_ppm = L2adj(gamma_true, L_ppm)
  
  res_tb[res_tb$method == "PPM","MCC"] = MCC_ppm
  res_tb[res_tb$method == "PPM","ARI"] = ARI_ppm
  
  
  L_ppm_5 = L_window(L_true, L_ppm, n, window = 5)
  gamma_ppm_5 = rep(0, n)
  gamma_ppm_5[L_ppm_5] = 1
  MCC_ppm_5 = L2mcc(gamma_true, L_ppm_5)
  ARI_ppm_5 = L2adj(gamma_true, L_ppm_5)
  
  res_tb[res_tb$method == "PPM_window5","MCC"] = MCC_ppm_5
  res_tb[res_tb$method == "PPM_window5","ARI"] = ARI_ppm_5
  
  
  L_ppm_10 = L_window(L_true, L_ppm, n, window = 10)
  gamma_ppm_10 = rep(0, n)
  gamma_ppm_10[L_ppm_10] = 1
  MCC_ppm_10 = L2mcc(gamma_true, L_ppm_10)
  ARI_ppm_10 = L2adj(gamma_true, L_ppm_10)
  
  res_tb[res_tb$method == "PPM_window10","MCC"] = MCC_ppm_10
  res_tb[res_tb$method == "PPM_window10","ARI"] = ARI_ppm_10
  
  ####################
  ###ppi
  ####################
  L_ppi_ci = sort(data$L_ppi)
  
  ppi_df = NULL
  flag = F
  for( i in 1:length(L_ppi_ci)){
    if(flag == FALSE){
      flag = T
      low = L_ppi_ci[i]
    }
    if(flag == T & (i == length(L_ppi_ci) | L_ppi_ci[i] +1 != L_ppi_ci[min(i+1,length(L_ppi_ci)) ] )){
      flag = F
      high = L_ppi_ci[i]
      ppi_df = rbind(ppi_df, c(low, high))
    }
  }
  
  if(ppi_df[1,1] ==1 & ppi_df[nrow(ppi_df), 2] == n){
    ppi_df[1, 1] = ppi_df[nrow(ppi_df), 1]
    ppi_df = ppi_df[-nrow(ppi_df),]
  }
  
  L_ppi = Ci2Bin(ppi_df, L_true)
  MCC_ppi = L2mcc(gamma_true, L_ppi)
  ARI_ppi = L2adj(gamma_true, L_ppi)
  
  
  res_tb[res_tb$method == "PPI","MCC"] = MCC_ppi
  res_tb[res_tb$method == "PPI","ARI"] = ARI_ppi
  
  ######
  
  L_ppi_5 = L_window(L_true,  L_ppi, n, window = 5)
  gamma_ppi_5 = rep(0, n)
  gamma_ppi_5[L_ppi_5] = 1
  MCC_ppi_5 = L2mcc(gamma_true, L_ppi_5)
  ARI_ppi_5 = L2adj(gamma_true, L_ppi_5)
  
  res_tb[res_tb$method == "PPI_window5","MCC"] = MCC_ppi_5
  res_tb[res_tb$method == "PPI_window5","ARI"] = ARI_ppi_5
  
  
  L_ppi_10 = L_window(L_true,  L_ppi, n, window = 10)
  gamma_ppi_10 = rep(0, n)
  gamma_ppi_10[L_ppi_10] = 1
  MCC_ppi_10 = L2mcc(gamma_true, L_ppi_10)
  ARI_ppi_10 = L2adj(gamma_true, L_ppi_10)
  
  res_tb[res_tb$method == "PPI_window10","MCC"] = MCC_ppi_10
  res_tb[res_tb$method == "PPI_window10","ARI"] = ARI_ppi_10
  
  
  ####################
  ########JASA method
  ####################
  
  
  L_jasa_post = sort( L_jasa_post)
  MCC_jasa = L2mcc(gamma_true, L_jasa_post)
  ARI_jasa = L2adj(gamma_true, L_jasa_post)
  
  res_tb[res_tb$method == "ALDUQ","MCC"] = MCC_jasa
  res_tb[res_tb$method == "ALDUQ","ARI"] = ARI_jasa
  
  L_jasa_5 = sort(L_window(L_true, L_jasa_post, n, window = 5))
  gamma_jasa_5 = rep(0, n)
  gamma_jasa_5[L_jasa_5] = 1
  MCC_jasa_5 = L2mcc(gamma_true, L_jasa_5)
  ARI_jasa_5 = L2adj(gamma_true, L_jasa_5)
  
  res_tb[res_tb$method == "ALDUQ_window5","MCC"] = MCC_jasa_5
  res_tb[res_tb$method == "ALDUQ_window5","ARI"] = ARI_jasa_5
  
  
  L_jasa_10 = L_window(L_true, L_jasa_post, n, window = 10)
  gamma_jasa_10 = rep(0, n)
  gamma_jasa_10[L_jasa_10] = 1
  MCC_jasa_10 = L2mcc(gamma_true, L_jasa_10)
  ARI_jasa_10 = L2adj(gamma_true, L_jasa_10)
  
  res_tb[res_tb$method == "ALDUQ_window10","MCC"] = MCC_jasa_10
  res_tb[res_tb$method == "ALDUQ_window10","ARI"] = ARI_jasa_10
  
  
  
  ###############
  ##convex hull
  #############
  
  L_chull = sort(unique(hpts))
  L_chull[L_chull == (n+1)] = 1
  L_chull = sort(L_chull)
  gamma_chull = rep(0, n)
  
  gamma_chull[L_chull] = 1
  MCC_chull = L2mcc(gamma_true, L_chull)
  ARI_chull = L2adj(gamma_true, L_chull)
  
  res_tb[res_tb$method == "chull","MCC"] = MCC_chull
  res_tb[res_tb$method == "chull","ARI"] = ARI_chull
  
  
  L_chull_5 = L_window(L_true, L_chull, n, window = 5)
  gamma_chull_5 = rep(0, n)
  gamma_chull_5[L_chull_5] = 1
  MCC_chull_5 = L2mcc(gamma_true, L_chull_5)
  ARI_chull_5 = L2adj(gamma_true, L_chull_5)
  
  res_tb[res_tb$method == "chull_window5","MCC"] = MCC_chull_5
  res_tb[res_tb$method == "chull_window5","ARI"] = ARI_chull_5
  
  
  L_chull_10 = L_window(L_true, L_chull, n, window = 10)
  gamma_chull_10 = rep(0, n)
  gamma_chull_10[L_chull_10] = 1
  MCC_chull_10 = L2mcc(gamma_true, L_chull_10)
  ARI_chull_10 = L2adj(gamma_true, L_chull_10)
  
  res_tb[res_tb$method == "chull_window10","MCC"] = MCC_chull_10
  res_tb[res_tb$method == "chull_window10","ARI"] = ARI_chull_10
  
  
  write.csv(res_tb, 
            file = file.path(output, "res_sum", paste0(fname, "_summary_coef.csv")),
            quote = F, row.names = F)
}


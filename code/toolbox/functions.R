# ========================================================================================
# get credible interval using ppi
# ========================================================================================

get_cred_interval <- function(L, ppi, len, pp = 0.01) {
  L_cred_interval <- list()
  L <- sort(L)
  for (pos in 1:length(L)) {
    n <- L[pos]
    L_up <- n
    L_lwr <- n
    while ((pos == length(L) & ((n < L[1]) | (n >= L[pos]))) || (n < L[ifelse(pos < length(L), pos + 1, 1)])) {
      l <- n
      n <- n + 1
      if (n > len) {
        n <- 1
      }

      cr <- cor.test(ppi[, l], ppi[, n], alternative = "less")
      if ((!is.na(cr$p.value)) && (cr$p.value <= pp)) {
        L_up <- n
      } else {
        break
      }
    }

    m <- L[pos]
    while ((pos == 1 & ((m >= 1) | (m > L[pos]))) || (m > L[pos - 1])) {
      l <- m
      m <- m - 1
      if (m < 1) {
        m <- len
      }

      cr <- cor.test(ppi[, l], ppi[, m], alternative = "less")
      if ((!is.na(cr$p.value)) && (cr$p.value <= pp)) {
        L_lwr <- m
      } else {
        break
      }
    }

    L_cred_interval <- c(L_cred_interval, list(data.frame(L = L[pos], lwr = L_lwr, upr = L_up)))
  }

  L_cred_interval <- do.call("rbind.data.frame", L_cred_interval)
  L_cred_interval
}


# ========================================================================================
# find the closest point for a given location on a curve
# Inputs
#  pos: position on curve [0, 1]
#  n: total number of vertices in original polygon chain
# Outputs
#  cp: number of vertice to the given position
# ========================================================================================

close_point <- Vectorize(function(pos, n) {
  cp <- pos * n + 1
  if (cp %% 1 >= 0.5) {
    cp <- ceiling(cp)
  } else {
    cp <- floor(cp)
  }
  if (cp > n) cp <- cp - n
  return(cp)
})

# ========================================================================================
# generate random initial gamma for MCMC algorithm
# Inputs
#  n: total number of vertices in original polygon chain
#  k: pre-defiend number of landmark points
# Outputs
#  gamma: binary vector, gamma[i] = 1 indicates vertix[i] is a landmark point while 0 means not.
# ========================================================================================

generate_gamma <- function(n, k) {
  gamma <- rep(0, n)
  idx <- sample(1:n, k)
  gamma[idx] <- 1
  return(gamma)
}


# ========================================================================================
# get the point-wise distance from reduced polygon to original polygon
# Inputs
#  pc: original polygon chain with n 2-dimensional vertices
#  pc_s: reduced polygon chain with n 2-dimensional vertices
#  open: whether the polygon is open,
#        when polygon is closed, di[1] and di[n] are duplicated and only di[1] is kept
# Outputs
#  di: numerical vector of point-wise distance between pc and pc_s, negative value means reduced point
#      on pc_s is outside of pc while positive value means reduced point inside pc
# Created by Cong on 2020/06/23
# use "sf" package to help decided with the point is inside/outside pc
# ========================================================================================
library(sf)
get_di <- function(pc, pc_s, open = F) {
  di <- sqrt(rowSums((pc - pc_s)^2))
  a_poly <- st_polygon(list(pc))
  a <- st_sfc(a_poly)
  p_multi <- st_multipoint(x = pc_s)
  p <- st_cast(st_sfc(p_multi), "POINT")
  ## original polygon in/outside new polygon, 1: inside pc; -1: outside -1
  sign <- ifelse(st_within(p, a, sparse = FALSE)[, 1], 1, -1)
  di <- di * sign
  if (!open) {
    return(di[-length(di)])
  } else {
    return(di)
  }
}

# ========================================================================================
# pc_reducer drops the redudant vertices on a polygon chain
# Inputs
#  pc: a polygon chain with n 2-dimensional vertices
# Outputs
#  The reduced polygon chain
# Created by Qiwei Li on 2020/06/22
# ========================================================================================
pc_reducer <- function(pc) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  closed <- pc_is.closed(pc)
  if (closed) {
    n <- n - 1
  }
  pc <- pc[1:n, ]
  angle <- rep(NA, n)
  for (i in 1:n) {
    if (i == 1) {
      angle[i] <- ta(pc[n, ], pc[i, ], pc[i + 1, ])
    } else if (i == n) {
      angle[i] <- ta(pc[i - 1, ], pc[i, ], pc[1, ])
    } else {
      angle[i] <- ta(pc[i - 1, ], pc[i, ], pc[i + 1, ])
    }
  }
  pc_reduced <- pc[which(angle != pi), ]
  if (closed) {
    pc_reduced <- rbind(pc_reduced, pc_reduced[1, ])
  }
  return(pc_reduced)
}
ta <- function(x, y, z) {
  return(acos(round(sum((x - y) * (z - y)) / sqrt(sum((x - y)^2)) / sqrt(sum((z - y)^2)), 3)))
}

# ========================================================================================
# BayFDR calculates the PPI cutoff that corresponds to a desired false discovery rate
# Inputs
#  PPI: a PPI vector, each element of which is a probability
#  alpha: a desired false discovery rate, with the default value being 0.05
# Outputs
#  The PPI cutoff
# Created by Qiwei Li on 2020/06/22
# ========================================================================================
BayFDR <- function(PPI, alpha = 0.05) {
  PPI_sorted <- sort(PPI, decreasing = TRUE)
  k <- 1
  fdr <- 0
  while (fdr < alpha) {
    fdr <- mean(1 - PPI_sorted[1:k])
    k <- k + 1
  }
  return(PPI_sorted[k])
}

# ========================================================================================
# pc_is.closed checks if a polygon chain is closed or not
# Inputs
#  pc: a polygon chain with n p-dimensional vertices
# Outputs
#  A boolean value
# Created by Qiwei Li on 2020/06/22
# ========================================================================================
pc_is.closed <- function(pc) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  if (sum(pc[1, ] == pc[n, ]) == p && p > 1) {
    closed <- TRUE
  } else {
    closed <- FALSE
  }
  return(closed)
}



# ========================================================================================
# pc_normalizor rescales a polygon chain so that its length equals to a desired value
# Inputs
#  pc: a polygon chain with n p-dimensional vertices
#  length: the desired value of polygon chain length, with the default being one
# Outputs
#  pc: the normalized polygon chain
#  center: the center of the original polygon chain
#  length: the length of the original polygon chain
# Created by Qiwei Li on 2020/06/22
# ========================================================================================
pc_normalizor <- function(pc, length = 1) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  temp <- 0
  for (i in 2:n) {
    temp_2 <- 0
    for (j in 1:p) {
      temp_2 <- temp_2 + (pc[i, j] - pc[i - 1, j])^2
    }
    temp <- temp + sqrt(temp_2)
  }
  temp <- temp / length
  if (pc_is.closed(pc)) {
    # closed
    temp_3 <- colMeans(pc[-n, ])
  } else {
    # open
    temp_3 <- colMeans(pc)
  }
  for (j in 1:p) {
    pc[, j] <- (pc[, j] - temp_3[j]) / temp
  }
  return(list("pc" = pc, "center" = temp_3, "length" = temp))
}



# ========================================================================================
# dist_calculator computes the distances between each adjacent vertices on a polygon chain
# Inputs
#  pc: a polygon chain with n p-dimensional vertices
#  cumsum: a boolean value indicating if cumulative sum of distances is needed
# Outputs
#  dist: the (cumulative) distances between each adjacent vertices
# Created by Qiwei Li on 2020/06/22
# ========================================================================================
dist_calculator_r <- function(pc, cumsum = FALSE) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  dist <- rep(0, n)
  dist[1] <- 0
  for (i in 2:n) {
    temp_2 <- 0
    for (j in 1:p) {
      temp_2 <- temp_2 + (pc[i, j] - pc[i - 1, j])^2
    }
    dist[i] <- sqrt(temp_2)
  }
  if (cumsum) {
    return(cumsum(dist))
  } else {
    return(dist)
  }
}












pc_gsmoother <- function(pc, dist, sigma = 0.05, approx = TRUE) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  pc_temp <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    for (j in 1:p) {
      pc_temp[i, j] <- sum(pc[, j] * wt_calculator(pc = pc, dist = dist, i = i, kernel = "gaussian", sigma = sigma, approx = approx))
    }
  }
  if (pc_is.closed(pc)) {
    pc_temp[n, ] <- pc_temp[1, ]
  }
  return(pc_temp)
}
pc_nsmoother <- function(pc, dist, K = 100, approx = TRUE) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  pc_temp <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    for (j in 1:p) {
      pc_temp[i, j] <- sum(pc[, j] * wt_calculator(pc = pc, dist = dist, i = i, kernel = "nearest", K = K, approx = approx))
    }
  }
  if (pc_is.closed(pc)) {
    pc_temp[n, ] <- pc_temp[1, ]
  }
  return(pc_temp)
}
wt_calculator <- function(pc, dist, ii, kernel = c("gaussian", "nearest"), sigma = 0.1, K = 90, approx = TRUE) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  dist_max <- max(dist)
  closed <- pc_is.closed(pc)
  wt <- rep(NA, n)
  index <- rep(1:n, 2)[ii:(ii + n - 1)]
  count <- 0
  for (i in index) {
    temp <- abs(dist[i] - dist[ii])
    count <- count + 1
    if (closed) {
      temp <- min(temp, dist_max - temp)
    }
    if (approx) {
      if (kernel == "gaussian" && temp > 3 * sigma) {
        break
      }
      if (kernel == "nearest" && count > K) {
        break
      }
    }
    wt[i] <- temp
  }
  if (approx) {
    count <- 0
    for (i in rev(index)) {
      temp <- abs(dist[i] - dist[ii])
      count <- count + 1
      if (closed) {
        temp <- min(temp, dist_max - temp)
      }
      if (kernel == "gaussian" && temp > 3 * sigma) {
        break
      }
      if (kernel == "nearest" && count > K) {
        break
      }
      wt[i] <- temp
    }
  }
  if (closed) {
    wt[n] <- NA
  }
  if (kernel == "gaussian") {
    wt <- -wt^2 / 2 / sigma^2
    wt <- exp(wt - max(wt, na.rm = TRUE))
  } else if (kernel == "nearest") {
    wt <- 1 / rank(wt)
    if (approx) {
      wt[which(wt < 1 / K)] <- NA
    }
    if (closed) {
      wt[n] <- NA
    }
  } else {
    stop(paste0("Check spelling or not support kernel: ", kernel))
  }
  wt <- wt / sum(wt, na.rm = TRUE)
  wt[is.na(wt)] <- 0
  return(wt)
}

srvf_calculator <- function(pc, dist) {
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  closed <- pc_is.closed(pc)
  pc1 <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    pc1[1:(n - 1), j] <- diff(pc[, j]) / diff(dist)
  }
  if (closed) {
    for (j in 1:p) {
      pc1[n, j] <- pc1[1, j]
    }
  } else {
    pc1[n, ] <- NA
  }
  temp <- sqrt(apply(pc1^2, 1, sum))
  for (i in 1:n) {
    if (temp[i] != 0 && !is.na(temp[i])) {
      pc1[i, ] <- pc1[i, ] / temp[i]
    } else {
      pc1[i, ] <- 0
    }
  }
  return(pc1)
}

curv_calculator <- function(pc, dist) {
  # only applicable if p = 1 or p = 2
  n <- dim(pc)[1]
  p <- dim(pc)[2]
  closed <- pc_is.closed(pc)
  curv <- NULL
  if (p == 1) {
    pc1 <- diff(pc) / diff(dist)
    pc2 <- diff(pc1) / diff(dist)[1:(n - 2)]
    curv <- abs(pc2) / (1 + pc1[1:(n - 2)]^2)^(3 / 2)
    curv <- c(curv, NA, NA)
  } else if (p == 2) {
    pc1_x <- diff(pc[, 1]) / diff(dist)
    pc1_y <- diff(pc[, 2]) / diff(dist)
    if (closed) {
      pc1_x <- c(pc1_x, pc1_x[1])
      pc1_y <- c(pc1_y, pc1_y[1])
    } else {
      pc1_x <- c(pc1_x, NA)
      pc1_y <- c(pc1_y, NA)
    }
    pc2_x <- diff(pc1_x) / diff(dist)
    pc2_y <- diff(pc1_y) / diff(dist)
    if (closed) {
      pc2_x <- c(pc2_x, pc2_x[1])
      pc2_y <- c(pc2_y, pc2_y[1])
    } else {
      pc2_x <- c(pc2_x, NA)
      pc2_y <- c(pc2_y, NA)
    }
    curv <- abs(pc1_x * pc2_y - pc1_y * pc2_x) / (pc1_x^2 + pc1_y^2)^(3 / 2)
  } else {
    stop("Only support data with p = 1 or p = 2")
  }
  return(curv)
}

# ========================================================================================
# calculae MCMC
# ========================================================================================
cal_mcc <- function(TP, TN, FP, FN) {
  return((TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
}

# ========================================================================================
# calculae MCMC using landmard points
# ========================================================================================

L2mcc <- function(gamma_true, L_pred) {
  n <- length(gamma_true)
  gamma_pred <- rep(0, n)
  gamma_pred[L_pred] <- 1
  TP <- sum(gamma_pred == 1 & gamma_true == 1)
  TN <- sum(gamma_pred == 0 & gamma_true == 0)
  FP <- sum(gamma_pred == 1 & gamma_true == 0)
  FN <- sum(gamma_pred == 0 & gamma_true == 1)
  return((TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
}

L2adj <- function(gamma_true, L_pred) {
  n <- length(gamma_true)
  gamma_pred <- rep(0, n)
  gamma_pred[L_pred] <- 1
  tab <- table(gamma_true, gamma_pred)
  f2 <- function(n) n * (n - 1) / 2
  sum.f2 <- function(v) sum(f2(v))
  marg.1 <- apply(tab, 1, sum)
  marg.2 <- apply(tab, 2, sum)
  n <- sum(tab)
  prod <- sum.f2(marg.1) * sum.f2(marg.2) / f2(n)
  num <- (sum.f2(as.vector(tab)) - prod)
  den <- 0.5 * (sum.f2(marg.1) + sum.f2(marg.2)) - prod
  return(num / den)
}

# ========================================================================================
# round function
# ========================================================================================
rround <- Vectorize(function(x) {
  dig <- x %% 1
  if (dig < 0.5) {
    return(floor(x))
  } else {
    return(ceiling(x))
  }
})

# ========================================================================================
# get landmark point from CI
# ========================================================================================
Ci2Bin <- function(df, L_true) {
  L_ci <- NULL
  for (i in 1:nrow(df)) {
    interv1 <- c(max(df[i, ]):n, 1:min(df[i, ]))
    interv2 <- min(df[i, ]):max(df[i, ])
    if (df[i, 1] > df[i, 2]) {
      interv <- interv1
    } else {
      interv <- interv2
    }
    if (length(intersect(L_true, interv)) > 0) {
      L_ci <- c(L_ci, intersect(L_true, interv))
    } else {
      L_ci <- c(L_ci, interv[rround(length(interv) / 2)])
    }
  }
  return(sort(unique(L_ci)))
}

# ========================================================================================
# get window around landmark
# ========================================================================================
L_window <- function(L_true, L_pred, n, window = 5) {
  L_pred <- sort(L_pred)
  L_window <- NULL
  for (l in L_pred) {
    ls <- l + -window:window
    ls[ls < 1] <- ls[ls < 1] + n
    ls[ls > n] <- ls[ls > n] - n

    if (length(intersect(L_true, ls)) == 1) {
      L_window <- c(L_window, intersect(L_true, ls))
      L_true <- setdiff(L_true, intersect(L_true, ls))
    } else if (length(intersect(L_true, ls)) == 0) {
      L_window <- c(L_window, l)
    } else {
      ll <- intersect(L_true, ls)
      ll <- ll[which.min(abs(ll - l))]
      L_window <- c(L_window, ll)
      L_true <- setdiff(L_true, ll)
    }
  }

  return(L_window)
}


# TPR=TP/(TP + FP)
L2TPR <- function(gamma_true, L_pred) {
  n <- length(gamma_true)
  gamma_pred <- rep(0, n)
  gamma_pred[L_pred] <- 1
  TP <- sum(gamma_pred == 1 & gamma_true == 1)
  TN <- sum(gamma_pred == 0 & gamma_true == 0)
  FP <- sum(gamma_pred == 1 & gamma_true == 0)
  FN <- sum(gamma_pred == 0 & gamma_true == 1)
  return((TP / (TP + FP)))
}

# FPR=FP/(FP + TN)
L2FPR <- function(gamma_true, L_pred) {
  n <- length(gamma_true)
  gamma_pred <- rep(0, n)
  gamma_pred[L_pred] <- 1
  TP <- sum(gamma_pred == 1 & gamma_true == 1)
  TN <- sum(gamma_pred == 0 & gamma_true == 0)
  FP <- sum(gamma_pred == 1 & gamma_true == 0)
  FN <- sum(gamma_pred == 0 & gamma_true == 1)
  return((FP / (FP + TN)))
}

round_distance <- Vectorize(function(a, b, interval = 1) {
  res1 <- a - b
  if (a < b) {
    res2 <- (a + interval - b)
  } else {
    res2 <- a - b - interval
  }
  return(min(abs(res1), abs(res2)))
})

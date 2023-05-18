##########################
##generate simple polygon
##Cong Zhang
##Feb 10
##add output of gamma 0114
##add pn: total # of random points
##add sigma ranges 0123, sigma on each segment will be sampled from the given sigma vector
##Feb 10 add guassian process
##Mar 22 add jitter to prevent ill_conditioned S matrix
##April 03, add rotation and edge length
##########################

####usage:
### num: number of shape changing point
### equil: if equilateral or not
### sigma: distance of random points to their belong shape segments are from N(0, sigma),sigma can be a scalar or vector
### pn: total num of points = pn*num, # of points on each segments are proportional to seg length
### kernal: normal or gaussian (needs provided l, trim)
### trim: no random points on on the head and tail trim region of each segment
### seed: random seed for generating polygons
###
### output: (num+1)-by-3 matrix, col 1: x-cood, col 2: y-cool, col 3: gamma 
#
#  1. Normal kernel, white noise
# sim_randon_polygon_generator(num = 3, equil = T, sigma = c(0.001, 0.01, 0.5), pn = 100, seed = 1)
#
# 2. Gaussian process
# sim_randon_polygon_generator(3, T, sigma = 1, pn = 300, kernel = "gaussian", l = 0.001, trim = 0.05, seed = 1)



## DEG_MIN < exterior angles < DEG_MAX
#####################
DEG_MAX = 170
DEG_MIN = 15
#####################


#library(sp)
#library(sf)
#library(raster)
#library(lwgeom)
#library(matrixcalc)


#library(foreach)
#library(doParallel)
#doParallel::registerDoParallel(4)

library(dplyr)
#library(mvnfast, package = "mvnfast")
library(mvtnorm)
library(Rcpp)
sourceCpp("code/landmark_detection/MCMC_shape.cpp")
source("code/toolbox/functions.R")

is.square.matrix <- function( x )
{
  ###
  ### determines if the given matrix is a square matrix
  ###
  ### arguments
  ### x = a matrix object
  ###
  if ( !is.matrix( x ) )
    stop( "argument x is not a matrix" )
  return( nrow(x) == ncol(x) )
}

is.symmetric.matrix <- function( x )
{
  ###
  ### this function determines if the matrix is symmetric
  ###
  ### argument
  ### x = a numeric matrix object
  ###
  if ( !is.matrix( x ) ) {
    stop( "argument x is not a matrix" )
  }
  if ( !is.numeric( x ) ) {
    stop( "argument x is not a numeric matrix" )
  }    
  if ( !is.square.matrix( x ) )
    stop( "argument x is not a square numeric matrix" )
  return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
}

is.positive.definite <- function( x, tol=1e-8 )
{
  ###
  ### this function determines if the given real symmetric matrix is positive definite
  ###
  ### parameters
  ### x = a square numeric matrix object
  ### tol = tolerance level for zero
  ###
  if ( !is.square.matrix( x ) )
    stop( "argument x is not a square matrix" )
  if ( !is.symmetric.matrix( x ) )
    stop( "argument x is not a symmetric matrix" )
  if ( !is.numeric( x ) )
    stop( "argument x is not a numeric matrix" )
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow( x )
  for ( i in 1: n ) {
    if ( abs( eigenvalues[i] ) < tol ) {
      eigenvalues[i] <- 0
    }
  }    
  if ( any( eigenvalues <= 0 ) ) {
    return( FALSE )
  }
  return( TRUE )
}

#####dist
### point to points
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
##point p to the line define by b ,c
dist_p2l <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
}

##point p to the line define by b ,c
dist_p2l_direction <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- det(m)/sqrt(sum(v1*v1))
  return(d)
}

##min dist from point p to line set, exclude 0
dist_p2lineset <- function(p, line_set){
  dp <- apply(line_set, 1, function(x) dist_p2l(p,x[1:2], x[3:4]))
  return(min(dp[dp!= 0]))
}


##calculate perimeter of polygon
perimeter <- function(pc){
  n <- dim(pc)[1];
  p <- dim(pc)[2];
  perim = 0
  for(i in 1:n-1){
    perim = perim + euc.dist(pc[i,], pc[i+1,])
  }
  return(perim)
}

####rotation
rotation_vec = function(x,y, theta){
  xp = x*cos(theta)-y*sin(theta)
  yp = x*sin(theta)+y*cos(theta)
  return(c(xp, yp))
}
rotation_shape = function(dat, theta){
  datp = apply(dat, 1, function(x) rotation_vec(x[1], x[2], theta))
  return(t(datp))
}


# intersect ---------------------------------------------------------------


###find intersect point between one segment and a line through a given point with known angle
##point p = (x0, y0)
##angel theta
##let r = (cos theta, sin theta)
##segment q(x1, y1) to q+s (x2, y2), then s  = (x2-x1, y2-y1)
##find t u such that p+tr = q + us

###corss prod for vec
crossvec <- function(a, b){
  return(a[1]*b[2]-a[2]*b[1])
}

####
## aline through p with slope that, intersect with (q, q2)
##return intersect points
intersent_segment_line <- function(p, theta, q, q2){
  r = c(cos(theta), sin(theta))
  s = q2 - q
  if( crossvec(r, s) == 0) return(c(NA, NA))
  u = crossvec((q-p), r)/crossvec(r, s)
  t = crossvec((q-p), s)/crossvec(r, s)
  if(u >= 0 & u <= 1){
    return(q+u*s)
  }else{
    return(c(NA, NA))
  }
}

##return true or false
## if seg(p1, p2) and (q1, q2) intersect
intersent_segments <- function(p1, p2, q1, q2){
  r = p2-p1
  s = q2 - q1
  if( crossvec(r, s) == 0) return(FALSE)
  u = crossvec((q1-p1), r)/crossvec(r, s)
  t = crossvec((q1-p1), s)/crossvec(r, s)
  if(u > 0 & u < 1 & t > 0 & t < 1){
    return(TRUE)
  }else {
    return(FALSE)
  }
}


# Generate simple poly by circle angles ------------------------------------------
theta2corrd <- Vectorize(function(theta,len){
  return(c(cos(theta)*len, sin(theta)*len))
})


##cal angle between 2 vectore
angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

seg_angle <- function(polyg){
  num <- nrow(polyg) -1
  df <- lead(polyg)[1:num,] - polyg[1:num,]
  degree <- apply(cbind(df[1:num,], df[c(2:num, 1),]), 1, function(x) angle(t(x[1:2]), x[3:4]))
  return(degree/pi*180)
}

gen_rand_polyg <- function(num, equil = T){
  if(num < 3) return(NULL)
  degs = 0
  inte = 0
  while(sum(degs <= DEG_MIN |degs >= DEG_MAX) > 0 | inte > 0){
    if(equil){
      theta <- seq(0, 2*pi, length.out = num+1)[1:num]
      len <- 50
    }else{
      theta <- runif(num,0, 2*pi)
      theta <- sort(theta)
      len <- runif(num, 0.1, 1)
      len = round(len*(50/max(len)))
    }
    polyg <- t(theta2corrd(theta, len))  
    polyg <- rbind(polyg, polyg[1,])
    degs <- seg_angle(polyg)
    
    if(equil) inte = 0
    else{
      df <- cbind(polyg[1:num,],lead(polyg)[1:num,])
      inte <- apply(df, 1, function(x) apply(df, 1, function(y) intersent_segments(y[1:2], y[3:4], x[1:2], x[3:4])))
      inte <- sum(inte)
    }
  }
  
  #plot(polyg, type = "n")
  #polygon(polyg)
  #text(polyg[1:num,])
  return(polyg)
}


# add random points along segment -----------------------------------------


####Nearest 2 intersect points
Neast2Intersect <- function(p, theta,line_set) {
  ##find interssct point
  df <- apply(line_set, 1, function(x) intersent_segment_line(p, theta, x[1:2], x[3:4]))
  df <- na.omit(t(df))
  # points(df, pch  = 20, col = "red")
  ##calculate distance
  if(theta == 0){
    side = df[,1] - p[1]
  }else if (theta == pi/2){
    side = df[,2] - p[2]
  }else{
    a <- -1/tan(theta)
    b = p[2]-p[1]*a
    #abline(b, a, col = "blue")
    side <- sign(apply(df, 1, function(x) a*x[1]+b-x[2]))
  }
  
  dt <- apply(df, 1, function(x) euc.dist(x, p))
  ##return nearest two
  up <- which(side >= 0)
  down <- which(side < 0)
  idx.up <- which(dt == min(dt[up]))
  idx.down <- which(dt == min(dt[down]))
  return(c(df[idx.up,], df[idx.down,]))
}

####generate p_new and d from truncated norm
gen_dpts_norm <- function(mean = 0, sigma = 0.01,  p, theta){
  p_new <- p
  d <- rnorm(1, mean, sigma)
  if(theta == pi/2){
    p_new[2] <- p_new[2] + d
  }else if(theta == 0){
    p_new[1] <- p_new[1] + d
  }else{
    dx = d*cos(theta)
    dy = d*sin(theta)
    p_new[1] <- p_new[1] + dx
    p_new[2] <- p_new[2] + dy
  }
  return(p_new)
}

trunct_norm_pts <- function(mean = 0, sigma = 0.01,  p, theta, seg_l){
  if(euc.dist(seg_l[1:2], seg_l[3:4]) < 0.0005) return(p)
  r = -1
  while( (r <= 0) | (r >= 1)){
    p_new <- gen_dpts_norm(mean, sigma, p, theta)
    #points(p_new[1], p_new[2], pch = 20, col = "blue")
    if(theta == 0){
      r = (p_new[1] - seg_l[1])/(seg_l[3]-seg_l[1])
    }else if(theta == 2/pi){
      r = (p_new[2] - seg_l[2])/(seg_l[4]-seg_l[2])
    }else{
      r <- (p_new - seg_l[1:2])/(seg_l[3:4]-seg_l[1:2])
      r <- mean(r, na.rm = T)
    }
  }
  return(p_new)
}

####given p and line set, generate new p perpendicular to original segment while not intersect with rest line set
generate_rand_p <- function(p, line_set, theta, mean = 0, sigma = 0.01){
  seg_l <- Neast2Intersect(p, theta, line_set)
  #points(rbind(seg_l[1:2], seg_l[3: 4]), col =  c("orange", "purple"), pch = 20)
  #cat(sigma)
  p_new <- trunct_norm_pts(mean, sigma,  p, theta, seg_l)
  return(p_new)
}
#p <- rand_p[1,]
#points(p[1], p[2], pch = 19)



add_rand_p_along_seg <- function(num = 10, seg, line_set, sigma = 0.001){
  vs = seg[2,] -seg[1,]
  rand_p = matrix(rep(seg[1,], each = num), ncol = 2) + seq(trim, 1-trim, length.out = num) %*% t(vs)
  #points(rand_p, pch = 20)
  ##get slope and thetaxdiff(seg[,1])
  if(line_slope == 0){
    theta = pi/2
  }else if(line_slope > 10^10){
    theta = 0
  }else{
    m <- -1/line_slope
    theta = atan(m)
  }
  
  pts <- apply(rand_p, 1, function(x) generate_rand_p(x, line_set, theta, mean = 0, sigma))
  return(pts)
}


add_rand_p_along_seg_gaussan = function(num, seg, line_set, sigma, l = l, nu, alpha, trim = 0.05, kernel){
  #cat("num ", num, " \n")
  #cat("sigma ", sigma, " \n")
  #cat("l ", l, " \n")
  
  vs = seg[2,] -seg[1,]
  rand_p = matrix(rep(seg[1,], each = num), ncol = 2) + seq(trim, 1-trim, length.out = num) %*% t(vs)
  
  #points(rand_p, col = "blue", pch = 20)
  line_slope <- diff(seg[,2])/diff(seg[,1])
  if(line_slope == 0){
    theta = pi/2
  }else if(line_slope > 10^10){
    theta = 0
  }else{
    m <- -1/line_slope
    theta = atan(m)
  }
  
  
  if(kernel == "normal"){
    di = rnorm(num, mean = 0, sigma)
    if(theta == 0){
      new_ps = rand_p
      new_ps[,1] =  new_ps[,1] + di
    }else if(theta == 2/pi){
      new_ps = rand_p
      new_ps[,2] =  new_ps[,2] + di
    }else{
      new_ps = rand_p + as.matrix(di) %*% c(cos(theta), sin(theta))
    }
    idx.r = 1
    count = 0
    while(length(idx.r)>0){
      n_pps = rbind(seg[1,], new_ps, seg[2,])
      #points(n_pps, pch = 20, col = "green", cex = 0.5)
      idx.intersct = numeric(num+1)
      for(i in 1:(num+1)){
        idx.intersct[i] = sum(apply(line_set, 1, function(x) intersent_segments(x[1:2], x[3:4], n_pps[i,], n_pps[i+1,])))
      }
      idx.rn = which(idx.intersct != 0)
      if(length(idx.rn) == 0) break;
      
      if( length(idx.rn) == length(idx.r) & all(idx.rn == idx.r)){
        count= count+1
      }else{
        count = 0
      }
      
      if(count > 50){
        new_ps[idx.r,] = rand_p[idx.r,]
      }else{
        idx.rn = c(idx.rn, idx.rn-1, idx.rn+1)
        idx.rn[idx.rn> num] = num
        idx.rn[idx.rn< 1] = 1
        idx.rn = sort(unique(idx.rn))
        idx.r = idx.rn
      }
      if(length(idx.r) >1){
        pts <- apply(rand_p[idx.r,], 1, function(x) generate_rand_p(x, line_set, theta, mean = 0, sigma))
        new_ps[idx.r,] = t(pts)
      }else{
        new_ps[idx.r,] = generate_rand_p(rand_p[idx.r,], line_set, theta, mean = 0, sigma)
      }
      
    }
    #points(new_ps, pch = 20, col = "purple")
    return(new_ps)
    
  }else{
    ##get covariance
    D = compute_dmat(rbind(rbind( seg[1,], rand_p),seg[2,]))
    if(kernel == "gaussian"){
      S = compute_Smat(D, l)
    }else if(kernel == "Matern"){
      S = compute_Matern(D, nu, l)
    }else if(kernel == "RationalQuadratic"){
      S = compute_RationalQuadratic(D, l, alpha)
    }
    
    ##
    #eps <- sqrt(.Machine$double.eps) 
    #S = S+ diag(eps, nrow(S))
    #if(! is.positive.definite(S)) {
    #  stop('change l to making S p.d.!!!')
    #}
    #S_cond = condi_mvn_cov(S, 1:num, c(0, c(num+1)))
    S_cond =symMat2pd(S)[2:(num+1),2:(num+1)]
    intersect_num = 1
    while(intersect_num != 0){
      #cat("intersect with previous line sets\n")
      #cat("sigma=", sigma, " \n")
      #di = rmvn(1, mu = rep(0, num), sigma*S_cond)
      di = rmvnorm(1, mean = rep(0, num), sigma = sigma*S_cond)
      if(theta == 0){
        new_ps = rand_p
        new_ps[,1] =  new_ps[,1] + di
      }else if(theta == 2/pi){
        new_ps = rand_p
        new_ps[,2] =  new_ps[,2] + di
      }else{
        new_ps = rand_p + t(di) %*% c(cos(theta), sin(theta))
      }
      #points(new_ps, pch = 20, col = "red")
      #intersect_num = foreach(i = 1:num, .combine = "+") %do% {
      #  sum(apply(line_set, 1, function(x) intersent_segments(x[1:2], x[3:4], new_ps[i,], rand_p[i,])))
      #}
      
      intersect_num = 0
      n_pps = rbind(seg[1,], new_ps, seg[2,])
      for(i in 1:(num+1)){
        intersect_num = intersect_num + sum(apply(line_set, 1, function(x) intersent_segments(x[1:2], x[3:4], n_pps[i,], n_pps[i+1,])))
      }
    }
    return(new_ps)
  }
  
  
}

add_rand_polyg <- function(polyg, sigma, border, pn = 100, kernel = "normal", l = 0.01,  nu, alpha,trim = 0.05){
  
  polyg <- cbind(polyg, rep(1, nrow(polyg)), rep(0, nrow(polyg)))
  out_tmp <- matrix(polyg[1,], nrow = 1, ncol = 4)
  new_pts <- NULL
  
  perim_l = 0
  for(i in 2:nrow(polyg)){
    perim_l = perim_l + euc.dist(polyg[i-1,], polyg[i,])
  }
  
  #cat("l ", l ," \n")
  #cat("sigma ", sigma ," \n")
  for(i in 1:(nrow(polyg)-1)){
    #cat(i)
    seg_i <- nrow(out_tmp)
    df_tmp <- rbind(out_tmp[,c(1, 2)], polyg[(i+1):nrow(polyg),c(1, 2)])
    sets <- cbind(df_tmp[1: (nrow(df_tmp)-1),],
                  df_tmp[2: nrow(df_tmp),])
    seg <- sets[seg_i,]
    seg <- as.matrix(rbind(seg[1:2], seg[3:4]))
    #segments(seg[1, 1], seg[1, 2], seg[2, 1], seg[2, 2], col = "red")
    rest_line <- sets[-seg_i,]
    rest_line <- rbind(rest_line, border)
    #segments(rest_line[,1], rest_line[,2], rest_line[,3], rest_line[,4])
    ##calculate points on segments (propoertinal to length)
    n = round(euc.dist(seg[1,], seg[2,])/perim_l*(pn-nrow(polyg)))
    if(length(sigma) >1){
      s = sample(sigma, 1)
    }else{
      s = sigma
    }
    #cat("s = ", s, "\n")
    new_pts <- add_rand_p_along_seg_gaussan(num = n, seg, line_set = rest_line, sigma = s, l = l,  nu = nu, alpha = alpha,
                                            trim = trim, kernel = kernel)
    new_pts = cbind(new_pts, rep(0, n), rep(s, n))
    
    out_tmp <- rbind(out_tmp, new_pts, polyg[i+1,])
  }
  #out_tmp[,c(1, 2)] <- scale_perim(out_tmp[,c(1, 2)])
  plot(out_tmp[,c(1, 2)], type = "n", asp = 1)
  polygon(out_tmp[,c(1, 2)])
  colnames(out_tmp) = c("x", "y", "gamma", "sigma")
  rownames(out_tmp) = 1:nrow(out_tmp)
  return(out_tmp)
}

# scale perimeter ---------------------------------------------------------
#*
#scale_perim <- function(polyg){
#  poly.st <- st_polygon(list(polyg))
#  perim <- st_perimeter(poly.st)
#  nz_centroid_sfc = st_centroid(poly.st )
#  nz_scale = (poly.st  - nz_centroid_sfc) * 1/perim 
#  return(nz_scale[[1]])
#}
#

# Gen polygon whole function ----------------------------------------------
##num number of landmark points
##equil
##sigma random points away from segments

sim_randon_polygon_generator <- function(num, equil = T, sigma, pn = 100, edgel = 50, kernel = "normal", l = 1,  nu = Inf, alpha = Inf, trim = 0.05, rotation = F,seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  sim_polyg <- gen_rand_polyg(num, equil)
  if(rotation ){
    theta = runif(1, -1, 1)*pi
    sim_polyg = rotation_shape(sim_polyg, theta)
  }
  sim_polyg = pc_normalizor(sim_polyg, edgel*num)$pc
  #sim_polyg <- scale_perim(sim_polyg)
  #add rand points to segments of polygon
  min_x <- -max(abs(sim_polyg[,1]))-5
  min_y <- -max(abs(sim_polyg[,2]))-5
  max_x <- max(abs(sim_polyg[,1]))+5
  max_y <- max(abs(sim_polyg[,2]))+5
  border <- rbind(c(min_x, min_y, min_x, max_y),
                  c(min_x,min_y, max_x, min_y),
                  c(min_x, max_y, max_x, max_y),
                  c(max_x,min_y, max_x, max_y))
  #plot(sim_polyg, type = "n")
  #polygon(sim_polyg, lty = 2)
  #text(sim_polyg[1:num,])
  #cat("sigma outer ", sigma, "\n")
  set.seed(seed)
  res <- add_rand_polyg(sim_polyg, sigma = sigma, border, pn, kernel = kernel, l = l,  nu = nu, alpha = alpha, trim = trim)
  #add_rand_polyg(polyg= sim_polyg, sigma = sigma, border = border, pn = pn, kernel = kernel, l = l)
  return(res)
}

#
##########
##TEST RUN
##########
#res = sim_randon_polygon_generator(num = 8, equil = F, sigma = 2, pn = 100,
#                                   kernel = "gaussian", l = 1, trim = 0.1, seed = 1)
#res = sim_randon_polygon_generator(3, T, 2, pn = 50, kernel = "gaussian", l = 0.01)
#res = sim_randon_polygon_generator(3, T, 0.0001, pn = 50, kernel = "gaussian", l = 0.0005)
#res = sim_randon_polygon_generator(3, T, c(0.0001, 0.0002), pn = 100, kernel = "gaussian", l = 0.001)
#res = sim_randon_polygon_generator(3, T, 0.0001, pn = 50, kernel = "gaussian", l = 0.0005)
#sim_randon_polygon_generator(3, T, 0.02)
#sim_randon_polygon_generator(4, T, 0.001)
#sim_randon_polygon_generator(4, T, 0.005)
#sim_randon_polygon_generator(4, F, 0.005)

#sim_randon_polygon_generator(4, F, 0.01)


#res = sim_randon_polygon_generator(3, T, c(1, 4), pn = 100)


//////////////
// MCMC
// Cong Zhang
// MCMC for shape analysis
// usage: MCMC_shape(NumericMatrix dat, int iter, int estK = 4, Rcpp::Nullable<Rcpp::IntegerVector> gamma_i = R_NilValue, 
//                      NumericVector updateGp = NumericVector::create(0.8, 0.1, 0.1),
//                        double alpha_sigma = 2, double beta_sigma = 0.01, 
//                        String kernel = "normal", 
//                        bool open = false, 
//                        bool ppm_store = false)
////////////////////////////////
// to get posterior given gamma
// usage: compute_posterior(NumericMatrix dat,  IntegerVector gamma_0, int estK = 4,
//                  double alpha_sigma = 3, double beta_sigma = 0.02)
////////////////////////////////
// added feature:
//    1. check if proposed gamma make polygon intersect
//    2. from L to ppm more quickly
//    3. add option of ppm_store for if calculate and return ppm during MCMC
////////////////////////////////


#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
double dist_point2seg(NumericVector p, NumericVector s1, NumericVector s2);
//NumericVector computeSegLogTerm(int Lk1, int Lk2, NumericVector di);

//compute posterior given dat and gamma
double compute_posterior(NumericMatrix dat,  IntegerVector gamma_0, int estK,
                         double alpha_sigma, double beta_sigma, String kernel, bool open);




/////dist
/// point to points
// [[Rcpp::export]]
double p2pdist(NumericVector s1, NumericVector s2){
  return ( sqrt(pow((s1[0]-s2[0]), 2)+pow((s1[1]-s2[1]), 2)));
}
//add direction
// [[Rcpp::export]]
double dist_point2seg(NumericVector p, NumericVector s1, NumericVector s2)
{
  double a = s2[1] - s1[1]; 
  double b = s1[0] - s2[0];
  double c = s2[0] * s1[1] - s1[0] * s2[1];
  
  return (a * p[0] + b * p[1] + c) / sqrt(a * a + b * b);
}


//calcaulate angel beween three points
// [[Rcpp::export]]
int GetAngleABC( NumericVector a, NumericVector b, NumericVector c )
{
  NumericVector ab = { b[0] - a[0], b[1] - a[1] };
  NumericVector cb = { b[0] - c[0], b[1] - c[1] };
  
  float dot = (ab[0] * cb[0] + ab[1] * cb[1]); // dot product
  float cross = (ab[0] * cb[1] - ab[1] * cb[0]); // cross product
  
  float alpha = atan2(cross, dot);
  
  return (int) floor(alpha * 180. / M_PI + 0.5);
}


//project points onto reduced polygon, (closed)
//return matrix
// [[Rcpp::export]]
NumericMatrix PointOnReducedP(IntegerVector L, const NumericMatrix &dat){
  int n = dat.nrow();
  int k = L.size();
  NumericMatrix new_p(n,2);
  L.push_front(L[k-1]);
  L.push_back(L[1]);
  int j = 0;
  for(int i = 0;i<n;i++){
    if(i>=L[j+1] & j <k){j++;}
    double a = dat(L[j+1],1) - dat(L[j],1); 
    double b = dat(L[j+1],0) - dat(L[j],0);
    double c = dat(L[j+1],0) * dat(L[j],1) - dat(L[j],0) * dat(L[j+1],1);
    new_p(i,0) = ( b*b*dat(i,0)+a*b*dat(i,1)-a*c)/(a*a+b*b);
    new_p(i,1) = ( a*a*dat(i,1)+a*b*dat(i,0)+b*c)/(a*a+b*b);
  }
  
  return(new_p);
}


// [[Rcpp::export]]
IntegerVector gamma2L(IntegerVector gamma){
  int n = gamma.size();
  //get landmark points
  IntegerVector L;
  for(int i=0; i<n; i++){
    if(gamma[i] == 1) L.push_back(i);
  }
  return L;
}

// [[Rcpp::export]]
IntegerVector gamma2Z(IntegerVector gamma){
  int n = gamma.size();
  int K = sum(gamma);
  //get landmark points
  IntegerVector Z = cumsum(gamma);
  for(int i=0; i<n; i++){
    if(Z[i] == 0){
      Z[i] = K;
    }else{
      break;
    }
  }
  return Z;
}

// [[Rcpp::export]]
IntegerMatrix gamma2seg(IntegerVector gamma, bool open = false){
  IntegerVector L = gamma2L(gamma);
  int n_L = L.size();
  if(open){
    L.push_back(gamma.size());
    IntegerMatrix seg(n_L,2);
    seg(_,0) = L[seq(0, n_L-1)];
    seg(_,1) = L[seq(1, n_L)];
    return seg;
  }else{
    IntegerMatrix seg(n_L,2);
    seg(_,0) = L;
    L.push_back(L[0]);
    seg(_,1) = L[seq(1, n_L)];
    return seg;
  }
}

//subset rows of matrix by condition
// [[Rcpp::export]]
IntegerMatrix submat_condition(IntegerMatrix X, LogicalVector condition) { 
  int n=X.nrow(), k=X.ncol();
  IntegerMatrix out(sum(condition),k);
  for (int i = 0, j = 0; i < n; i++) {
    if(condition[i]) {
      out(j,_) = X(i,_);
      j = j+1;
    }
  }
  return(out);
}

//subset rows of matrix by condition
// [[Rcpp::export]]
NumericMatrix submat_condition_n(NumericMatrix X, LogicalVector condition) { 
  int n=X.nrow(), k=X.ncol();
  NumericMatrix out(sum(condition),k);
  for (int i = 0, j = 0; i < n; i++) {
    if(condition[i]) {
      out(j,_) = X(i,_);
      j = j+1;
    }
  }
  return(out);
}

// [[Rcpp::export]]
Rcpp::List diff_seg(IntegerVector gamma1, IntegerVector gamma2, bool open = false){
  IntegerMatrix seg1 = gamma2seg(gamma1, open = open);
  IntegerMatrix seg2 = gamma2seg(gamma2, open = open);
  int n1 = seg1.nrow();
  int n2 = seg2.nrow();
  IntegerVector ssame1 = rep(0, n1);
  IntegerVector ssame2 = rep(0, n2);
  for(int i =0; i<n1;i++){
    //Rcout<<i<< "  i   \n";
    for(int j =0; j<n2;j++){
      //Rcout<<j<< "  j   \n";
      if(seg1(i, 0) == seg2(j, 0) & seg1(i, 1) == seg2(j, 1)){
        ssame1[i] = 1;
        ssame2[j] = 1;
        break;
      }
    }
  }
  IntegerMatrix d = submat_condition(seg1, ssame1 == 0);
  IntegerMatrix a = submat_condition(seg2, ssame2 == 0);
  List L = List::create(Named("deleted") = d , _["added"] = a);
  return L;
}




// [[Rcpp::export]]
double computeSegLogTerm(int Lk1, int Lk2, NumericMatrix dat, double alpha_sigma, double beta_sigma, String kernel = "normal", double phi = 1, bool fixedSigma = false,
                         bool open = false){
  int n = dat.nrow();
  int n_k;
  if(Lk1<Lk2){
    n_k = Lk2 - Lk1;
  }else{
    n_k = n - Lk1 + Lk2;
  }
  //Rcout<< n_k << " n_k\n";
  //Rcout << Lk1  << " Lk1\n";
  //Rcout << dat(Lk1,0) <<  " "<<dat(Lk1,1) << " dat_lk1\n";
  //Rcout << dat(Lk2,0)  <<  " "<<dat(Lk2,1) << " dat_lk2\n";
  NumericVector di(n_k);
  int j = 0;
  if(Lk1<Lk2){
    if(open & (Lk2 ==n)){
      for(int i = Lk1; i<Lk2;i++){
        di[j] = dist_point2seg(dat(i,_), dat(Lk1,_), dat(Lk2-1,_));
        j++;
      }
    }else{
      for(int i = Lk1; i<Lk2;i++){
        di[j] = dist_point2seg(dat(i,_), dat(Lk1,_), dat(Lk2,_));
        j++;
      }
    }
    
  }else{
    for(int i = 0; i<Lk2;i++){
      di[j] = dist_point2seg(dat(i,_), dat(Lk1,_), dat(Lk2,_));
      j++;
    }
    for(int i = Lk1; i<n;i++){
      di[j] = dist_point2seg(dat(i,_), dat(Lk1,_), dat(Lk2,_));
      j++;
    }
  }
  double a_term, b_term, term;
  if(fixedSigma){
    b_term = sum(pow(di, 2.0))/(2.0*pow(phi, 2) );
    term = -(n_k)/2.0*log(2*M_PI) -(n_k)*log(phi) - b_term;
  }else if(kernel == "normal"){
    a_term = alpha_sigma+n_k/2.0;
    b_term = beta_sigma+sum(pow(di, 2.0))/2.0;
    term = (lgamma(a_term)-a_term*log(b_term));
  }else if(kernel == "laplace"){
    a_term = alpha_sigma+n_k;
    b_term = beta_sigma+sqrt(2)*sum(abs(di));
    term = (lgamma(a_term)-a_term*log(b_term));
  }
  //Rcout<< a_term << " a_term\n";
  //Rcout<< b_term << " b_term\n";
  return term;
}



// [[Rcpp::export]]
Rcpp::List updategamma(IntegerVector gamma, NumericVector p = NumericVector::create(0.8, 0.1, 0.1)) {
  //Rcout<< p << "  p   \n";
  //int flag = sample(LogicalVector::create(true, false), 1, false, NumericVector::create(p, (1-p)/2.0, (1-p)/2.0 ))[0];
  int flag = sample(IntegerVector::create(1, 2, 3), 1, false, p)[0];
  //Rcout<< flag_add_delete << "  flag   \n";
  int cnt = 0 ;
  int n = gamma.size();
  IntegerVector newgamma = clone(gamma);
  IntegerVector idx = seq_len(n);
  // add/delete
  if((flag == 1) | (sum(gamma) == n)){
    while(cnt < 3){
      int candidate = sample(idx, 1, false)[0] ;
      newgamma[candidate-1] = 1- newgamma[candidate-1];
      cnt = sum(newgamma);
    }
  // swap
  }else if(flag == 2){
    IntegerVector L = gamma2L(gamma);
    int p1 = sample(L, 1, false)[0];
    int p0;
    int c1, c2;
    int i = 1;
    //Rcout << "p1 " << p1 << " \n";
    while(i < n){
      c1 = ( (p1+i < n)? (p1+i): (i-1));
      c2 = ( (p1-i >=0)? (p1-i): (p1+n-i));
      //Rcout << "i " << i << " \n";
      //Rcout << "c1 " << c1 << " \n";
      //Rcout << "c2 " << c2 << " \n";
      if( (gamma[c1] == 1) &  (gamma[c2] == 1)){
        i++;
        continue;
      }
      if(gamma[c1] == 1){
        p0 = c2;
        break;
      }else if(gamma[c2] == 1){
        p0 = c1;
        break;
      }else{
        p0 = sample(IntegerVector::create(c1, c2), 1, false)[0];
        break;
      }
    }
    //Rcout << "p0 " << p0 << " \n";
    newgamma[p1] = 1- newgamma[p1];
    newgamma[p0] = 1- newgamma[p0];
    //shift
  }else if(flag == 3){
    bool left = sample(LogicalVector::create(true, false), 1, false)[0];
    if(left){
      //Rcout << "left" << left << "\n";
      gamma.push_back(gamma[0]);
      newgamma = gamma[seq(1, n)];
    }else{
      //Rcout << "right" << left << "\n";
      gamma.push_front(gamma[n-1]);
      newgamma = gamma[seq(0, n-1)];
    }
  }
  List L = List::create(Named("gamma") = newgamma , _["flag"] = flag);
  return L ;
}


//////////////////////
//check if proposed gamma make polygon intersect
//corss prod for vec
// [[Rcpp::export]]
double crossvec_rcpp(NumericVector a, NumericVector b){
  return(a[0]*b[1]-a[1]*b[0]);
}

//return true or false
// if seg(p1, p2) and (q1, q2) intersect
// [[Rcpp::export]]
bool intersect_segments_rcpp(NumericVector p1, NumericVector p2, NumericVector q1, NumericVector q2){
  NumericVector r = { p2[0] - p1[0], p2[1] - p1[1] };
  NumericVector s = { q2[0] - q1[0], q2[1] - q1[1] };
  
  if( crossvec_rcpp(r, s) == 0) return(false);
  double u = crossvec_rcpp((q1-p1), r)/crossvec_rcpp(r, s);
  double t = crossvec_rcpp((q1-p1), s)/crossvec_rcpp(r, s);
  if(u > 0 & u < 1 & t > 0 & t < 1){
    return(true);
  }else {
    return(false);
  }
}


//return true or false
// check if there is intersected segment in polygonal chain
// [[Rcpp::export]]
bool poly_intersect_segments_rcpp(NumericMatrix &dat){
  int n = dat.nrow();
  for(int i =0;i<n-2;i++){
    for(int j = 0;j<n-1;j++){
      if(intersect_segments_rcpp(dat(i,_), dat(i+1,_), dat(j,_), dat(j+1,_))){
        return(true);
      }
    }
  }
  return(false);
}

//return true or false
// check if there is intersected segment in polygonal chain
// [[Rcpp::export]]
bool intersent_shape(NumericMatrix &dat, IntegerVector gamma, bool open = false){
  
  NumericMatrix datr = submat_condition_n(dat, gamma == 1);
  int n = datr.nrow();
  for(int i =0;i<n-2;i++){
    for(int j = 0;j<n-1;j++){
      if(intersect_segments_rcpp(datr(i,_), datr(i+1,_), datr(j,_), datr(j+1,_))){
        return(true);
      }
    }
  }
  if(! open){
    for(int i =0;i<n-1;i++){
      if(intersect_segments_rcpp(datr(i,_), datr(i+1,_), datr(n-1,_), datr(0,_))){
        return(true);
      }
    }
  }
  return(false);
}




// [[Rcpp::export]]
NumericVector cal_di(NumericMatrix dat, IntegerVector L){
  int n = dat.nrow();
  IntegerVector L2 = clone(L);
  L2.push_front(L2[L2.size()-1]);
  L2.push_back(L2[1]);
  //Rcout << "L2 " << L2 << "\n\n";
  int j = 0;
  NumericVector d(n);
  for(int i =0; i<n;i++){
    //Rcout << "i " << i << "\n";
    if( (i >= L2[j+1]) & (j <L2.size() -2)){
      j++;
    }
    //Rcout << "j " << j << "\n";
    //Rcout << "L " << L2[j] << " "  << L2[j+1] << "\n";
    d[i] = dist_point2seg(dat(i,_), dat(L2[j], _), dat(L2[j+1], _));
    //Rcout << "d" << d << "\n";
  }
  return(d);
}


// [[Rcpp::export]]
IntegerVector L2gamma(IntegerVector L, int n){
  IntegerVector gamma(n);
  for(int i = 0;i < L.size();i++){
    gamma[L[i]] = 1;
  }
  return gamma;
}

// [[Rcpp::export]]
void update_ppm_count(NumericMatrix &ppm, IntegerVector gamma){
  IntegerVector Z = gamma2Z(gamma);
  int n = gamma.size();
  for(int i = 0; i < (n-1);i++){
    //Rcout << "i " << i << " \n";
    for(int j = i+1;j < n; j++){
      //Rcout << "j " << j << " \n";
      if(Z[i] == Z[j]){
        ppm(i,j)++;
        ppm(j,i)++;
      }
    }
  }
}


// [[Rcpp::export]]
void update_ppm_count2(mat &ppm, IntegerVector L){
  int K = L.size();
  int n = ppm.n_cols;
  for(int i = 0;i<K-1;i++){
    ppm.submat(L[i], L[i], L[i+1]-1, L[i+1]-1) = ppm.submat(L[i], L[i], L[i+1]-1, L[i+1]-1)+1;
  }
  ivec a;
  if(L[0] != 0){
    a = seq(0, L[0]-1);
  }
  ivec b = seq(L[K-1], n-1);
  uvec idx =  arma::conv_to<arma::uvec>::from(join_cols(a,b));
  ppm.submat(idx, idx) = ppm.submat(idx, idx)+1;
}


// [[Rcpp::export]]
double LogOmegaterm(int n, int K, double alpha, double beta){
  return( lgamma(alpha +K) + lgamma(beta+n-K));
}

// [[Rcpp::export]]
Rcpp::List MCMC_shape(NumericMatrix dat, int iter, int estK = 4, 
                      Rcpp::Nullable<Rcpp::IntegerVector> gamma_i = R_NilValue, 
                      NumericVector updateGp = NumericVector::create(0.8, 0.1, 0.1),
                      double alpha_sigma = 2, double beta_sigma = 0.01,
                      String kernel = "normal", bool open = false,
                      bool ppm_store = false) {
  
  int burn = iter/2;
  int n = dat.nrow();
  
  //hyper-parameter 
  double alpha_w = 2.0*estK/(n-1.0);
  double beta_w = 2.0 - alpha_w;
  
  //Initialize gamma
  IntegerVector gamma_0(n);
  if(gamma_i.isNotNull()){
    gamma_0 = gamma_i;
  }else{
    gamma_0 = rep(0, n);
    IntegerVector idx = seq(0, n-1);
    IntegerVector p1 = sample(idx, estK, false) ;
    gamma_0[p1] = 1;
    
    while(intersent_shape(dat, gamma_0, open)){
      gamma_0 = rep(0, n);
      IntegerVector idx = seq(0, n-1);
      IntegerVector p1 = sample(idx, estK, false) ;
      gamma_0[p1] = 1;
    }
    
    
  }
  if(open){  //first point shall be landmark points, last point should be included in the last segment
    gamma_0[0] = 1;
    gamma_0[n-1] = 0;
    //gamma_0.push_back(gamma_0[n-1]);
    //gamma_0[n] = 1;
  }
  
  
  //get used variables
  IntegerVector L_0 = gamma2L(gamma_0);
  int K_0 = L_0.size();
  //if(open){
  //  K_0 = K_0 -1;
  //}
  
  
  //hyper-parameter terms
  double log_K_term = alpha_sigma*log(beta_sigma) - lgamma(alpha_sigma);
  //double Log_Gamma_term = LogOmegaterm(n, K_0, alpha_w, beta_w);
  double Log_Gabw_gamma1_term = lgamma(alpha_w +1) + lgamma(beta_w +1 -1);
  double Log_Gabw_gamma0_term = lgamma(alpha_w +0) + lgamma(beta_w +1 -0);
  
  //get segment
  IntegerMatrix df_0 = gamma2seg(gamma_0, open = open);
  NumericVector df_0_terms(K_0);
  for(int i = 0;i<K_0;i++){
    df_0_terms[i] = computeSegLogTerm(df_0(i,0), df_0(i,1), dat, alpha_sigma, beta_sigma, kernel);
  }
  
  double post_0 = K_0*log_K_term + sum(df_0_terms) +(K_0)*Log_Gabw_gamma1_term + (n-K_0)*Log_Gabw_gamma0_term;
  //double post_0 = K_0*log_K_term + sum(df_0_terms) +Log_Gamma_term;
  //record intermediate values
  NumericVector hastings(iter);
  double posterior_max = 0;
  int gamma_map_index = 0;
  Rcpp::List Llist;
  IntegerVector Ks(iter);
  //NumericVector PPI (n);
  NumericMatrix PPI(iter-burn+1,n);
  IntegerVector gamma_map (n);
  //ppm
  mat ppm(n, n, fill::zeros);
  
  //
  int count = 0;
  double accept_ad = 0;
  double accept_swap = 0;
  double accept_shift = 0;
  double gamma_ad = 0;
  double gamma_swap = 0;
  double gamma_shift = 0;
  int larger_lklh = 0;
  
  for(int ii = 0; ii < iter; ii++){
    ///////////////////
    //propose new gamma
    ///////////////////
    Rcpp::List newgamma = updategamma(gamma_0,updateGp);
    IntegerVector gamma_tmp = newgamma["gamma"];
    
    while(intersent_shape(dat, gamma_tmp, open)){
      newgamma = updategamma(gamma_0,updateGp);
      gamma_tmp = newgamma["gamma"];
    }
    
    int gamma_flag = newgamma["flag"];
    switch(gamma_flag){
    case 1:
      gamma_ad++;
      break;
    case 2:
      gamma_swap++;
      break;
    case 3:
      gamma_shift++;
      break;
    }
    ////////////////////////
    
    if(open){  //first and added last one to be landmark points
      gamma_tmp[0] = 1;
      gamma_tmp[n-1] = 0;
    }
    
    //get used variables
    IntegerVector L_tmp = gamma2L(gamma_tmp);
    
    int K_tmp = L_tmp.size();
    //if(open){
    //  K_tmp = K_tmp -1;
    //}
    
    //IntegerVector z_tmp = cumsum(gamma_tmp);
    //z_tmp[z_tmp==0] <- K_tmp;
    
    
    //changed segments
    //IntegerMatrix df_tmp = gamma2seg(gamma_tmp);
    Rcpp::List seg_dif = diff_seg(gamma_0, gamma_tmp, open = open);
    IntegerMatrix add_seg_terms = seg_dif["added"];
    IntegerMatrix del_seg_terms = seg_dif["deleted"];
    
    NumericVector added_terms(add_seg_terms.nrow());
    for(int i = 0;i < add_seg_terms.nrow();i++){
      added_terms[i] = computeSegLogTerm(add_seg_terms(i,0), add_seg_terms(i,1), dat, alpha_sigma, beta_sigma, kernel);
    }
    
    NumericVector del_terms(del_seg_terms.nrow());
    for(int i = 0;i < del_seg_terms.nrow();i++){
      del_terms[i] = computeSegLogTerm(del_seg_terms(i,0), del_seg_terms(i,1), dat, alpha_sigma, beta_sigma, kernel);
    }
    
    //gamma omega
    double Log_Gamma_term_tmp = LogOmegaterm(n, K_tmp, alpha_w, beta_w);
    
    double hasting_diff = (K_tmp - K_0)*log_K_term + sum(added_terms) - sum(del_terms) + 
      (K_tmp - K_0)*Log_Gabw_gamma1_term + (K_0-K_tmp)*Log_Gabw_gamma0_term;
    //double hasting_diff = (K_tmp - K_0)*log_K_term + sum(added_terms) - sum(del_terms) + 
    //  Log_Gamma_term_tmp - Log_Gamma_term;
    
    if(hasting_diff > 0)  larger_lklh =  larger_lklh +1;
    
    double pa = 1;
    if(exp(hasting_diff) < 1) pa = exp(hasting_diff);
    bool accept = sample(LogicalVector::create(true, false), 1, false, NumericVector::create(pa, 1-pa))[0];
    if(accept){
      gamma_0 = gamma_tmp;
      L_0 = L_tmp;
      K_0 = K_tmp;
      hastings[ii] = hasting_diff;
      switch(gamma_flag){
      case 1:
        accept_ad++;
        break;
      case 2:
        accept_swap++;
        break;
      case 3:
        accept_shift++;
        break;
      }
    }
    
    if( sum(hastings) > posterior_max) {
      posterior_max = sum(hastings);
      gamma_map_index = ii;
      gamma_map = clone(gamma_tmp);
    }
    
    Llist.push_back(L_0);
    Ks[ii] = K_0;
    
    if((ii >= burn) & ppm_store){
      //PPI = PPI + as<NumericVector>(gamma_0);
      PPI(ii-burn,_) = gamma_0;
      update_ppm_count2(ppm, L_0);
    }
    
    // Monitor the process
    if(ii*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
  }
  
  double accept_r_ad = accept_ad/gamma_ad;
  double accept_r_swap = accept_swap/gamma_swap;
  double accept_r_shift = accept_shift/gamma_shift;
  NumericVector posteriors = cumsum(hastings);
  posteriors = posteriors + post_0;
  
  //PPI = PPI/(iter-burn+1);
  if(ppm_store){
    ppm = ppm/(iter-burn);
    return(Rcpp::List::create(Named("iter") = iter, _["burn"] = burn, _["gamma_map"] = gamma_map,
                              _["gamma_map_index"]= gamma_map_index,  _["Llist"] = Llist,
                              _["hastings"] = hastings, _["posteriors"] = posteriors, _["ppm"] = ppm,_["PPI"] = PPI,
                                _["Ks"] = Ks, _["accept_rate_ad"] = accept_r_ad, _["accept_rate_swap"] = accept_r_swap,
                                  _["accept_rate_shift"] = accept_r_shift, _["larger_posterior"] = larger_lklh));
  }else{
    return(Rcpp::List::create(Named("iter") = iter, _["burn"] = burn, _["gamma_map"] = gamma_map,
                              _["gamma_map_index"]= gamma_map_index,  _["Llist"] = Llist,
                              _["hastings"] = hastings, _["posteriors"] = posteriors, 
                              _["Ks"] = Ks, _["accept_rate_ad"] = accept_r_ad, _["accept_rate_swap"] = accept_r_swap,
                                _["accept_rate_shift"] = accept_r_shift, _["larger_posterior"] = larger_lklh));
  }
  
}



// [[Rcpp::export]]
Rcpp::List MCMC_shape_fixedSigma(NumericMatrix dat, int iter, int estK = 4, 
                                 Rcpp::Nullable<Rcpp::IntegerVector> gamma_i = R_NilValue, 
                                 NumericVector updateGp = NumericVector::create(0.8, 0.1, 0.1),
                                 double alpha_sigma = 2, double beta_sigma = 0.01, String kernel = "normal", bool open = false, bool ppm_store = false,
                                 double phi = 1, bool fixedSigma = false) {
  
  int burn = iter/2;
  int n = dat.nrow();
  
  //hyper-parameter 
  double alpha_w = 2.0*estK/(n-1.0);
  double beta_w = 2.0 - alpha_w;
  
  
  
  //Initialize gamma
  IntegerVector gamma_0(n);
  if(gamma_i.isNotNull()){
    gamma_0 = gamma_i;
  }else{
    gamma_0 = rep(0, n);
    IntegerVector idx = seq(0, n-1);
    IntegerVector p1 = sample(idx, estK, false) ;
    gamma_0[p1] = 1;
    
    while(intersent_shape(dat, gamma_0, open)){
      gamma_0 = rep(0, n);
      IntegerVector idx = seq(0, n-1);
      IntegerVector p1 = sample(idx, estK, false) ;
      gamma_0[p1] = 1;
    }
    
    
  }
  if(open){  //first point shall be landmark points, last point should be included in the last segment
    gamma_0[0] = 1;
    gamma_0[n-1] = 0;
    //gamma_0.push_back(gamma_0[n-1]);
    //gamma_0[n] = 1;
  }
  
  
  //get used variables
  IntegerVector L_0 = gamma2L(gamma_0);
  int K_0 = L_0.size();
  //if(open){
  //  K_0 = K_0 -1;
  // }
  
  //hyper-parameter terms
  double log_K_term = alpha_sigma*log(beta_sigma) - lgamma(alpha_sigma);
  if(fixedSigma){
    log_K_term = 0;
  }
  //double Log_Gamma_term = LogOmegaterm(n, K_0, alpha_w, beta_w);
  
  double Log_Gabw_gamma1_term = lgamma(alpha_w +1) + lgamma(beta_w +1 -1);
  double Log_Gabw_gamma0_term = lgamma(alpha_w +0) + lgamma(beta_w +1 -0);
  
  //Rcout << "GET SEGMENTS\n";
  //get segment
  IntegerMatrix df_0 = gamma2seg(gamma_0, open = open);
  //Rcout << "GET SEGMENTS DONE\n";
  //Rcout << df_0 << std::endl;
  NumericVector df_0_terms(K_0);
  for(int i = 0;i<K_0;i++){
    df_0_terms[i] = computeSegLogTerm(df_0(i,0), df_0(i,1), dat, alpha_sigma, beta_sigma, kernel, phi, fixedSigma, open);
  }
  //Rcout << "LKLH DONE\n";
  //double post_0 = K_0*log_K_term + sum(df_0_terms) + Log_Gamma_term;
  
  double post_0 = K_0*log_K_term + sum(df_0_terms) +(K_0)*Log_Gabw_gamma1_term + (n-K_0)*Log_Gabw_gamma0_term;
  
  //record intermediate values
  NumericVector hastings(iter);
  double posterior_max = 0;
  int gamma_map_index = 0;
  Rcpp::List Llist;
  IntegerVector Ks(iter);
  //NumericVector PPI (n);
  NumericMatrix PPI(iter-burn+1,n);
  IntegerVector gamma_map (n);
  //ppm
  mat ppm(n, n, fill::zeros);
  
  //
  int count = 0;
  double accept_ad = 0;
  double accept_swap = 0;
  double accept_shift = 0;
  double gamma_ad = 0;
  double gamma_swap = 0;
  double gamma_shift = 0;
  int larger_lklh = 0;
  
  //Rcout << "START MCMC RUN!!\n";
  
  for(int ii = 0; ii < iter; ii++){
    ///////////////////
    //propose new gamma
    ///////////////////
    Rcpp::List newgamma = updategamma(gamma_0,updateGp);
    IntegerVector gamma_tmp = newgamma["gamma"];
    
    while(intersent_shape(dat, gamma_tmp, open)){
      newgamma = updategamma(gamma_0,updateGp);
      gamma_tmp = newgamma["gamma"];
    }
    
    int gamma_flag = newgamma["flag"];
    switch(gamma_flag){
    case 1:
      gamma_ad++;
      break;
    case 2:
      gamma_swap++;
      break;
    case 3:
      gamma_shift++;
      break;
    }
    ////////////////////////
    
    if(open){  //first and added last one to be landmark points
      gamma_tmp[0] = 1;
      gamma_tmp[n-1] = 0;
    }
    
    //get used variables
    IntegerVector L_tmp = gamma2L(gamma_tmp);
    
    int K_tmp = L_tmp.size();
    //if(open){
    //  K_tmp = K_tmp -1;
    //}
    
    //IntegerVector z_tmp = cumsum(gamma_tmp);
    //z_tmp[z_tmp==0] <- K_tmp;
    
    
    //changed segments
    //IntegerMatrix df_tmp = gamma2seg(gamma_tmp);
    Rcpp::List seg_dif = diff_seg(gamma_0, gamma_tmp, open = open);
    IntegerMatrix add_seg_terms = seg_dif["added"];
    IntegerMatrix del_seg_terms = seg_dif["deleted"];
    
    NumericVector added_terms(add_seg_terms.nrow());
    for(int i = 0;i < add_seg_terms.nrow();i++){
      added_terms[i] = computeSegLogTerm(add_seg_terms(i,0), add_seg_terms(i,1), dat, alpha_sigma, beta_sigma, kernel, phi, fixedSigma, open);
    }
    
    NumericVector del_terms(del_seg_terms.nrow());
    for(int i = 0;i < del_seg_terms.nrow();i++){
      del_terms[i] = computeSegLogTerm(del_seg_terms(i,0), del_seg_terms(i,1), dat, alpha_sigma, beta_sigma, kernel, phi, fixedSigma, open);
    }
    
    //double Log_Gamma_term_tmp = LogOmegaterm(n, K_tmp, alpha_w, beta_w);
    //double hasting_diff = (K_tmp - K_0)*log_K_term + sum(added_terms) - sum(del_terms) + 
    //  Log_Gamma_term_tmp - Log_Gamma_term;
    
    double hasting_diff = (K_tmp - K_0)*log_K_term + sum(added_terms) - sum(del_terms) + 
      (K_tmp - K_0)*Log_Gabw_gamma1_term + (K_0-K_tmp)*Log_Gabw_gamma0_term;
    
    
    if(hasting_diff > 0)  larger_lklh =  larger_lklh +1;
    
    double pa = 1;
    if(exp(hasting_diff) < 1) pa = exp(hasting_diff);
    bool accept = sample(LogicalVector::create(true, false), 1, false, NumericVector::create(pa, 1-pa))[0];
    if(accept){
      gamma_0 = gamma_tmp;
      L_0 = L_tmp;
      K_0 = K_tmp;
      hastings[ii] = hasting_diff;
      switch(gamma_flag){
      case 1:
        accept_ad++;
        break;
      case 2:
        accept_swap++;
        break;
      case 3:
        accept_shift++;
        break;
      }
    }
    
    if( sum(hastings) > posterior_max) {
      posterior_max = sum(hastings);
      gamma_map_index = ii;
      gamma_map = clone(gamma_tmp);
    }
    
    Llist.push_back(L_0);
    Ks[ii] = K_0;
    
    if((ii >= burn) & ppm_store){
      //PPI = PPI + as<NumericVector>(gamma_0);
      PPI(ii-burn,_) = gamma_0;
      update_ppm_count2(ppm, L_0);
    }
    
    // Monitor the process
    if(ii*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
  }
  
  double accept_r_ad = accept_ad/gamma_ad;
  double accept_r_swap = accept_swap/gamma_swap;
  double accept_r_shift = accept_shift/gamma_shift;
  NumericVector posteriors = cumsum(hastings);
  posteriors = posteriors + post_0;
  
  //PPI = PPI/(iter-burn+1);
  if(ppm_store){
    ppm = ppm/(iter-burn);
    return(Rcpp::List::create(Named("iter") = iter, _["burn"] = burn, _["gamma_map"] = gamma_map,
                              _["gamma_map_index"]= gamma_map_index,  _["Llist"] = Llist,
                              _["hastings"] = hastings, _["posteriors"] = posteriors, _["ppm"] = ppm,_["PPI"] = PPI,
                                _["Ks"] = Ks, _["accept_rate_ad"] = accept_r_ad, _["accept_rate_swap"] = accept_r_swap,
                                  _["accept_rate_shift"] = accept_r_shift, _["larger_posterior"] = larger_lklh));
  }else{
    return(Rcpp::List::create(Named("iter") = iter, _["burn"] = burn, _["gamma_map"] = gamma_map,
                              _["gamma_map_index"]= gamma_map_index,  _["Llist"] = Llist,
                              _["hastings"] = hastings, _["posteriors"] = posteriors, 
                              _["Ks"] = Ks, _["accept_rate_ad"] = accept_r_ad, _["accept_rate_swap"] = accept_r_swap,
                                _["accept_rate_shift"] = accept_r_shift, _["larger_posterior"] = larger_lklh));
  }
  
}


// [[Rcpp::export]]
double compute_posterior(NumericMatrix dat,  IntegerVector gamma_0, int estK = 4,
                         double alpha_sigma = 3, double beta_sigma = 0.01, String kernel ="normal", bool open = false,double phi = 1, bool fixedSigma = false){
  int n = dat.nrow();
  if(open){
    gamma_0[0] = 1;
    gamma_0[n-1] = 0;
  }
  //hyper-parameter 
  double alpha_w = 2.0*estK/(n-1.0);
  double beta_w = 2.0 - alpha_w;
  
  //get used variables
  IntegerVector L_0 = gamma2L(gamma_0);
  int K_0 = sum(gamma_0);
  
  //hyper-parameter terms
  double log_K_term = alpha_sigma*log(beta_sigma) - lgamma(alpha_sigma);
  double Log_Gamma_term = LogOmegaterm(n, K_0, alpha_w, beta_w);
  //get segment
  IntegerMatrix df_0 = gamma2seg(gamma_0, open);
  NumericVector df_0_terms(K_0);
  
  for(int i = 0;i<K_0;i++){
    df_0_terms[i] = computeSegLogTerm(df_0(i,0), df_0(i,1), dat, alpha_sigma, beta_sigma, kernel, phi, fixedSigma, open);
  }
  double post_0 = K_0*log_K_term + sum(df_0_terms) +Log_Gamma_term;
  return(post_0);
}


#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]

NumericVector EvaluateSimilarity(NumericMatrix Mat,  NumericVector cluster) {
  
  int n = Mat.nrow();
  int q = Mat.ncol();
  
  NumericVector res(2, 0.0);
  
  for(int i = 0; i < n - 1; ++i){
    for(int j = i + 1; j < n; ++j){
      double dis = 0.0;
      for(int k = 0; k < q; ++k){
        dis += pow(Mat(i, k) - Mat(j, k), 2);
      }
      if(cluster(i) == cluster(j)){
        res(0) += dis;
      } else {
        res(1) += dis;
      }
    }
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat sampleFromMND(int n, int p, arma::mat M, arma::mat U, arma::mat V){
  arma::mat X = arma::randn<arma::mat>(n, p);
  return M + arma::chol(U)*X*arma::chol(V);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat riwishart(int nu, arma::mat V){
  
  int m = V.n_rows;
  
  arma::mat T = arma::zeros<arma::mat>(m, m);
  
  for(int i = 0; i < m; i++) {
    T(i, i) = sqrt(R::rchisq(nu-i));
  }

  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i, j) = R::norm_rand();
    }
  }
  
  arma::mat C = arma::trans(trimatl(T)) * arma::chol(V);
  arma::mat CI = inv(C);
  
  return CI * arma::trans(CI);

}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec dmvnrmArma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
        out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
      
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double dmvnrmRowArma(arma::rowvec x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 

    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    int xdim = x.n_cols;
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    arma::vec z = rooti * arma::trans( x - mean) ;    
    double out = constants - 0.5 * arma::sum(z%z) + rootisum;     
  
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec mvrnormArma(arma::rowvec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::rowvec Y = arma::randn<arma::rowvec>(ncols);
   return mu + Y * arma::chol(sigma);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double callpois(int x, double alpha, bool log){
  Environment stats("package:stats");
  Function dpois = stats["dpois"];
  return as<NumericVector>(wrap(dpois(x, alpha, log)))(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double harmonic(int q){
 double seqx = 0;
  for (int i = 0; i < q; ++i) {
    seqx = seqx + 1/(double(i) +1);
  }
  return seqx; 
}  


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec colmean(arma::mat x) {
  int n = x.n_rows;
  int p = x.n_cols;
  arma::rowvec out = arma::zeros<arma::rowvec>(p);
  for(int i = 0; i < p; ++i){
    out(i) = sum(x.col(i))/n;
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec callKmeans(arma::mat x, int k){
  Environment myEnv("package:stats");
  Function mykmeans = myEnv["kmeans"];
  Rcpp::List kmeansRes = wrap(mykmeans(x, k));
  return as<NumericVector>(kmeansRes["cluster"]);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double Cquantile(arma::vec x, double q) {
  arma::vec y = x;
  sort(y);
  return y(floor(x.n_elem*q));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat fastInverse(arma::mat psi, arma::mat lambda, bool diag) {
  arma::mat psi_inv = psi;
  if(diag == TRUE){
    psi_inv.diag() = 1/(psi.diag()+0.000001);
  } else {
    psi_inv = inv(psi);
  }
  int k = lambda.n_rows;
  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.eye();
  
  arma::mat term_1 = psi_inv*lambda.t();
  arma::mat term_2 = inv(I + lambda*term_1);
  arma::mat term_3 = lambda*psi_inv;
  arma::mat res = psi_inv - term_1*term_2*term_3;
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat initialize_lambda(arma::mat Z_norm, arma::mat X) {
  arma::mat term_1 = inv(Z_norm.t()*Z_norm);
  arma::mat term_2 = Z_norm.t()*X;
  arma::mat lambda = term_1*term_2;
  return lambda;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::uvec reorderMatF(arma::mat F){
  int k = F.n_rows;
  arma::uvec neworder = arma::zeros<arma::uvec>(k);
  if(k > 1){
    arma::vec amount(k);
    for(int i = 0; i < k; ++i){
      amount(i) = sum(F.row(i));      
    }
   neworder = sort_index(amount, "descend");
   arma::vec order_amount = amount.elem(neworder);
   return neworder(arma::find(order_amount > 0));
  } else {
   return neworder;
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double logratio(int active, int q) {
  double ratioP = 0;
  if(active == 0) ratioP = std::log(0.01/(q-1)); 
  if(active == q-1) ratioP = 2;
  if(active < q - 1 && active > 0) ratioP = std::log(active) - std::log(q-1-active);
  return ratioP;
}          

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat subsetMat(arma::mat X, arma::uvec neworder, bool byrow){
  int k = X.n_rows;
  int q = X.n_cols;
  if(byrow == TRUE){
    arma::uvec qq = arma::zeros<arma::uvec>(q);
    return X.submat(neworder, arma::find(qq == 0));
  } else {
    arma::uvec qq = arma::zeros<arma::uvec>(k);
    return X.submat(arma::find(qq == 0), neworder);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List initilizationIBP(arma::mat X, int k) {
  double iniSd = stddev(vectorise(X))/2;
  int q = X.n_cols;
  arma::vec seqx = arma::zeros<arma::vec>(2);
  seqx(1) = 1;
  arma::mat Fini = arma::reshape(RcppArmadillo::sample(seqx, q*k, TRUE), k, q);
  arma::mat Vini = arma::randn<arma::mat>(k, q)*iniSd;
  arma::mat lambda = Fini%Vini; 
  arma::mat I(k, k);
  I.eye();
  arma::mat term_1 = X*lambda.t();
  arma::mat term_2 = inv(lambda*lambda.t() + I);
  arma::mat Qini = term_1*term_2;
  return Rcpp::List::create(Rcpp::Named("inilambda") = lambda, 
    Rcpp::Named("iniF") = Fini, 
    Rcpp::Named("iniV") = Vini, 
    Rcpp::Named("iniQ") = Qini);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSinitilizationIBP(arma::mat Y, arma::mat X, int k1, int k2) {
  // initialization 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Y);
  
  arma::mat matd = arma::zeros<arma::mat>(k1, k1);
  matd.diag() = sqrt(s(arma::span(0, k1-1)));
  arma::mat Z_Y = U.cols(0, k1-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambda_x;
  
  double iniSd = stddev(vectorise(lambda_x));
  int q = X.n_cols;
  
  arma::vec seqx = arma::zeros<arma::vec>(2);
  seqx(1) = 1;
  arma::mat Fini = arma::reshape(RcppArmadillo::sample(seqx, q*k2, TRUE), k2, q);
  arma::mat Vini = arma::randn<arma::mat>(k2, q)*iniSd;
  arma::mat lambda_u = Fini%Vini; 
  arma::mat I(k2, k2);
  I.eye();
  arma::mat term_1 = res_x*lambda_u.t();
  arma::mat term_2 = inv(lambda_u*lambda_u.t() + I);
  arma::mat Qini = term_1*term_2;
  return Rcpp::List::create(Rcpp::Named("lambdaY") = lambda_y, 
    Rcpp::Named("lambdaX") = lambda_x, 
    Rcpp::Named("lambdaU") = lambda_u, 
    Rcpp::Named("Z") = Z_Y_norm, 
    Rcpp::Named("iniF") = Fini, 
    Rcpp::Named("iniV") = Vini, 
    Rcpp::Named("iniQ") = Qini);
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List DirichletSpikeModel(arma::mat X, int K, int iniL, int TruncateL, int iter, 
                               int nu0, double sigma, double r, double s, double alpha, arma::rowvec mu0,
                               arma::mat Sigma0, int kappa0, double m, double g, double h, double c, double d, 
                               double s1, double s2, int iter_to_average){
  
  int n = X.n_rows; 
  int q = X.n_cols;
  arma::mat U;
  arma::vec S;
  arma::mat P;
  
  bool SVDsuccess = false;
  while(SVDsuccess == false) {
    SVDsuccess = svd(U, S, P, X);
    if(SVDsuccess == false){
      X += 1e-4;
    }
  }
  
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(S(arma::span(0, K-1)));
  arma::mat Z_X = U.cols(0, K-1)*matd;
  arma::mat Q = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
  arma::vec C = callKmeans(Q, iniL) - 1; 
  arma::mat lamU = initialize_lambda(Q, X); 
  
  arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
  arma::mat F = arma::conv_to<arma::mat>::from(a);
  arma::mat lambdaU = F%lamU;
  
  arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
  arma::mat Mu = arma::zeros<arma::mat>(TruncateL, K);  
  arma::cube Sigma = arma::zeros<arma::cube>(K, K, TruncateL); 
  arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
  
  for(int l = 0; l < iniL; ++l){
    arma::uvec IDs = find(C == l);
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(Q.rows(IDs));
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(Q.rows(IDs));
    seqx(l) = l;
  }
  
  arma::mat diagK = arma::zeros<arma::mat>(K, K);
  diagK.eye();
  
  arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
  
  for(int l = iniL; l < TruncateL; ++l){
    // Rcpp::Rcout << "test2" << std::endl;
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
    seqx(l) = l;
  }
  
  arma::vec probF(2);
  probF.zeros();
  arma::vec seqF = arma::zeros<arma::vec>(2);
  seqF(1) = 1;
  
  double sigma2inv = 1/(sigma*sigma);
  
  arma::mat psi_x = arma::zeros<arma::mat>(q, q);
  psi_x.eye();
  arma::mat psi_x_inv = psi_x;
  arma::vec psi_x_vector = psi_x.diag();
  arma::vec prob(TruncateL);
  prob.zeros();
  
  arma::mat E_x = X - Q*lambdaU;
  arma::vec v = arma::zeros<arma::vec>(TruncateL);
  arma::vec rho = arma::zeros<arma::vec>(K);
  
  for(int i = 0; i < K; ++i){
    rho(i) = R::rbeta(r*s, s*(1-r));
  }
  
  arma::vec vProd = v;
  arma::vec vProd_minus = v;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    // Sample V_l
    
    for(int l = 0; l < TruncateL; ++l) {
      arma::uvec ID1s = find(C == l);
      arma::uvec ID2s = find(C > l);
      if(any(C == l)){
        v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
      } else {
        v(l) = R::rbeta(1, alpha);
      }	
      if(l == 0){
        vProd_minus(l) = std::log(1-v(l)); 
        vProd(l) = std::log(v(l)); 
      } else {
        vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
        vProd(l) = vProd_minus(l-1) + std::log(v(l));
      }
    }
    
    
    // Sample C_i
    
    for(int i = 0; i < n; ++i) {
      arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
      for(int l = 0; l < TruncateL; ++l) {
        likelihood(l) = vProd(l) + dmvnrmRowArma(Q.row(i), Mu.row(l), Sigma.slice(l), TRUE);
      }
      likelihood = likelihood + abs(max(likelihood));
      likelihood = exp(likelihood);
      prob = likelihood/sum(likelihood);
      C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
    }
    
    // Sample Q_i 
    psi_x_inv.diag() = 1/(psi_x_vector+0.000001);
    arma::mat term_1 = psi_x_inv*lambdaU.t();
    
    for(int i = 0; i < n; ++i) {
      arma::mat inv_Sigma_i = inv(Sigma.slice(C(i)));
      arma::mat Dinv = inv(lambdaU*term_1 + inv_Sigma_i);
      arma::mat part_1 = X.row(i)*term_1 + Mu.row(C(i))*inv_Sigma_i;
      arma::rowvec Wi = part_1*Dinv;
      Q.row(i) = mvrnormArma(Wi, Dinv);
    }
    
    // Sample Mu_l, Sigma_l 
    
    int nl = 0;
    arma::vec currentL(TruncateL);
    currentL.zeros();
    for(int l = 0; l < TruncateL; ++l) {
      arma::mat sumTerm = arma::zeros<arma::mat>(K, K); 
      arma::rowvec muBar = arma::zeros<arma::rowvec>(K);
      if(any(C == l)){
        currentL(l) = 1;
        arma::uvec IDs = find(C == l);
        muBar = colmean(Q.rows(IDs));
        nl = IDs.n_elem;
        for(int ii = 0; ii < nl; ii++){
          arma::rowvec ll = Q.row(IDs(ii)) - muBar;
          sumTerm = sumTerm + ll.t()*ll;
        }
      } 		
      Sigma.slice(l) = riwishart(nu0 + nl, Sigma0 + sumTerm + kappa0*nl*(muBar - mu0).t()*(muBar - mu0)/(kappa0 + nl));	
      Mu.row(l) = mvrnormArma((kappa0*mu0 + nl*muBar)/(kappa0 + nl), Sigma.slice(l)/(kappa0 + nl));
    }	
    
    // Rcpp::Rcout << "2" << std::endl;
    
    // Sample F    
    for(int j = 0; j < q; ++j){    
      
      arma::mat B = arma::zeros<arma::mat>(n, K);
      
      for(int k = 0; k < K; ++k){           
        
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < K; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*lambdaU(l, j);
            } 
          }  
          B(i, k) = X(i, j) - QFVsum;
        }
        
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        double ratiop = std::log(m*rho(k)/(1 - m*rho(k)));
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + ratiop;
        probF(0) = 1/(std::exp(ratio)+1);
        probF(1) = 1 - probF(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
        if(F(k, j) == 1){
          lambdaU(k, j) = R::rnorm(mukj, sqrt(sigmakj));
        } else {
          lambdaU(k, j) = 0;
        }
      }  
    }
    // Rcpp::Rcout << sum(vectorise(F)) << std::endl;
    
    // Sample rho_k 
    for(int k = 0; k < K; ++k){ 
      rho(k) = R::rbeta(r*s + sum(F.row(k)), s*(1-r) + q - sum(F.row(k)));        
    }    
    
    // Sample psi_x 
    E_x = X - Q*lambdaU;
    
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
    }
    psi_x.diag() = psi_x_vector;
    h = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(F))/2, 1/(d + sum(pow(lambdaU(arma::find(F == 1)), 2))/2));;	
    sigma = sqrt(1/sigma2inv);
    
    d = R::rgamma(1 + K, 1/(1 + sigma2inv*K));		
    
    double vProd_minus_add = 0;
    double counter = 0;
    for(int l = 0; l < TruncateL; ++l) {
      if(currentL(l) == 1){
        vProd_minus_add = vProd_minus_add + std::log(1-v(l));
        counter = counter + 1;
      }
    }
    
    // Sample alpha
    alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
    
    
    if(iter - kk <= iter_to_average){
      for(int i = 0; i < n; ++i) {
        AveC(C(i), i) = AveC(C(i), i) + 1;
      }
    }
  }
  
  AveC = AveC/iter_to_average;   
  
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q,
    Rcpp::Named("F") = F,
    Rcpp::Named("C") = C,
    Rcpp::Named("AveC") = AveC,
    Rcpp::Named("lambdaU") = lambdaU, 
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("psi") = psi_x_vector); 
}

// A clustering model with sparse factors with spike and slab prior
// model dropout

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List DirichletSpikeModelZero(arma::mat X, arma::mat H, arma::vec tau, double kappa_int, arma::vec kappagrid, 
                    int K, int iniL, int TruncateL, int iter, int nu0, double sigma, double r, double s, 
                    double alpha, arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, double h, 
                    double c, double d, double s1, double s2, int iter_to_average){
  
  int n = X.n_rows; 
  int q = X.n_cols;
  arma::mat U;
  arma::vec S;
  arma::mat P;
  double kappa = kappa_int;
  int kappan = kappagrid.n_elem;
  
  bool SVDsuccess = false;
  while(SVDsuccess == false) {
    SVDsuccess = svd(U, S, P, X);
    if(SVDsuccess == false){
      X += 1e-4;
    }
  }
  
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(S(arma::span(0, K-1)));
  arma::mat Z_X = U.cols(0, K-1)*matd;
  arma::mat Q = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
  arma::vec C = callKmeans(Q, iniL) - 1; 
  arma::mat lamU = initialize_lambda(Q, X); 
  
  arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
  arma::mat F = arma::conv_to<arma::mat>::from(a);
  arma::mat lambdaU = F%lamU;
  
  arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
  arma::mat Mu = arma::zeros<arma::mat>(TruncateL, K);  
  arma::cube Sigma = arma::zeros<arma::cube>(K, K, TruncateL); 
  arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
  
  for(int l = 0; l < iniL; ++l){
    arma::uvec IDs = find(C == l);
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(Q.rows(IDs));
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(Q.rows(IDs));
    seqx(l) = l;
  }
  
  arma::mat diagK = arma::zeros<arma::mat>(K, K);
  diagK.eye();
  
  arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
  
  for(int l = iniL; l < TruncateL; ++l){
    // Rcpp::Rcout << "test2" << std::endl;
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
    seqx(l) = l;
  }
  
  arma::vec probF(2);
  probF.zeros();
  arma::vec seqF = arma::zeros<arma::vec>(2);
  seqF(1) = 1;
  
  double sigma2inv = 1/(sigma*sigma);
  
  arma::mat psi_x = arma::zeros<arma::mat>(q, q);
  psi_x.eye();
  arma::mat psi_x_inv = psi_x;
  arma::vec psi_x_vector = psi_x.diag();
  arma::vec prob(TruncateL);
  prob.zeros();
  
  arma::mat E_x = X - Q*lambdaU;
  arma::vec v = arma::zeros<arma::vec>(TruncateL);
  arma::vec rho = arma::zeros<arma::vec>(K);
  
  for(int i = 0; i < K; ++i){
    rho(i) = R::rbeta(r*s, s*(1-r));
  }
  
  arma::vec vProd = v;
  arma::vec vProd_minus = v;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    // Sample V_l
    
    for(int l = 0; l < TruncateL; ++l) {
      arma::uvec ID1s = find(C == l);
      arma::uvec ID2s = find(C > l);
      if(any(C == l)){
        v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
      } else {
        v(l) = R::rbeta(1, alpha);
      }	
      if(l == 0){
        vProd_minus(l) = std::log(1-v(l)); 
        vProd(l) = std::log(v(l)); 
      } else {
        vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
        vProd(l) = vProd_minus(l-1) + std::log(v(l));
      }
    }
    
    
    // Sample C_i
    
    for(int i = 0; i < n; ++i) {
      arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
      for(int l = 0; l < TruncateL; ++l) {
        likelihood(l) = vProd(l) + dmvnrmRowArma(Q.row(i), Mu.row(l), Sigma.slice(l), TRUE);
      }
      likelihood = likelihood + abs(max(likelihood));
      likelihood = exp(likelihood);
      prob = likelihood/sum(likelihood);
      C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
    }
    
    // Sample Q_i 
    psi_x_inv.diag() = 1/(psi_x_vector+0.000001);
    arma::mat term_1 = psi_x_inv*lambdaU.t();
    
    for(int i = 0; i < n; ++i) {
      arma::mat inv_Sigma_i = inv(Sigma.slice(C(i)));
      arma::mat Dinv = inv(lambdaU*term_1 + inv_Sigma_i);
      arma::mat part_1 = X.row(i)*term_1 + Mu.row(C(i))*inv_Sigma_i;
      arma::rowvec Wi = part_1*Dinv;
      Q.row(i) = mvrnormArma(Wi, Dinv);
    }
    
    // Sample Mu_l, Sigma_l 
    
    int nl = 0;
    arma::vec currentL(TruncateL);
    currentL.zeros();
    for(int l = 0; l < TruncateL; ++l) {
      arma::mat sumTerm = arma::zeros<arma::mat>(K, K); 
      arma::rowvec muBar = arma::zeros<arma::rowvec>(K);
      if(any(C == l)){
        currentL(l) = 1;
        arma::uvec IDs = find(C == l);
        muBar = colmean(Q.rows(IDs));
        nl = IDs.n_elem;
        for(int ii = 0; ii < nl; ii++){
          arma::rowvec ll = Q.row(IDs(ii)) - muBar;
          sumTerm = sumTerm + ll.t()*ll;
        }
      } 		
      Sigma.slice(l) = riwishart(nu0 + nl, Sigma0 + sumTerm + kappa0*nl*(muBar - mu0).t()*(muBar - mu0)/(kappa0 + nl));	
      Mu.row(l) = mvrnormArma((kappa0*mu0 + nl*muBar)/(kappa0 + nl), Sigma.slice(l)/(kappa0 + nl));
    }	
    
    // Rcpp::Rcout << "2" << std::endl;
    
    // Sample F    
    for(int j = 0; j < q; ++j){    
      
      arma::mat B = arma::zeros<arma::mat>(n, K);
      
      for(int k = 0; k < K; ++k){           
        
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < K; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*lambdaU(l, j);
            } 
          }  
          B(i, k) = X(i, j) - QFVsum;
        }
        
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        double ratiop = std::log(m*rho(k)/(1 - m*rho(k)));
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + ratiop;
        probF(0) = 1/(std::exp(ratio)+1);
        probF(1) = 1 - probF(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
        if(F(k, j) == 1){
          lambdaU(k, j) = R::rnorm(mukj, sqrt(sigmakj));
        } else {
          lambdaU(k, j) = 0;
        }
      }  
    }
    // Rcpp::Rcout << sum(vectorise(F)) << std::endl;
    
    // Sample rho_k 
    for(int k = 0; k < K; ++k){ 
      rho(k) = R::rbeta(r*s + sum(F.row(k)), s*(1-r) + q - sum(F.row(k)));        
    }    
    
    // Sample psi_x 
    E_x = X - Q*lambdaU;
    
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
    }
    psi_x.diag() = psi_x_vector;
    h = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(F))/2, 1/(d + sum(pow(lambdaU(arma::find(F == 1)), 2))/2));;	
    sigma = sqrt(1/sigma2inv);
    
    d = R::rgamma(1 + K, 1/(1 + sigma2inv*K));		
    
    double vProd_minus_add = 0;
    double counter = 0;
    for(int l = 0; l < TruncateL; ++l) {
      if(currentL(l) == 1){
        vProd_minus_add = vProd_minus_add + std::log(1-v(l));
        counter = counter + 1;
      }
    }
    
    // Sample alpha
    alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
    
    arma::vec kappa_likelihood = arma::zeros<arma::vec>(kappan);
    
    // Sample X_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < q; ++j){
        if(H(i, j) == 1){
          double term_ij = 2*kappa*psi_x_vector(j) + 1;
          double mu_ij = (as_scalar(Q.row(i)*lambdaU.col(j)) - 2*kappa*psi_x_vector(j)*tau(j))/term_ij ;
          double sigma_ij = psi_x_vector(j)/term_ij;
          X(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(X(i, j) + tau(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(X(i, j) + tau(j), 2)));
          }
        }
        
      }
    }
    
   kappa_likelihood = kappa_likelihood + abs(max(kappa_likelihood));
   kappa_likelihood = exp(kappa_likelihood);
   arma::vec kappa_prob = kappa_likelihood/sum(kappa_likelihood);
   kappa = as_scalar(RcppArmadillo::sample(kappagrid, 1, FALSE, as<NumericVector>(wrap(kappa_prob))));
    
    if(iter - kk <= iter_to_average){
      for(int i = 0; i < n; ++i) {
        AveC(C(i), i) = AveC(C(i), i) + 1;
      }
    }
  }
  
  AveC = AveC/iter_to_average;   
  
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q,
    Rcpp::Named("F") = F,
    Rcpp::Named("C") = C,
    Rcpp::Named("AveC") = AveC,
    Rcpp::Named("lambdaU") = lambdaU, 
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("psi") = psi_x_vector,
    Rcpp::Named("kappa") = kappa); 
}


// a latent factor model using an IBP prior 

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List IBPfactormodel(arma::mat X, int numk, int iter, double sigma, double alpha, 
                          double kappa, double g, double h, double c, double d, int limit){
  
  Rcpp::List res = initilizationIBP(X, numk); 
  arma::mat iniF = res["iniF"];
  arma::mat iniV = res["iniV"];
  arma::mat iniQ = res["iniQ"];
  arma::mat inilambda = res["inilambda"];
  
  int activeK = numk; 
  int n = X.n_rows;
  int q = X.n_cols;
  double har = harmonic(q);
  
  arma::mat recordF = arma::zeros<arma::mat>(limit, q);
  arma::mat recordV = arma::zeros<arma::mat>(limit, q);
  arma::mat recordlambda = arma::zeros<arma::mat>(limit, q);
  arma::mat recordQ = arma::zeros<arma::mat>(n, limit);
  
  recordF.rows(0, activeK-1) = iniF;
  recordV.rows(0, activeK-1) = iniV;
  recordlambda.rows(0, activeK-1) = inilambda;   
  recordQ.cols(0, activeK - 1) = iniQ;
  
  arma::mat psi = arma::zeros<arma::mat>(q, q);
  psi.eye();
  double sigma2inv = 1/(sigma*sigma);
  arma::vec prob(2);
  prob.zeros();
  arma::vec seqx = arma::zeros<arma::vec>(2);
  seqx(1) = 1;
  
  arma::mat E_x = X - iniQ*inilambda;
  arma::vec psi_x_vector = psi.diag();
  
  for (int kk = 0; kk < iter; ++kk) { 
    //Sample F;
    
    for(int j = 0; j < q; ++j){  
      
      arma::mat F = recordF.rows(0, activeK - 1);
      arma::mat V = recordV.rows(0, activeK - 1);
      arma::mat lambda = recordlambda.rows(0, activeK - 1);
      arma::mat Q = recordQ.cols(0, activeK - 1);
      arma::mat B = arma::zeros<arma::mat>(n, activeK);
      
      for(int k = 0; k < activeK; ++k){           
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < activeK; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
            } 
          }  
          B(i, k) = X(i, j) - QFVsum;
        }
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        arma::uvec allselect = find(F.row(k) == 1);
        int active = 0;
        if( F(k, j) == 1 ){
          active = allselect.n_elem - 1;
        } else {
          active = allselect.n_elem;
        }
        
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
        prob(0) = 1/(std::exp(ratio)+1);
        prob(1) = 1 - prob(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob)))); 
        
        if(F(k, j) == 1){
          	V(k, j) = R::rnorm(mukj, sqrt(sigmakj));
        }  
        lambda(k, j) = V(k, j)*F(k, j);
        recordF(k, j) = F(k, j);
        recordV(k, j) = V(k, j);
        recordlambda(k, j) = lambda(k, j); 
      }  
    
        // Generate new features
      int knew = R::rpois(alpha/(q-1));
		    
      if(knew != 0){
        if(knew + activeK > limit) break;
        double ap = 0;
        if(kappa == 0){
          ap = - callpois(knew, alpha/(q-1), TRUE);
        } else {
          ap = callpois(knew, alpha/(q-1), TRUE) - callpois(knew, kappa*alpha/(q-1), TRUE); 
        }
        arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
        FAugmented.col(j).fill(1);
        arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
        arma::mat lambdaAugmented = FAugmented%VAugmented;
        arma::mat knewI = arma::zeros<arma::mat>(knew, knew);
        knewI.eye(); 
        arma::mat term_1 = inv(psi)*lambdaAugmented.t();
        arma::mat D = lambdaAugmented*term_1 + knewI;
        double CDCsum = 0; 
        arma::mat resQ = X - Q*lambda; 
        arma::mat tt = term_1*inv(D);
        for(int i = 0; i < n; ++i){
          arma::mat ct = resQ.row(i)*tt;
          CDCsum = CDCsum + as_scalar(ct*D*ct.t());
        }
        double lratio = - n/2*std::log(det(D)) + 0.5*CDCsum + ap;
                
        int kfinal = 0;  
        if(std::exp(lratio) > 1){
          kfinal = knew;
        }   
      	if(kfinal != 0){  
          recordF.rows(activeK, activeK + kfinal - 1) = FAugmented;
          recordV.rows(activeK, activeK + kfinal - 1) = VAugmented;
          recordlambda.rows(activeK, activeK + kfinal - 1) = lambdaAugmented;		
          arma::mat sharedterm = fastInverse(psi, lambdaAugmented, TRUE);
          arma::mat term_2 = lambdaAugmented*sharedterm;
          arma::mat E_Q_w = (term_2*resQ.t()).t();
          arma::mat Var_Q_w = knewI - term_2*lambdaAugmented.t();
		      arma::mat Qalter = arma::zeros<arma::mat>(n, knew);
          for(int i = 0; i < n; ++i){
					   Qalter.row(i) = mvrnormArma(E_Q_w.row(i), Var_Q_w);
				  }
          recordQ.cols(activeK, activeK + kfinal-1) = Qalter;
          activeK = activeK + kfinal;
      	}
      }         
    } 	 
    
            
    arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
    arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
    arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat lambdanew = subsetMat(recordlambda.rows(0, activeK - 1), dseq, TRUE); 
  
    activeK = Fnew.n_rows;

    arma::mat sharedterm = fastInverse(psi, lambdanew, TRUE);
    arma::mat term_3 = lambdanew*sharedterm;
    arma::mat E_Q_w = (term_3*X.t()).t();
    arma::mat activeI = arma::zeros<arma::mat>(activeK, activeK);
    activeI.eye();
    arma::mat Var_Q_w = activeI - term_3*lambdanew.t();
  	arma::mat Qnew = arma::zeros<arma::mat>(n, activeK);
    for(int i = 0; i < n; ++i){
				Qnew.row(i) = mvrnormArma(E_Q_w.row(i), Var_Q_w);
		}
     
    E_x = X - Qnew*lambdanew;
    //Rcpp::Rcout << sum(vectorise(E_x.t()*E_x)) << std::endl;
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
	  }
	  psi.diag() = psi_x_vector;
	  h = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(Vnew(arma::find(Fnew == 1)), 2))/2));;	
	  sigma = sqrt(1/sigma2inv);
	  d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));
    
    //Sampling alpha 
   	alpha = R::rgamma(1 + activeK, 1/(1 + har));
     
    recordF.rows(0, activeK-1) = Fnew;
    recordV.rows(0, activeK-1) = Vnew;
    recordlambda.rows(0, activeK-1) = lambdanew;   
    recordQ.cols(0, activeK-1) = Qnew;
     
    //Rcpp::Rcout << stddev(Vnew(arma::find(Fnew == 1))) << std::endl;
    //Rcpp::Rcout << psi_x_vector << std::endl;
  }  
  
  return Rcpp::List::create(Rcpp::Named("lambda") = recordlambda.rows(0, activeK-1), 
          Rcpp::Named("F") = recordF.rows(0, activeK-1), 
          Rcpp::Named("V") = recordV.rows(0, activeK-1), 
          Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
          Rcpp::Named("activeK") = activeK);
}  

//  A mixture factor model using a dirichlet process prior on the number of clusters (mixture components)
// and an IBP prior on the number of latent factors.
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  
  Rcpp::List DirichletIBPModel(arma::mat X, int K, int iniL, int TruncateL, int iter, int maxK, 
                               int nu0, double sigma, double r, double s, double alpha, double alpha2, 
                               arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, double h, 
                               double c, double d, double kappa, double s1, double s2, int iter_to_average){
    int n = X.n_rows; 
    int q = X.n_cols;
    arma::mat U;
    arma::vec S;
    arma::mat P;
    svd(U, S, P, X);
    
    double har = harmonic(q);
    
    arma::mat matd = arma::zeros<arma::mat>(K, K);
    matd.diag() = sqrt(S(arma::span(0, K-1)));
    arma::mat Z_X = U.cols(0, K-1)*matd;
    arma::mat iniQ = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
    arma::vec C = callKmeans(iniQ, iniL) - 1; 
    arma::mat lamU = initialize_lambda(iniQ, X);
    
    arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
    arma::mat iniF = arma::conv_to<arma::mat>::from(a);
    arma::mat inilambdaU = iniF%lamU;
    arma::mat iniV = inilambdaU;
    
    arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
    arma::mat Mu = arma::zeros<arma::mat>(TruncateL, maxK);  
    arma::cube Sigma = arma::zeros<arma::cube>(maxK, maxK, TruncateL); 
    arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
    
    for(int l = 0; l < iniL; ++l){
      arma::uvec IDs = find(C == l);
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(iniQ.rows(IDs));
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(iniQ.rows(IDs));
      seqx(l) = l;
    }
    
    arma::mat diagK = arma::eye(K, K);
    
    arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
    
    for(int l = iniL; l < TruncateL; ++l){
      
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
      seqx(l) = l;
    }
    
    arma::vec probF(2);
    probF.zeros();
    arma::vec seqF = arma::zeros<arma::vec>(2);
    seqF(1) = 1;
    
    double sigma2inv = 1/(sigma*sigma);
    
    arma::mat psi_x = arma::eye(q, q);
    arma::mat psi_x_inv = psi_x;
    arma::vec psi_x_vector = psi_x.diag();
    arma::vec prob(TruncateL);
    prob.zeros();
    
    arma::mat E_x = X - iniQ*inilambdaU;
    arma::vec v = arma::zeros<arma::vec>(TruncateL);
    
    arma::vec vProd = v;
    arma::vec vProd_minus = v;
    
    arma::mat recordF = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordV = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordlambdaU = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordQ = arma::zeros<arma::mat>(n, maxK);
    
    int activeK = K;
    recordF.rows(0, activeK-1) = iniF;
    recordV.rows(0, activeK-1) = iniV;
    recordlambdaU.rows(0, activeK-1) = inilambdaU;   
    recordQ.cols(0, activeK - 1) = iniQ;
    
    for (int kk = 0; kk < iter; ++kk) { 
      
      // Sample F    
      for(int j = 0; j < q; ++j){    
        
        arma::mat F = recordF.rows(0, activeK - 1);
        arma::mat V = recordV.rows(0, activeK - 1);
        arma::mat lambdaU = recordlambdaU.rows(0, activeK - 1);
        arma::mat Q = recordQ.cols(0, activeK - 1);
        arma::mat B = arma::zeros<arma::mat>(n, activeK);
        
        for(int k = 0; k < activeK; ++k){           
          
          for(int i = 0; i < n; ++i){          
            double QFVsum = 0;
            for(int l = 0; l < activeK; ++l){
              if(l != k){
                QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
              } 
            }  
            B(i, k) = X(i, j) - QFVsum;
          }
          
          double sumQB = 0;
          for(int i = 0; i < n; ++i){
            sumQB = sumQB + Q(i, k)*B(i, k);
          }
          
          double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
          double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
          double mukj = sumQB*sigmakj/psi_x_vector(j);
          arma::uvec allselect = find(F.row(k) == 1);
          int active = 0;
          if( F(k, j) == 1 ){
            active = allselect.n_elem - 1;
          } else {
            active = allselect.n_elem;
          }
          
          double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
          
          probF(0) = 1/(std::exp(ratio)+1);
          probF(1) = 1 - probF(0); 
          F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
          
          if(F(k, j) == 1){
            V(k, j)  = R::rnorm(mukj, sqrt(sigmakj));
          } 
          
          lambdaU(k, j) = V(k, j)*F(k, j);
          recordF(k, j) = F(k, j);
          recordV(k, j) = V(k, j);
          recordlambdaU(k, j) = lambdaU(k, j); 
        }  
        // Generate new features
        int knew = R::rpois(alpha2/(q-1));
        
        if(knew != 0){
          if(knew + activeK > maxK) break;
          double ap = 0;
          if(kappa == 0){
            ap = - callpois(knew, alpha2/(q-1), TRUE);
          } else {
            ap = callpois(knew, alpha2/(q-1), TRUE) - callpois(knew, kappa*alpha2/(q-1), TRUE); 
          }
          arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
          FAugmented.col(j).fill(1);
          arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
          arma::mat lambdaUAugmented = FAugmented%VAugmented;
          
          psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
          
          arma::mat term_4 = psi_x_inv*lambdaUAugmented.t();
          
          arma::cube H_i = arma::zeros<arma::cube>(knew, knew, TruncateL); 
          arma::cube G_i_firstterm = arma::zeros<arma::cube>(q, knew, TruncateL); 
          arma::mat G_i_secondterm = arma::zeros<arma::mat>(TruncateL, knew); 
          
          arma::mat H_i_firstterm = lambdaUAugmented*term_4;
          arma::vec mu_sigma_mu = arma::zeros<arma::vec>(TruncateL);
          arma::cube Sigma_augmented_list = arma::zeros<arma::cube>(knew, knew, TruncateL);
          arma::mat Mu_augmented_list = arma::zeros<arma::mat>(TruncateL, knew);
          
          
          for(int l = 0; l < TruncateL; ++l){
            arma::mat knewI = arma::eye(knew, knew);
            Sigma_augmented_list.slice(l) = riwishart(knew + nu0, knewI); 
            Mu_augmented_list.row(l) = mvrnormArma(arma::zeros<arma::rowvec>(knew), Sigma_augmented_list.slice(l)/kappa0);
            if(any(C==l)){
              arma::rowvec Mu_augmented = Mu_augmented_list.row(l);
              arma::mat Sigma_augmented = Sigma_augmented_list.slice(l);
              arma::mat inv_Sigma_augmented = inv(Sigma_augmented);
              arma::mat term_2 = Mu_augmented*inv_Sigma_augmented;
              mu_sigma_mu(l) = as_scalar(term_2*Mu_augmented.t());
              H_i.slice(l) = H_i_firstterm + inv_Sigma_augmented;
              arma::mat term_3 = inv(H_i.slice(l));
              G_i_secondterm.row(l) = term_2*term_3;
              G_i_firstterm.slice(l) = term_4*term_3; 		     
            }	
          }
          
          E_x = X - Q*lambdaU;
          double GHG_sum = 0; 
          double mu_sigma_mu_sum = 0;
          double det_H = 0;
          
          for(int i = 0; i < n; ++i){ 
            arma::rowvec gt = E_x.row(i)*G_i_firstterm.slice(C(i)) + G_i_secondterm.row(C(i));
            GHG_sum = GHG_sum + as_scalar(gt*H_i.slice(C(i))*gt.t());
            mu_sigma_mu_sum = mu_sigma_mu_sum + mu_sigma_mu(C(i));
            det_H = det_H + std::log(det(H_i.slice(C(i)))) + log(det(Sigma_augmented_list.slice(C(i))));
          }					
          double lratio  = - 0.5*det_H + 0.5*(GHG_sum - mu_sigma_mu_sum) + ap;
          
          int kfinal = 0;  
          if(std::exp(lratio) > 1){
            kfinal = knew;
          }   
          if(kfinal != 0){  
            int newK = activeK + kfinal;
            recordF.rows(activeK, newK - 1) = FAugmented;
            recordV.rows(activeK, newK - 1) = VAugmented;
            recordlambdaU.rows(activeK, newK - 1) = lambdaUAugmented;		
            
            arma::mat Q_alter = arma::zeros<arma::mat>(n, knew);
            arma::mat term_5 = psi_x_inv*lambdaUAugmented.t();
            arma::mat D_inv_firstterm = lambdaUAugmented*term_5;
            
            for(int i = 0; i < n; ++i){ 
              arma::mat term_6 = inv(Sigma_augmented_list.slice(C(i)));
              arma::mat D_inv = inv(D_inv_firstterm + term_6);
              arma::mat sum_0 = E_x.row(i)*term_5 + Mu_augmented_list.row(C(i))*term_6;
              arma::rowvec W_i = sum_0*D_inv;
              Q_alter.row(i) = mvrnormArma(W_i, D_inv);
            }
            
            recordQ.cols(activeK, newK - 1) = Q_alter;
            Mu.cols(activeK, newK - 1) = Mu_augmented_list;               
            activeK = newK;
          }
        }         
      }
      //Rcpp::Rcout << "testf" << std::endl;
      arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
      arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
      arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
      arma::mat lambdaUnew = subsetMat(recordlambdaU.rows(0, activeK - 1), dseq, TRUE); 
      arma::mat Qnew = subsetMat(recordQ.cols(0, activeK - 1), dseq, FALSE);
      
      activeK = Fnew.n_rows;
      
      arma::mat Munew = arma::zeros<arma::mat>(TruncateL, activeK);
      arma::cube Sigmanew = arma::zeros<arma::cube>(activeK, activeK, TruncateL); 
      
      //Rcpp::Rcout << "test!" << std::endl; 
      // Sample Mu_l, Sigma_l  
      int nl = 0;
      arma::vec currentL(TruncateL);
      currentL.zeros();
      for(int l = 0; l < TruncateL; ++l) {
        arma::mat sumTerm = arma::zeros<arma::mat>(activeK, activeK); 
        arma::rowvec muBar = arma::zeros<arma::rowvec>(activeK);
        if(any(C == l)){
          currentL(l) = 1;
          arma::uvec IDs = find(C == l);
          muBar = colmean(Qnew.rows(IDs));
          nl = IDs.n_elem;
          for(int ii = 0; ii < nl; ii++){
            arma::rowvec ll = Qnew.row(IDs(ii)) - muBar;
            sumTerm = sumTerm + ll.t()*ll;
          }  
        }  
        Sigmanew.slice(l) = riwishart(nu0 + nl, arma::eye(activeK, activeK) + sumTerm + kappa0*nl*muBar.t()*muBar/(kappa0 + nl));  
        Munew.row(l) = mvrnormArma(nl*muBar/(kappa0 + nl), Sigmanew.slice(l)/(kappa0 + nl)); //mu0 = 0
      }	
      
      // Sample V_l
      
      for(int l = 0; l < TruncateL; ++l) {
        arma::uvec ID1s = find(C == l);
        arma::uvec ID2s = find(C > l);
        if(any(C == l)){
          v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
        } else {
          v(l) = R::rbeta(1, alpha);
        }	
        if(l == 0){
          vProd_minus(l) = std::log(1-v(l)); 
          vProd(l) = std::log(v(l)); 
        } else {
          vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
          vProd(l) = vProd_minus(l-1) + std::log(v(l));
        }
      }
      
      // Sample Q_i 
      psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
      
      //Rcpp::Rcout << activeK << std::endl;
      arma::mat term_7 = psi_x_inv*lambdaUnew.t();
      for(int i = 0; i < n; ++i) {
        arma::mat inv_Sigmanew_i = inv(Sigmanew.slice(C(i)));
        arma::mat Dinv = inv(lambdaUnew*term_7 + inv_Sigmanew_i);
        arma::mat term_8 = X.row(i)*term_7 + Munew.row(C(i))*inv_Sigmanew_i;
        arma::rowvec Wi = term_8*Dinv;
        Qnew.row(i) = mvrnormArma(Wi, Dinv);
      }
      
      // Sample C_i
      for(int i = 0; i < n; ++i) {
        arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
        for(int l = 0; l < TruncateL; ++l) {
          likelihood(l) = vProd(l) + dmvnrmRowArma(Qnew.row(i), Munew.row(l), Sigmanew.slice(l), TRUE);
        }
        likelihood = likelihood + abs(max(likelihood));
        likelihood = exp(likelihood);
        prob = likelihood/sum(likelihood);
        C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
        
      }
      
      E_x = X - Qnew*lambdaUnew;
      
      for(int j = 0; j < q; ++j){
        psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
        if(psi_x_vector(j) < 0.000001){
          psi_x_vector(j) = 0.000001;
        }
      }
      
      psi_x.diag() = psi_x_vector;
      h = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
      
      //Sample sigma
      sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(lambdaUnew(arma::find(Fnew == 1)), 2))/2));	
      sigma = sqrt(1/sigma2inv);
      
      d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));		
      
      double vProd_minus_add = 0;
      double counter = 0;
      for(int l = 0; l < TruncateL; ++l) {
        if(currentL(l) == 1){
          vProd_minus_add = vProd_minus_add + std::log(1-v(l));
          counter = counter + 1;
        }
      }
      
      // Sample alpha
      alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
      
      // Sample alpha2 
      alpha2 = R::rgamma(1 + activeK, 1/(1 + har));
      
      recordF.rows(0, activeK-1) = Fnew;
      recordV.rows(0, activeK-1) = Vnew;
      recordlambdaU.rows(0, activeK-1) = lambdaUnew;   
      recordQ.cols(0, activeK-1) = Qnew; 
      Mu.cols(0, activeK-1) = Munew; 
      
      for(int l = 0; l < TruncateL; ++l){
        Sigma.slice(l).submat(arma::span(0, activeK-1), arma::span(0, activeK-1)) = Sigmanew.slice(l);
      }
      
      if(iter - kk <= iter_to_average){
        for(int i = 0; i < n; ++i) {
          AveC(C(i), i) = AveC(C(i), i) + 1;
        }
      }
    }    
    
    AveC = AveC/iter_to_average;   
    
    return Rcpp::List::create(
      Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
      Rcpp::Named("F") = recordF.rows(0, activeK-1),
      Rcpp::Named("C") = C,
      Rcpp::Named("AveC") = AveC,
      Rcpp::Named("lambdaU") = recordlambdaU.rows(0, activeK-1),
      Rcpp::Named("Mu") = Mu,
      Rcpp::Named("Sigma") = Sigma,
      Rcpp::Named("sigma") = sigma,
      Rcpp::Named("psi") = psi_x_vector); 
  }


//  A mixture factor model using a dirichlet process prior on the number of clusters (mixture components)
// and an IBP prior on the number of latent factors, also add the modeling of dropout.

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List DirichletIBPModelZero(arma::mat X, arma::mat H, arma::vec tau, double kappa_int, arma::vec kappagrid, int K, int iniL, int TruncateL, int iter, int maxK, 
                             int nu0, double sigma, double r, double s, double alpha, double alpha2, 
                             arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, double h, 
                             double c, double d, double kappa_ibp, double s1, double s2, int iter_to_average){
  int n = X.n_rows; 
  int q = X.n_cols;
  arma::mat U;
  arma::vec S;
  arma::mat P;
  svd(U, S, P, X);
  
  double kappa = kappa_int;
  int kappan = kappagrid.n_elem;
  
  double har = harmonic(q);
  
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(S(arma::span(0, K-1)));
  arma::mat Z_X = U.cols(0, K-1)*matd;
  arma::mat iniQ = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
  arma::vec C = callKmeans(iniQ, iniL) - 1; 
  arma::mat lamU = initialize_lambda(iniQ, X);
  
  arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
  arma::mat iniF = arma::conv_to<arma::mat>::from(a);
  arma::mat inilambdaU = iniF%lamU;
  arma::mat iniV = inilambdaU;
  
  arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
  arma::mat Mu = arma::zeros<arma::mat>(TruncateL, maxK);  
  arma::cube Sigma = arma::zeros<arma::cube>(maxK, maxK, TruncateL); 
  arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
  
  for(int l = 0; l < iniL; ++l){
    arma::uvec IDs = find(C == l);
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(iniQ.rows(IDs));
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(iniQ.rows(IDs));
    seqx(l) = l;
  }
  
  arma::mat diagK = arma::eye(K, K);
  
  arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
  
  for(int l = iniL; l < TruncateL; ++l){
    
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
    seqx(l) = l;
  }
  
  arma::vec probF(2);
  probF.zeros();
  arma::vec seqF = arma::zeros<arma::vec>(2);
  seqF(1) = 1;
  
  double sigma2inv = 1/(sigma*sigma);
  
  arma::mat psi_x = arma::eye(q, q);
  arma::mat psi_x_inv = psi_x;
  arma::vec psi_x_vector = psi_x.diag();
  arma::vec prob(TruncateL);
  prob.zeros();
  
  arma::mat E_x = X - iniQ*inilambdaU;
  arma::vec v = arma::zeros<arma::vec>(TruncateL);
  
  arma::vec vProd = v;
  arma::vec vProd_minus = v;
  
  arma::mat recordF = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordV = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordlambdaU = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordQ = arma::zeros<arma::mat>(n, maxK);
  
  int activeK = K;
  recordF.rows(0, activeK-1) = iniF;
  recordV.rows(0, activeK-1) = iniV;
  recordlambdaU.rows(0, activeK-1) = inilambdaU;   
  recordQ.cols(0, activeK - 1) = iniQ;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    // Sample F    
    for(int j = 0; j < q; ++j){    
      
      arma::mat F = recordF.rows(0, activeK - 1);
      arma::mat V = recordV.rows(0, activeK - 1);
      arma::mat lambdaU = recordlambdaU.rows(0, activeK - 1);
      arma::mat Q = recordQ.cols(0, activeK - 1);
      arma::mat B = arma::zeros<arma::mat>(n, activeK);
      
      for(int k = 0; k < activeK; ++k){           
        
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < activeK; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
            } 
          }  
          B(i, k) = X(i, j) - QFVsum;
        }
        
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        arma::uvec allselect = find(F.row(k) == 1);
        int active = 0;
        if( F(k, j) == 1 ){
          active = allselect.n_elem - 1;
        } else {
          active = allselect.n_elem;
        }
        
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
        
        probF(0) = 1/(std::exp(ratio)+1);
        probF(1) = 1 - probF(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
        
        if(F(k, j) == 1){
          V(k, j)  = R::rnorm(mukj, sqrt(sigmakj));
        } 
        
        lambdaU(k, j) = V(k, j)*F(k, j);
        recordF(k, j) = F(k, j);
        recordV(k, j) = V(k, j);
        recordlambdaU(k, j) = lambdaU(k, j); 
      }  
      // Generate new features
      int knew = R::rpois(alpha2/(q-1));
      
      if(knew != 0){
        if(knew + activeK > maxK) break;
        double ap = 0;
        if(kappa_ibp == 0){
          ap = - callpois(knew, alpha2/(q-1), TRUE);
        } else {
          ap = callpois(knew, alpha2/(q-1), TRUE) - callpois(knew, kappa_ibp*alpha2/(q-1), TRUE); 
        }
        arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
        FAugmented.col(j).fill(1);
        arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
        arma::mat lambdaUAugmented = FAugmented%VAugmented;
        
        psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
        
        arma::mat term_4 = psi_x_inv*lambdaUAugmented.t();
        
        arma::cube H_i = arma::zeros<arma::cube>(knew, knew, TruncateL); 
        arma::cube G_i_firstterm = arma::zeros<arma::cube>(q, knew, TruncateL); 
        arma::mat G_i_secondterm = arma::zeros<arma::mat>(TruncateL, knew); 
        
        arma::mat H_i_firstterm = lambdaUAugmented*term_4;
        arma::vec mu_sigma_mu = arma::zeros<arma::vec>(TruncateL);
        arma::cube Sigma_augmented_list = arma::zeros<arma::cube>(knew, knew, TruncateL);
        arma::mat Mu_augmented_list = arma::zeros<arma::mat>(TruncateL, knew);
        
        
        for(int l = 0; l < TruncateL; ++l){
          arma::mat knewI = arma::eye(knew, knew);
          Sigma_augmented_list.slice(l) = riwishart(knew + nu0, knewI); 
          Mu_augmented_list.row(l) = mvrnormArma(arma::zeros<arma::rowvec>(knew), Sigma_augmented_list.slice(l)/kappa0);
          if(any(C==l)){
            arma::rowvec Mu_augmented = Mu_augmented_list.row(l);
            arma::mat Sigma_augmented = Sigma_augmented_list.slice(l);
            arma::mat inv_Sigma_augmented = inv(Sigma_augmented);
            arma::mat term_2 = Mu_augmented*inv_Sigma_augmented;
            mu_sigma_mu(l) = as_scalar(term_2*Mu_augmented.t());
            H_i.slice(l) = H_i_firstterm + inv_Sigma_augmented;
            arma::mat term_3 = inv(H_i.slice(l));
            G_i_secondterm.row(l) = term_2*term_3;
            G_i_firstterm.slice(l) = term_4*term_3; 		     
          }	
        }
        
        E_x = X - Q*lambdaU;
        double GHG_sum = 0; 
        double mu_sigma_mu_sum = 0;
        double det_H = 0;
        
        for(int i = 0; i < n; ++i){ 
          arma::rowvec gt = E_x.row(i)*G_i_firstterm.slice(C(i)) + G_i_secondterm.row(C(i));
          GHG_sum = GHG_sum + as_scalar(gt*H_i.slice(C(i))*gt.t());
          mu_sigma_mu_sum = mu_sigma_mu_sum + mu_sigma_mu(C(i));
          det_H = det_H + std::log(det(H_i.slice(C(i)))) + log(det(Sigma_augmented_list.slice(C(i))));
        }					
        double lratio  = - 0.5*det_H + 0.5*(GHG_sum - mu_sigma_mu_sum) + ap;
        
        int kfinal = 0;  
        if(std::exp(lratio) > 1){
          kfinal = knew;
        }   
        if(kfinal != 0){  
          int newK = activeK + kfinal;
          recordF.rows(activeK, newK - 1) = FAugmented;
          recordV.rows(activeK, newK - 1) = VAugmented;
          recordlambdaU.rows(activeK, newK - 1) = lambdaUAugmented;		
          
          arma::mat Q_alter = arma::zeros<arma::mat>(n, knew);
          arma::mat term_5 = psi_x_inv*lambdaUAugmented.t();
          arma::mat D_inv_firstterm = lambdaUAugmented*term_5;
          
          for(int i = 0; i < n; ++i){ 
            arma::mat term_6 = inv(Sigma_augmented_list.slice(C(i)));
            arma::mat D_inv = inv(D_inv_firstterm + term_6);
            arma::mat sum_0 = E_x.row(i)*term_5 + Mu_augmented_list.row(C(i))*term_6;
            arma::rowvec W_i = sum_0*D_inv;
            Q_alter.row(i) = mvrnormArma(W_i, D_inv);
          }
          
          recordQ.cols(activeK, newK - 1) = Q_alter;
          Mu.cols(activeK, newK - 1) = Mu_augmented_list;               
          activeK = newK;
        }
      }         
    }
    //Rcpp::Rcout << "testf" << std::endl;
    arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
    arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
    arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat lambdaUnew = subsetMat(recordlambdaU.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat Qnew = subsetMat(recordQ.cols(0, activeK - 1), dseq, FALSE);
    
    activeK = Fnew.n_rows;
    
    arma::mat Munew = arma::zeros<arma::mat>(TruncateL, activeK);
    arma::cube Sigmanew = arma::zeros<arma::cube>(activeK, activeK, TruncateL); 
    
    //Rcpp::Rcout << "test!" << std::endl; 
    // Sample Mu_l, Sigma_l  
    int nl = 0;
    arma::vec currentL(TruncateL);
    currentL.zeros();
    for(int l = 0; l < TruncateL; ++l) {
      arma::mat sumTerm = arma::zeros<arma::mat>(activeK, activeK); 
      arma::rowvec muBar = arma::zeros<arma::rowvec>(activeK);
      if(any(C == l)){
        currentL(l) = 1;
        arma::uvec IDs = find(C == l);
        muBar = colmean(Qnew.rows(IDs));
        nl = IDs.n_elem;
        for(int ii = 0; ii < nl; ii++){
          arma::rowvec ll = Qnew.row(IDs(ii)) - muBar;
          sumTerm = sumTerm + ll.t()*ll;
        }  
      }  
      Sigmanew.slice(l) = riwishart(nu0 + nl, arma::eye(activeK, activeK) + sumTerm + kappa0*nl*muBar.t()*muBar/(kappa0 + nl));  
      Munew.row(l) = mvrnormArma(nl*muBar/(kappa0 + nl), Sigmanew.slice(l)/(kappa0 + nl)); //mu0 = 0
    }	
    
    // Sample V_l
    
    for(int l = 0; l < TruncateL; ++l) {
      arma::uvec ID1s = find(C == l);
      arma::uvec ID2s = find(C > l);
      if(any(C == l)){
        v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
      } else {
        v(l) = R::rbeta(1, alpha);
      }	
      if(l == 0){
        vProd_minus(l) = std::log(1-v(l)); 
        vProd(l) = std::log(v(l)); 
      } else {
        vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
        vProd(l) = vProd_minus(l-1) + std::log(v(l));
      }
    }
    
    // Sample Q_i 
    psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
    
    //Rcpp::Rcout << activeK << std::endl;
    arma::mat term_7 = psi_x_inv*lambdaUnew.t();
    for(int i = 0; i < n; ++i) {
      arma::mat inv_Sigmanew_i = inv(Sigmanew.slice(C(i)));
      arma::mat Dinv = inv(lambdaUnew*term_7 + inv_Sigmanew_i);
      arma::mat term_8 = X.row(i)*term_7 + Munew.row(C(i))*inv_Sigmanew_i;
      arma::rowvec Wi = term_8*Dinv;
      Qnew.row(i) = mvrnormArma(Wi, Dinv);
    }
    
    // Sample C_i
    for(int i = 0; i < n; ++i) {
      arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
      for(int l = 0; l < TruncateL; ++l) {
        likelihood(l) = vProd(l) + dmvnrmRowArma(Qnew.row(i), Munew.row(l), Sigmanew.slice(l), TRUE);
      }
      likelihood = likelihood + abs(max(likelihood));
      likelihood = exp(likelihood);
      prob = likelihood/sum(likelihood);
      C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
      
    }
    
    E_x = X - Qnew*lambdaUnew;
    
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
    }
    
    psi_x.diag() = psi_x_vector;
    h = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(lambdaUnew(arma::find(Fnew == 1)), 2))/2));	
    sigma = sqrt(1/sigma2inv);
    
    d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));		
    
    double vProd_minus_add = 0;
    double counter = 0;
    for(int l = 0; l < TruncateL; ++l) {
      if(currentL(l) == 1){
        vProd_minus_add = vProd_minus_add + std::log(1-v(l));
        counter = counter + 1;
      }
    }
    
    // Sample alpha
    alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
    
    // Sample alpha2 
    alpha2 = R::rgamma(1 + activeK, 1/(1 + har));
    
    recordF.rows(0, activeK-1) = Fnew;
    recordV.rows(0, activeK-1) = Vnew;
    recordlambdaU.rows(0, activeK-1) = lambdaUnew;   
    recordQ.cols(0, activeK-1) = Qnew; 
    Mu.cols(0, activeK-1) = Munew; 
    
    for(int l = 0; l < TruncateL; ++l){
      Sigma.slice(l).submat(arma::span(0, activeK-1), arma::span(0, activeK-1)) = Sigmanew.slice(l);
    }
    
    arma::vec kappa_likelihood = arma::zeros<arma::vec>(kappan);
    
    // Sample X_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < q; ++j){
        if(H(i, j) == 1){
          double term_ij = 2*kappa*psi_x_vector(j) + 1;
          double mu_ij = (as_scalar(Qnew.row(i)*lambdaUnew.col(j)) - 2*kappa*psi_x_vector(j)*tau(j))/term_ij ;
          double sigma_ij = psi_x_vector(j)/term_ij;
          X(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(X(i, j) + tau(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(X(i, j) + tau(j), 2)));
          }
        }
        
      }
    }
    
    kappa_likelihood = kappa_likelihood + abs(max(kappa_likelihood));
    kappa_likelihood = exp(kappa_likelihood);
    arma::vec kappa_prob = kappa_likelihood/sum(kappa_likelihood);
    kappa = as_scalar(RcppArmadillo::sample(kappagrid, 1, FALSE, as<NumericVector>(wrap(kappa_prob))));
    
    if(iter - kk <= iter_to_average){
      for(int i = 0; i < n; ++i) {
        AveC(C(i), i) = AveC(C(i), i) + 1;
      }
    }
  }    
  
  AveC = AveC/iter_to_average;   
  
  return Rcpp::List::create(
    Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
    Rcpp::Named("F") = recordF.rows(0, activeK-1),
    Rcpp::Named("C") = C,
    Rcpp::Named("AveC") = AveC,
    Rcpp::Named("lambdaU") = recordlambdaU.rows(0, activeK-1),
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("psi") = psi_x_vector,
    Rcpp::Named("kappa") = kappa); 
}


// Using a partial least square framework to model ERCC spike-in
// A mixture factor model on rows using a dirichlet process prior 
// with unknown number of clusters and a fixed number of factors
  
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
  
Rcpp::List DirichletSpikePLSModel(arma::mat Y, arma::mat X, int k1, int K, int iniL, int TruncateL, 
                                    int iter, int nu0, double sigma, double r, double s, double alpha, 
                                    arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, 
                                    double c, double d, double diagH, double h1, double h2, 
                                    double s1, double s2, int iter_to_average){
    
    // initialization 
    arma::mat U1;
    arma::vec S1;
    arma::mat V1;
    
    bool SVD1success = false;
    while(SVD1success == false) {
      SVD1success = svd(U1, S1, V1, Y);
      if(SVD1success == false){
        Y += 1e-4;
      }
    }
    
    int p = Y.n_cols;
    int n = X.n_rows; 
    int q = X.n_cols;
    
    arma::mat matd1 = arma::zeros<arma::mat>(k1, k1);
    matd1.diag() = sqrt(S1(arma::span(0, k1-1)));
    arma::mat Z_Y = U1.cols(0, k1-1)*matd1;
    arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
    arma::mat lambdaY = initialize_lambda(Z_Y_norm, Y);
    arma::mat lambdaX = initialize_lambda(Z_Y_norm, X); 
    arma::mat res_x = X - Z_Y_norm*lambdaX;
    
    arma::mat U;
    arma::vec S;
    arma::mat P;
    
    bool SVD2success = false;
    while(SVD2success == false) {
      SVD2success = svd(U, S, P, res_x);
      if(SVD2success == false){
        res_x += 1e-4;
      }
    }
    
    arma::mat matd = arma::zeros<arma::mat>(K, K);
    matd.diag() = sqrt(S(arma::span(0, K-1)));
    arma::mat Z_X = U.cols(0, K-1)*matd;
    arma::mat Q = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
    arma::vec C = callKmeans(Q, iniL) - 1; 
    arma::mat lamU = initialize_lambda(Q, res_x);
    
    arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
    arma::mat F = arma::conv_to<arma::mat>::from(a);
    arma::mat lambdaU = F%lamU;
    
    arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
    arma::mat Mu = arma::zeros<arma::mat>(TruncateL, K);  
    arma::cube Sigma = arma::zeros<arma::cube>(K, K, TruncateL); 
    arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
    
    for(int l = 0; l < iniL; ++l){
      arma::uvec IDs = find(C == l);
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(Q.rows(IDs));
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(Q.rows(IDs));
      seqx(l) = l;
    }
    
    arma::mat diagK = arma::eye(K, K);
    
    arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
    
    for(int l = iniL; l < TruncateL; ++l){
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0, diagK);
      seqx(l) = l;
    }
    
    arma::vec probF(2);
    probF.zeros();
    arma::vec seqF = arma::zeros<arma::vec>(2);
    seqF(1) = 1;
    
    double sigma2inv = 1/(sigma*sigma);
    
    arma::mat k1I = arma::eye(k1, k1);
    
    arma::mat psi_x = arma::eye(q, q);
    arma::mat psi_x_inv = psi_x;
    arma::vec psi_x_vector = psi_x.diag();
    
    arma::mat psi_y = arma::eye(p, p);
    arma::vec psi_y_vector = psi_y.diag();
    
    arma::mat psi = arma::eye(p+q, p+q);
    
    arma::mat H = arma::zeros<arma::mat>(k1, k1);
    H.diag().fill(diagH);
    
    arma::vec prob(TruncateL);
    prob.zeros();
    
    arma::mat Z = Z_Y_norm;
    arma::mat lambda = join_rows(lambdaY, lambdaX);
    arma::mat E_y = Y - Z*lambdaY;
    arma::mat E_x = X - Z*lambdaX - Q*lambdaU;
    
    arma::mat W = join_rows(Y, X - Q*lambdaU);
    arma::mat A = X - Z*lambdaX;
    
    arma::vec v = arma::zeros<arma::vec>(TruncateL);
    arma::vec rho = arma::zeros<arma::vec>(K);
    
    for(int i = 0; i < K; ++i){
      rho(i) = R::rbeta(r*s, s*(1-r));
    }
    
    arma::vec vProd = v;
    arma::vec vProd_minus = v;
    
    for (int kk = 0; kk < iter; ++kk) { 
      
      // Sample V_l
      
      for(int l = 0; l < TruncateL; ++l) {
        arma::uvec ID1s = find(C == l);
        arma::uvec ID2s = find(C > l);
        if(any(C == l)){
          v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
        } else {
          v(l) = R::rbeta(1, alpha);
        }	
        if(l == 0){
          vProd_minus(l) = std::log(1-v(l)); 
          vProd(l) = std::log(v(l)); 
        } else {
          vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
          vProd(l) = vProd_minus(l-1) + std::log(v(l));
        }
      }
      
      // Sample Mu_l, Sigma_l 
      
      int nl = 0;
      arma::vec currentL(TruncateL);
      currentL.zeros();
      for(int l = 0; l < TruncateL; ++l) {
        arma::mat sumTerm = arma::zeros<arma::mat>(K, K); 
        arma::rowvec muBar = arma::zeros<arma::rowvec>(K);
        if(any(C == l)){
          currentL(l) = 1;
          arma::uvec IDs = find(C == l);
          muBar = colmean(Q.rows(IDs));
          nl = IDs.n_elem;
          for(int ii = 0; ii < nl; ii++){
            arma::rowvec ll = Q.row(IDs(ii)) - muBar;
            sumTerm = sumTerm + ll.t()*ll;
          }
        } 		
        Sigma.slice(l) = riwishart(nu0 + nl, Sigma0 + sumTerm + kappa0*nl*(muBar - mu0).t()*(muBar - mu0)/(kappa0 + nl));	
        Mu.row(l) = mvrnormArma((kappa0*mu0 + nl*muBar)/(kappa0 + nl), Sigma.slice(l)/(kappa0 + nl));
      }	  
      
      // Sample C_i
      
      for(int i = 0; i < n; ++i) {
        arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
        for(int l = 0; l < TruncateL; ++l) {
          likelihood(l) = vProd(l) + dmvnrmRowArma(Q.row(i), Mu.row(l), Sigma.slice(l), TRUE);
        }
        likelihood = likelihood + abs(max(likelihood));
        likelihood = exp(likelihood);
        prob = likelihood/sum(likelihood);
        C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
      }
      
      // Sample Q_i 
      A = X - Z*lambdaX;
      psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
      
      arma::mat term_1 = psi_x_inv*lambdaU.t();
      
      for(int i = 0; i < n; ++i) {
        arma::mat inv_Sigma_i = inv(Sigma.slice(C(i)));
        arma::mat Dinv = inv(lambdaU*term_1 + inv_Sigma_i);
        arma::mat sum_0 = A.row(i)*term_1 + Mu.row(C(i))*inv_Sigma_i;
        arma::rowvec Wi = sum_0*Dinv;
        Q.row(i) = mvrnormArma(Wi, Dinv);
      }
      
      // Sample F    
      A = X - Z*lambdaX; 
      
      for(int j = 0; j < q; ++j){    
        
        arma::mat B = arma::zeros<arma::mat>(n, K);
        
        for(int k = 0; k < K; ++k){           
          
          for(int i = 0; i < n; ++i){          
            double QFVsum = 0;
            for(int l = 0; l < K; ++l){
              if(l != k){
                QFVsum = QFVsum + Q(i, l)*lambdaU(l, j);
              } 
            }  
            B(i, k) = A(i, j) - QFVsum;
          }
          
          double sumQB = 0;
          for(int i = 0; i < n; ++i){
            sumQB = sumQB + Q(i, k)*B(i, k);
          }
          
          double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
          double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
          double mukj = sumQB*sigmakj/psi_x_vector(j);
          double ratiop = std::log(m*rho(k)/(1 - m*rho(k)));
          double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + ratiop;
          probF(0) = 1/(std::exp(ratio)+1);
          probF(1) = 1 - probF(0); 
          F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
          if(F(k, j) == 1){
            lambdaU(k, j) = R::rnorm(mukj, sqrt(sigmakj));
          } else {
            lambdaU(k, j) = 0;
          }
        }  
      }
      
      // Sample rho_k 
      for(int k = 0; k < K; ++k){ 
        rho(k) = R::rbeta(r*s + sum(F.row(k)), s*(1-r) + q - sum(F.row(k)));        
      }    
      
      //  Sample Z 
      W = join_rows(Y, X - Q*lambdaU);
      arma::mat sharedterm2 = fastInverse(psi, lambda, TRUE);
      arma::mat term_2 = lambda*sharedterm2;
      arma::mat E_Z_w = (term_2*W.t()).t();
      arma::mat Var_Z_w = k1I - term_2*lambda.t();
      for(int i = 0; i < n; ++i){
        Z.row(i) = mvrnormArma(E_Z_w.row(i), Var_Z_w);
      }
      
      arma::mat sharedHterm = fastInverse(H, Z, TRUE);
      
      // Sample lambdaY 
      arma::mat term_3 = Y.t()*Z; 
      arma::mat term_4 = term_3*sharedHterm;
      lambdaY = sampleFromMND(k1, p, term_4.t(), sharedHterm, psi_y);
      
      // Sample lambda_x
      arma::mat term_5 = (X - Q*lambdaU).t();
      arma::mat term_6 = Z*sharedHterm;
      arma::mat term_7 = term_5*term_6;
      lambdaX = sampleFromMND(k1, q, term_7.t(), sharedHterm, psi_x);
      
      lambda = join_rows(lambdaY, lambdaX);
      
      // Sample psi_y, psi_x 
      E_y = Y - Z*lambdaY;
      for(int j = 0; j < p; ++j){
        psi_y_vector(j) = 1/(R::rgamma(g + n/2, 1/(h1 + sum(vectorise(E_y.col(j)%E_y.col(j)))/2)));
        if(psi_y_vector(j) < 0.000001){
          psi_y_vector(j) = 0.000001;
        }
      }
      psi_y.diag() = psi_y_vector;
      h1 = R::rgamma(1 + g*p, 1/(1 + sum(1/psi_y_vector)));
      
      E_x = X - Z*lambdaX - Q*lambdaU;
      //Rcpp::Rcout << sum(vectorise(E_x.t()*E_x)) << std::endl;
      for(int j = 0; j < q; ++j){
        psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h2 + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
        if(psi_x_vector(j) < 0.000001){
          psi_x_vector(j) = 0.000001;
        }
      }
      psi_x.diag() = psi_x_vector;
      h2 = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
      
      psi.diag() = join_cols(psi_y_vector, psi_x_vector);
      
      //Sample sigma
      sigma2inv = R::rgamma(c + sum(vectorise(F))/2, 1/(d + sum(pow(lambdaU(arma::find(F == 1)), 2))/2));	
      sigma = sqrt(1/sigma2inv);
      d = R::rgamma(1 + K, 1/(1 + sigma2inv*K));		
      
      double vProd_minus_add = 0;
      double counter = 0;
      for(int l = 0; l < TruncateL; ++l) {
        if(currentL(l) == 1){
          vProd_minus_add = vProd_minus_add + std::log(1-v(l));
          counter = counter + 1;
        }
      }
      
      // Sample alpha
      alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
      
      if(iter - kk <= iter_to_average){
        for(int i = 0; i < n; ++i) {
          AveC(C(i), i) = AveC(C(i), i) + 1;
        }
      }
    }    
    
    AveC = AveC/iter_to_average;
    
    return Rcpp::List::create(
      Rcpp::Named("Q") = Q,
      Rcpp::Named("F") = F,
      Rcpp::Named("C") = C,
      Rcpp::Named("AveC") = AveC,
      Rcpp::Named("Z") = Z,
      Rcpp::Named("lambdaX") = lambdaX,
      Rcpp::Named("lambdaY") = lambdaY,
      Rcpp::Named("lambdaU") = lambdaU, 
      Rcpp::Named("Mu") = Mu,
      Rcpp::Named("Sigma") = Sigma,
      Rcpp::Named("sigma") = sigma,
      Rcpp::Named("psi") = join_cols(psi_y_vector, psi_x_vector)); 
  }



// Using a partial least square framework to model ERCC spike-in
// A mixture factor model on rows using a dirichlet process prior 
// with unknown number of clusters and a fixed number of factors

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List DirichletSpikePLSModelZero(arma::mat Y, arma::mat X, arma::mat Hind, arma::mat Gind, arma::vec tau1, arma::vec tau2, 
                                  double kappa_int, arma::vec kappagrid, int k1, int K, int iniL, int TruncateL, 
                                  int iter, int nu0, double sigma, double r, double s, double alpha, 
                                  arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, 
                                  double c, double d, double diagH, double h1, double h2, 
                                  double s1, double s2, int iter_to_average){
  
  // initialization 
  arma::mat U1;
  arma::vec S1;
  arma::mat V1;
  
  bool SVD1success = false;
  while(SVD1success == false) {
    SVD1success = svd(U1, S1, V1, Y);
    if(SVD1success == false){
      Y += 1e-4;
    }
  }
  
  int p = Y.n_cols;
  int n = X.n_rows; 
  int q = X.n_cols;
  
  double kappa = kappa_int;
  int kappan = kappagrid.n_elem;
  
  arma::mat matd1 = arma::zeros<arma::mat>(k1, k1);
  matd1.diag() = sqrt(S1(arma::span(0, k1-1)));
  arma::mat Z_Y = U1.cols(0, k1-1)*matd1;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambdaY = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambdaX = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambdaX;
  
  arma::mat U;
  arma::vec S;
  arma::mat P;
  
  bool SVD2success = false;
  while(SVD2success == false) {
    SVD2success = svd(U, S, P, res_x);
    if(SVD2success == false){
      res_x += 1e-4;
    }
  }
  
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(S(arma::span(0, K-1)));
  arma::mat Z_X = U.cols(0, K-1)*matd;
  arma::mat Q = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
  arma::vec C = callKmeans(Q, iniL) - 1; 
  arma::mat lamU = initialize_lambda(Q, res_x);
  
  arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
  arma::mat F = arma::conv_to<arma::mat>::from(a);
  arma::mat lambdaU = F%lamU;
  
  arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n); 
  arma::mat Mu = arma::zeros<arma::mat>(TruncateL, K);  
  arma::cube Sigma = arma::zeros<arma::cube>(K, K, TruncateL); 
  arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
  
  for(int l = 0; l < iniL; ++l){
    arma::uvec IDs = find(C == l);
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(Q.rows(IDs));
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(Q.rows(IDs));
    seqx(l) = l;
  }
  
  arma::mat diagK = arma::eye(K, K);
  
  arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
  
  for(int l = iniL; l < TruncateL; ++l){
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0, diagK);
    seqx(l) = l;
  }
  
  arma::vec probF(2);
  probF.zeros();
  arma::vec seqF = arma::zeros<arma::vec>(2);
  seqF(1) = 1;
  
  double sigma2inv = 1/(sigma*sigma);
  
  arma::mat k1I = arma::eye(k1, k1);
  
  arma::mat psi_x = arma::eye(q, q);
  arma::mat psi_x_inv = psi_x;
  arma::vec psi_x_vector = psi_x.diag();
  
  arma::mat psi_y = arma::eye(p, p);
  arma::vec psi_y_vector = psi_y.diag();
  
  arma::mat psi = arma::eye(p+q, p+q);
  
  arma::mat H = arma::zeros<arma::mat>(k1, k1);
  H.diag().fill(diagH);
  
  arma::vec prob(TruncateL);
  prob.zeros();
  
  arma::mat Z = Z_Y_norm;
  arma::mat lambda = join_rows(lambdaY, lambdaX);
  arma::mat E_y = Y - Z*lambdaY;
  arma::mat E_x = X - Z*lambdaX - Q*lambdaU;
  
  arma::mat W = join_rows(Y, X - Q*lambdaU);
  arma::mat A = X - Z*lambdaX;
  
  arma::vec v = arma::zeros<arma::vec>(TruncateL);
  arma::vec rho = arma::zeros<arma::vec>(K);
  
  for(int i = 0; i < K; ++i){
    rho(i) = R::rbeta(r*s, s*(1-r));
  }
  
  arma::vec vProd = v;
  arma::vec vProd_minus = v;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    // Sample V_l
    
    for(int l = 0; l < TruncateL; ++l) {
      arma::uvec ID1s = find(C == l);
      arma::uvec ID2s = find(C > l);
      if(any(C == l)){
        v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
      } else {
        v(l) = R::rbeta(1, alpha);
      }	
      if(l == 0){
        vProd_minus(l) = std::log(1-v(l)); 
        vProd(l) = std::log(v(l)); 
      } else {
        vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
        vProd(l) = vProd_minus(l-1) + std::log(v(l));
      }
    }
    
    // Sample Mu_l, Sigma_l 
    
    int nl = 0;
    arma::vec currentL(TruncateL);
    currentL.zeros();
    for(int l = 0; l < TruncateL; ++l) {
      arma::mat sumTerm = arma::zeros<arma::mat>(K, K); 
      arma::rowvec muBar = arma::zeros<arma::rowvec>(K);
      if(any(C == l)){
        currentL(l) = 1;
        arma::uvec IDs = find(C == l);
        muBar = colmean(Q.rows(IDs));
        nl = IDs.n_elem;
        for(int ii = 0; ii < nl; ii++){
          arma::rowvec ll = Q.row(IDs(ii)) - muBar;
          sumTerm = sumTerm + ll.t()*ll;
        }
      } 		
      Sigma.slice(l) = riwishart(nu0 + nl, Sigma0 + sumTerm + kappa0*nl*(muBar - mu0).t()*(muBar - mu0)/(kappa0 + nl));	
      Mu.row(l) = mvrnormArma((kappa0*mu0 + nl*muBar)/(kappa0 + nl), Sigma.slice(l)/(kappa0 + nl));
    }	  
    
    // Sample C_i
    
    for(int i = 0; i < n; ++i) {
      arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
      for(int l = 0; l < TruncateL; ++l) {
        likelihood(l) = vProd(l) + dmvnrmRowArma(Q.row(i), Mu.row(l), Sigma.slice(l), TRUE);
      }
      likelihood = likelihood + abs(max(likelihood));
      likelihood = exp(likelihood);
      prob = likelihood/sum(likelihood);
      C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
    }
    
    // Sample Q_i 
    A = X - Z*lambdaX;
    psi_x_inv.diag() = 1/(psi_x_vector +0.000001);
    
    arma::mat term_1 = psi_x_inv*lambdaU.t();
    
    for(int i = 0; i < n; ++i) {
      arma::mat inv_Sigma_i = inv(Sigma.slice(C(i)));
      arma::mat Dinv = inv(lambdaU*term_1 + inv_Sigma_i);
      arma::mat sum_0 = A.row(i)*term_1 + Mu.row(C(i))*inv_Sigma_i;
      arma::rowvec Wi = sum_0*Dinv;
      Q.row(i) = mvrnormArma(Wi, Dinv);
    }
    
    // Sample F    
    A = X - Z*lambdaX; 
    
    for(int j = 0; j < q; ++j){    
      
      arma::mat B = arma::zeros<arma::mat>(n, K);
      
      for(int k = 0; k < K; ++k){           
        
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < K; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*lambdaU(l, j);
            } 
          }  
          B(i, k) = A(i, j) - QFVsum;
        }
        
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        double ratiop = std::log(m*rho(k)/(1 - m*rho(k)));
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + ratiop;
        probF(0) = 1/(std::exp(ratio)+1);
        probF(1) = 1 - probF(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
        if(F(k, j) == 1){
          lambdaU(k, j) = R::rnorm(mukj, sqrt(sigmakj));
        } else {
          lambdaU(k, j) = 0;
        }
      }  
    }
    
    // Sample rho_k 
    for(int k = 0; k < K; ++k){ 
      rho(k) = R::rbeta(r*s + sum(F.row(k)), s*(1-r) + q - sum(F.row(k)));        
    }    
    
    //  Sample Z 
    W = join_rows(Y, X - Q*lambdaU);
    arma::mat sharedterm2 = fastInverse(psi, lambda, TRUE);
    arma::mat term_2 = lambda*sharedterm2;
    arma::mat E_Z_w = (term_2*W.t()).t();
    arma::mat Var_Z_w = k1I - term_2*lambda.t();
    for(int i = 0; i < n; ++i){
      Z.row(i) = mvrnormArma(E_Z_w.row(i), Var_Z_w);
    }
    
    arma::mat sharedHterm = fastInverse(H, Z, TRUE);
    
    // Sample lambdaY 
    arma::mat term_3 = Y.t()*Z; 
    arma::mat term_4 = term_3*sharedHterm;
    lambdaY = sampleFromMND(k1, p, term_4.t(), sharedHterm, psi_y);
    
    // Sample lambda_x
    arma::mat term_5 = (X - Q*lambdaU).t();
    arma::mat term_6 = Z*sharedHterm;
    arma::mat term_7 = term_5*term_6;
    lambdaX = sampleFromMND(k1, q, term_7.t(), sharedHterm, psi_x);
    
    lambda = join_rows(lambdaY, lambdaX);
    
    // Sample psi_y, psi_x 
    E_y = Y - Z*lambdaY;
    for(int j = 0; j < p; ++j){
      psi_y_vector(j) = 1/(R::rgamma(g + n/2, 1/(h1 + sum(vectorise(E_y.col(j)%E_y.col(j)))/2)));
      if(psi_y_vector(j) < 0.000001){
        psi_y_vector(j) = 0.000001;
      }
    }
    psi_y.diag() = psi_y_vector;
    h1 = R::rgamma(1 + g*p, 1/(1 + sum(1/psi_y_vector)));
    
    E_x = X - Z*lambdaX - Q*lambdaU;
    //Rcpp::Rcout << sum(vectorise(E_x.t()*E_x)) << std::endl;
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h2 + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
    }
    psi_x.diag() = psi_x_vector;
    h2 = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    psi.diag() = join_cols(psi_y_vector, psi_x_vector);
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(F))/2, 1/(d + sum(pow(lambdaU(arma::find(F == 1)), 2))/2));	
    sigma = sqrt(1/sigma2inv);
    d = R::rgamma(1 + K, 1/(1 + sigma2inv*K));		
    
    double vProd_minus_add = 0;
    double counter = 0;
    for(int l = 0; l < TruncateL; ++l) {
      if(currentL(l) == 1){
        vProd_minus_add = vProd_minus_add + std::log(1-v(l));
        counter = counter + 1;
      }
    }
    
    // Sample alpha
    alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
    
    
    arma::vec kappa_likelihood = arma::zeros<arma::vec>(kappan);
    
    // Sample X_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < q; ++j){
        if(Hind(i, j) == 1){
          double term_ij = 2*kappa*psi_x_vector(j) + 1;
          double mu_ij = (as_scalar(Q.row(i)*lambdaU.col(j)) + as_scalar(Z.row(i)*lambdaX.col(j)) 
                            - 2*kappa*psi_x_vector(j)*tau1(j))/term_ij ;
          double sigma_ij = psi_x_vector(j)/term_ij;
          X(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(X(i, j) + tau1(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(X(i, j) + tau1(j), 2)));
          }
        }
        
      }
    }
    
    // Sample Y_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < p; ++j){
        if(Gind(i, j) == 1){
          double term_ij = 2*kappa*psi_y_vector(j) + 1;
          double mu_ij = (as_scalar(Z.row(i)*lambdaY.col(j)) 
                            - 2*kappa*psi_y_vector(j)*tau2(j))/term_ij ;
          double sigma_ij = psi_y_vector(j)/term_ij;
          Y(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(Y(i, j) + tau2(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(Y(i, j) + tau2(j), 2)));
          }
        }
        
      }
    }
    
    kappa_likelihood = kappa_likelihood + abs(max(kappa_likelihood));
    kappa_likelihood = exp(kappa_likelihood);
    arma::vec kappa_prob = kappa_likelihood/sum(kappa_likelihood);
    kappa = as_scalar(RcppArmadillo::sample(kappagrid, 1, FALSE, as<NumericVector>(wrap(kappa_prob))));
    
    if(iter - kk <= iter_to_average){
      for(int i = 0; i < n; ++i) {
        AveC(C(i), i) = AveC(C(i), i) + 1;
      }
    }
  }    
  
  AveC = AveC/iter_to_average;
  
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q,
    Rcpp::Named("F") = F,
    Rcpp::Named("C") = C,
    Rcpp::Named("AveC") = AveC,
    Rcpp::Named("Z") = Z,
    Rcpp::Named("lambdaX") = lambdaX,
    Rcpp::Named("lambdaY") = lambdaY,
    Rcpp::Named("lambdaU") = lambdaU, 
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("psi") = join_cols(psi_y_vector, psi_x_vector),
    Rcpp::Named("kappa") = kappa); 
}


// Using a partial least square framework to model ERCC spike-in
// a latent factor model using an IBP prior 

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSIBPfactormodel(arma::mat Y, arma::mat X, int k1, int k2, int iter, double sigma, double alpha, double kappa, double diagH, double g, double h1, double h2, double c, double d, int limit){
  
  Rcpp::List res = PLSinitilizationIBP(Y, X, k1, k2);
  arma::mat iniF = res["iniF"];
  arma::mat iniV = res["iniV"];
  arma::mat iniQ = res["iniQ"];
  arma::mat lambdaY = res["lambdaY"];
  arma::mat lambdaX = res["lambdaX"];
  arma::mat inilambdaU = res["lambdaU"];
  arma::mat Z = res["Z"];
  
  arma::mat lambda = join_rows(lambdaY, lambdaX);
  
  int activeK = k2; 
  int n = X.n_rows;
  int q = X.n_cols;
  int p = Y.n_cols;
  double har = harmonic(q);
  
  arma::mat recordF = arma::zeros<arma::mat>(limit, q);
  arma::mat recordV = arma::zeros<arma::mat>(limit, q);
  arma::mat recordlambdaU = arma::zeros<arma::mat>(limit, q);
  arma::mat recordQ = arma::zeros<arma::mat>(n, limit);
  
  recordF.rows(0, activeK-1) = iniF;
  recordV.rows(0, activeK-1) = iniV;
  recordlambdaU.rows(0, activeK-1) = inilambdaU;   
  recordQ.cols(0, activeK - 1) = iniQ;
  
  arma::mat k1I = arma::zeros<arma::mat>(k1, k1);
  k1I.eye(); 
  arma::mat psi_x = arma::zeros<arma::mat>(q, q);
  psi_x.eye();  
  arma::mat psi_y = arma::zeros<arma::mat>(p, p);
  psi_y.eye();
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.eye();  
  arma::mat H = arma::zeros<arma::mat>(k1, k1);
  H.diag().fill(diagH);
  
  double sigma2inv = 1/(sigma*sigma);
  arma::vec prob(2);
  prob.zeros();
  arma::vec seqx = arma::zeros<arma::vec>(2);
  seqx(1) = 1;
  
  arma::mat E_x = X - Z*lambdaX - iniQ*inilambdaU;
  arma::mat E_y = Y - Z*lambdaY;
  arma::vec psi_x_vector = psi_x.diag();
  arma::vec psi_y_vector = psi_y.diag();
  arma::mat W = join_rows(Y, X - iniQ*inilambdaU);
  arma::mat A = X - Z*lambdaX;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    A = X - Z*lambdaX;
    //Sample F;
    for(int j = 0; j < q; ++j){  
      
      arma::mat F = recordF.rows(0, activeK - 1);
      arma::mat V = recordV.rows(0, activeK - 1);
      arma::mat lambdaU = recordlambdaU.rows(0, activeK - 1);
      arma::mat Q = recordQ.cols(0, activeK - 1);
      arma::mat B = arma::zeros<arma::mat>(n, activeK);
      
      for(int k = 0; k < activeK; ++k){           
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < activeK; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
            } 
          }  
          B(i, k) = A(i, j) - QFVsum;
        }
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        arma::uvec allselect = find(F.row(k) == 1);
        int active = 0;
        if( F(k, j) == 1 ){
          active = allselect.n_elem - 1;
        } else {
          active = allselect.n_elem;
        }
        
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
        prob(0) = 1/(std::exp(ratio)+1);
        prob(1) = 1 - prob(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob)))); 
        
        if(F(k, j) == 1){
            V(k, j) = R::rnorm(mukj, sqrt(sigmakj));
        }  
        lambdaU(k, j) = V(k, j)*F(k, j);
        recordF(k, j) = F(k, j);
        recordV(k, j) = V(k, j);
        recordlambdaU(k, j) = lambdaU(k, j); 
      }  
    
        // Generate new features
      int knew = R::rpois(alpha/(q-1));
		    
      if(knew != 0){
        if(knew + activeK > limit) break;
        double ap = 0;
        if(kappa == 0){
          ap = - callpois(knew, alpha/(q-1), TRUE);
        } else {
          ap = callpois(knew, alpha/(q-1), TRUE) - callpois(knew, kappa*alpha/(q-1), TRUE); 
        }
        arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
        FAugmented.col(j).fill(1);
        arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
        arma::mat lambdaAugmented = FAugmented%VAugmented;
        arma::mat knewI = arma::zeros<arma::mat>(knew, knew);
        knewI.eye(); 
        arma::mat term_1 = inv(psi_x)*lambdaAugmented.t();
        arma::mat D = lambdaAugmented*term_1 + knewI;
        double CDCsum = 0; 
        arma::mat resQ = A - Q*lambdaU; 
        arma::mat tt = term_1*inv(D);
        for(int i = 0; i < n; ++i){
          arma::mat ct = resQ.row(i)*tt;
          CDCsum = CDCsum + as_scalar(ct*D*ct.t());
        }
        double lratio = - n/2*std::log(det(D)) + 0.5*CDCsum + ap;
                
        int kfinal = 0;  
        if(std::exp(lratio) > 1){
          kfinal = knew;
        }   
      	if(kfinal != 0){  
          recordF.rows(activeK, activeK + kfinal - 1) = FAugmented;
          recordV.rows(activeK, activeK + kfinal - 1) = VAugmented;
          recordlambdaU.rows(activeK, activeK + kfinal - 1) = lambdaAugmented;		
          arma::mat sharedterm = fastInverse(psi_x, lambdaAugmented, TRUE);
          arma::mat term_2 = lambdaAugmented*sharedterm;
          arma::mat E_Q_w = (term_2*resQ.t()).t();
          arma::mat Var_Q_w = knewI - term_2*lambdaAugmented.t();
		      arma::mat Qalter = arma::zeros<arma::mat>(n, knew);
          for(int i = 0; i < n; ++i){
					   Qalter.row(i) = mvrnormArma(E_Q_w.row(i), Var_Q_w);
				  }
          recordQ.cols(activeK, activeK + kfinal-1) = Qalter;
          activeK = activeK + kfinal;
      	}
      }         
    } 	 
    
            
    arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
    arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
    arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat lambdaUnew = subsetMat(recordlambdaU.rows(0, activeK - 1), dseq, TRUE); 
  
    activeK = Fnew.n_rows;

    // Sample Q 
    arma::mat EQ = X - Z*lambdaX;
    arma::mat sharedterm = fastInverse(psi_x, lambdaUnew, TRUE);
    arma::mat term_3 = lambdaUnew*sharedterm;
    arma::mat E_Q_w = (term_3*EQ.t()).t();
    arma::mat activeI = arma::zeros<arma::mat>(activeK, activeK);
    activeI.eye();
    arma::mat Var_Q_w = activeI - term_3*lambdaUnew.t();
  	arma::mat Qnew = arma::zeros<arma::mat>(n, activeK);
    for(int i = 0; i < n; ++i){
				Qnew.row(i) = mvrnormArma(E_Q_w.row(i), Var_Q_w);
		}
    
   //  Sample Z 
    W = join_rows(Y, X - Qnew*lambdaUnew);
	  arma::mat sharedterm2 = fastInverse(psi, lambda, TRUE);
    arma::mat term_4 = lambda*sharedterm2;
	  arma::mat E_Z_w = (term_4*W.t()).t();
    
    arma::mat Var_Z_w = k1I - term_4*lambda.t();
    for(int i = 0; i < k1; i++){ 
      for(int j = 0; j < k1; j++){ 
        if(Var_Z_w(i, j) < 0.000001){
          Var_Z_w(i, j) = 0.000001;
        }
      }
    }
    
    for(int i = 0; i < n; ++i){
				Z.row(i) = mvrnormArma(E_Z_w.row(i), Var_Z_w);
		}
    
    arma::mat sharedHterm = fastInverse(H, Z, TRUE);
    // Sample lambdaY 
    arma::mat term_5 = Y.t()*Z;
	  lambdaY = sampleFromMND(k1, p, (term_5*sharedHterm).t(), sharedHterm, psi_y);
	
    // Sample lambda_x
    arma::mat term_6 = (X - Qnew*lambdaUnew).t();
    arma::mat term_7 = Z*sharedHterm;
    arma::mat term_8 = (term_6*term_7).t();
	  lambdaX = sampleFromMND(k1, q, term_8, sharedHterm, psi_x);
	
	  lambda = join_rows(lambdaY, lambdaX);
     
    // Sample psi_y, psi_x 
    E_y = Y - Z*lambdaY;
	  for(int j = 0; j < p; ++j){
      psi_y_vector(j) = 1/(R::rgamma(g + n/2, 1/(h1 + sum(vectorise(E_y.col(j)%E_y.col(j)))/2)));
      if(psi_y_vector(j) < 0.000001){
        psi_y_vector(j) = 0.000001;
      }
	  }
	  psi_y.diag() = psi_y_vector;
	  h1 = R::rgamma(1 + g*p, 1/(1 + sum(1/psi_y_vector)));
    
    E_x = X - Z*lambdaX - Qnew*lambdaUnew;
    //Rcpp::Rcout << sum(vectorise(E_x.t()*E_x)) << std::endl;
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h2 + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
	  }
	  psi_x.diag() = psi_x_vector;
	  h2 = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    psi.diag() = join_cols(psi_y_vector, psi_x_vector);
  
    //Sample sigma
	  sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(Vnew(arma::find(Fnew == 1)), 2))/2));
	  sigma = sqrt(1/sigma2inv);
	  d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));
    
    //Sampling alpha 
   	alpha = R::rgamma(1 + activeK, 1/(1 + har));
     
    recordF.rows(0, activeK-1) = Fnew;
    recordV.rows(0, activeK-1) = Vnew;
    recordlambdaU.rows(0, activeK-1) = lambdaUnew;   
    recordQ.cols(0, activeK-1) = Qnew; 
  }  
  
  return Rcpp::List::create(
          Rcpp::Named("lambdaX") = lambdaX,
          Rcpp::Named("lambdaY") = lambdaY,
          Rcpp::Named("lambdaU") = recordlambdaU.rows(0, activeK-1), 
          Rcpp::Named("Z") =  Z,
          Rcpp::Named("F") = recordF.rows(0, activeK-1), 
          Rcpp::Named("V") = recordV.rows(0, activeK-1), 
          Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
          Rcpp::Named("psi") = join_cols(psi_y_vector, psi_x_vector),
          Rcpp::Named("activeK") = activeK);
}    

// Using a partial least square framework to model ERCC spike-in
// A mixture factor model using a dirichlet process prior on the number of clusters (mixture components)
// and an IBP prior on the number of latent factors.
  
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
  
Rcpp::List DirichletIBPPLSModel(arma::mat Y, arma::mat X, int k1, int K, int iniL, int TruncateL, int iter, 
                                  int maxK, int nu0, double sigma, double r, double s, double alpha, double alpha2, 
                                  arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, double c, double d, 
                                  double kappa, double diagH, double h1, double h2, double s1, double s2, int iter_to_average){
    // initialization 
    arma::mat U1;
    arma::vec S1;
    arma::mat V1;
    
    bool SVD1success = false;
    while(SVD1success == false) {
      SVD1success = svd(U1, S1, V1, Y);
      if(SVD1success == false){
        Y += 1e-4;
      }
    }
    
    int p = Y.n_cols;
    int n = X.n_rows; 
    int q = X.n_cols;

    double har = harmonic(q);
    
    arma::mat matd1 = arma::zeros<arma::mat>(k1, k1);
    matd1.diag() = sqrt(S1(arma::span(0, k1-1)));
    arma::mat Z_Y = U1.cols(0, k1-1)*matd1;
    arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
    arma::mat lambdaY = initialize_lambda(Z_Y_norm, Y);
    arma::mat lambdaX = initialize_lambda(Z_Y_norm, X); 
    arma::mat res_x = X - Z_Y_norm*lambdaX;
    
    arma::mat U;
    arma::vec S;
    arma::mat P;
    
    bool SVD2success = false;
    while(SVD2success == false) {
      SVD2success = svd(U, S, P, res_x);
      if(SVD2success == false){
        res_x += 1e-4;
      }
    }
    
    arma::mat matd = arma::zeros<arma::mat>(K, K);
    matd.diag() = sqrt(S(arma::span(0, K-1)));
    arma::mat Z_X = U.cols(0, K-1)*matd;
    arma::mat iniQ = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
    arma::vec C = callKmeans(iniQ, iniL) - 1; 
    arma::mat lamU = initialize_lambda(iniQ, res_x); 
    
    arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
    arma::mat iniF = arma::conv_to<arma::mat>::from(a);
    arma::mat inilambdaU = iniF%lamU;
    arma::mat iniV = inilambdaU;
    
    arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n);
    arma::mat Mu = arma::zeros<arma::mat>(TruncateL, maxK);  
    arma::cube Sigma = arma::zeros<arma::cube>(maxK, maxK, TruncateL); 
    arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
    
    for(int l = 0; l < iniL; ++l){
      arma::uvec IDs = find(C == l);
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(iniQ.rows(IDs));
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(iniQ.rows(IDs));
      seqx(l) = l;
    }
    
    arma::mat diagK = arma::eye(K, K);
    arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
    
    for(int l = iniL; l < TruncateL; ++l){
      Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
      Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
      seqx(l) = l;
    }
    
    arma::vec probF(2);
    probF.zeros();
    arma::vec seqF = arma::zeros<arma::vec>(2);
    seqF(1) = 1;
    
    double sigma2inv = 1/(sigma*sigma);
    
    arma::mat k1I = arma::eye(k1, k1);
    
    arma::mat psi_x = arma::eye(q, q);
    arma::mat psi_x_inv = psi_x;
    arma::vec psi_x_vector = psi_x.diag();
    
    arma::mat psi_y = arma::eye(p, p);
    arma::vec psi_y_vector = psi_y.diag();
    
    arma::mat psi = arma::eye(p+q, p+q);
    
    arma::mat H = arma::zeros<arma::mat>(k1, k1);
    H.diag().fill(diagH);
    
    arma::vec prob(TruncateL);
    prob.zeros();
    
    arma::mat Z = Z_Y_norm;
    arma::mat lambda = join_rows(lambdaY, lambdaX);
    arma::mat E_y = Y - Z*lambdaY;
    arma::mat E_x = X - Z*lambdaX - iniQ*inilambdaU;
    
    arma::mat W = join_rows(Y, X - iniQ*inilambdaU);
    arma::mat A = X - Z*lambdaX;
    
    arma::vec v = arma::zeros<arma::vec>(TruncateL);
    arma::vec vProd = v;
    arma::vec vProd_minus = v;
    
    arma::mat recordF = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordV = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordlambdaU = arma::zeros<arma::mat>(maxK, q);
    arma::mat recordQ = arma::zeros<arma::mat>(n, maxK);
    
    int activeK = K;
    recordF.rows(0, activeK-1) = iniF;
    recordV.rows(0, activeK-1) = iniV;
    recordlambdaU.rows(0, activeK-1) = inilambdaU;   
    recordQ.cols(0, activeK - 1) = iniQ;
    
    for (int kk = 0; kk < iter; ++kk) { 
      
      // Sample F    
      for(int j = 0; j < q; ++j){    
        
        arma::mat F = recordF.rows(0, activeK - 1);
        arma::mat V = recordV.rows(0, activeK - 1);
        arma::mat lambdaU = recordlambdaU.rows(0, activeK - 1);
        arma::mat Q = recordQ.cols(0, activeK - 1);
        arma::mat B = arma::zeros<arma::mat>(n, activeK);
        
        A = X - Z*lambdaX; 
        
        for(int k = 0; k < activeK; ++k){           
          
          for(int i = 0; i < n; ++i){          
            double QFVsum = 0;
            for(int l = 0; l < activeK; ++l){
              if(l != k){
                QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
              } 
            }  
            B(i, k) = A(i, j) - QFVsum;
          }
          
          double sumQB = 0;
          for(int i = 0; i < n; ++i){
            sumQB = sumQB + Q(i, k)*B(i, k);
          }
          
          double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
          double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
          double mukj = sumQB*sigmakj/psi_x_vector(j);
          arma::uvec allselect = find(F.row(k) == 1);
          int active = 0;
          if( F(k, j) == 1 ){
            active = allselect.n_elem - 1;
          } else {
            active = allselect.n_elem;
          }
          
          double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
          
          probF(0) = 1/(std::exp(ratio)+1);
          probF(1) = 1 - probF(0); 
          F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
          
          if(F(k, j) == 1){
            V(k, j)  = R::rnorm(mukj, sqrt(sigmakj));
          } 
          
          lambdaU(k, j) = V(k, j)*F(k, j);
          recordF(k, j) = F(k, j);
          recordV(k, j) = V(k, j);
          recordlambdaU(k, j) = lambdaU(k, j); 
        }  
        // Generate new features
        int knew = R::rpois(alpha2/(q-1));
        
        if(knew != 0){
          if(knew + activeK > maxK) break;
          double ap = 0;
          if(kappa == 0){
            ap = - callpois(knew, alpha2/(q-1), TRUE);
          } else {
            ap = callpois(knew, alpha2/(q-1), TRUE) - callpois(knew, kappa*alpha2/(q-1), TRUE); 
          }
          arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
          FAugmented.col(j).fill(1);
          arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
          arma::mat lambdaUAugmented = FAugmented%VAugmented;
          
          psi_x_inv.diag() = 1/(psi_x_vector + 0.000001);
          
          arma::mat term_10 = psi_x_inv*lambdaUAugmented.t();
          arma::cube H_i = arma::zeros<arma::cube>(knew, knew, TruncateL); 
          arma::cube G_i_firstterm = arma::zeros<arma::cube>(q, knew, TruncateL); 
          arma::mat G_i_secondterm = arma::zeros<arma::mat>(TruncateL, knew);
          
          arma::mat H_i_firstterm = lambdaUAugmented*term_10;
          arma::vec mu_sigma_mu = arma::zeros<arma::vec>(TruncateL);
          arma::cube Sigma_augmented_list = arma::zeros<arma::cube>(knew, knew, TruncateL);
          arma::mat Mu_augmented_list = arma::zeros<arma::mat>(TruncateL, knew);
          
          for(int l = 0; l < TruncateL; ++l){
            arma::mat knewI = arma::eye(knew, knew);
            Sigma_augmented_list.slice(l) = riwishart(knew + nu0, knewI); 
            Mu_augmented_list.row(l) = mvrnormArma(arma::zeros<arma::rowvec>(knew), Sigma_augmented_list.slice(l)/kappa0);
            if(any(C==l)){
              arma::rowvec Mu_augmented = Mu_augmented_list.row(l);
              arma::mat Sigma_augmented = Sigma_augmented_list.slice(l);
              arma::mat inv_Sigma_augmented = inv(Sigma_augmented);
              arma::mat term_9 = Mu_augmented*inv_Sigma_augmented;
              mu_sigma_mu(l) = as_scalar(term_9*Mu_augmented.t());
              H_i.slice(l) = H_i_firstterm + inv_Sigma_augmented;
              arma::mat inv_H_l = inv(H_i.slice(l));
              G_i_secondterm.row(l) = term_9*inv_H_l;
              G_i_firstterm.slice(l) = term_10*inv_H_l; 		     
            }	
          }
          
          E_x = X - Z*lambdaX - Q*lambdaU;
          double GHG_sum = 0; 
          double mu_sigma_mu_sum = 0;
          double det_H = 0;
          
          for(int i = 0; i < n; ++i){ 
            arma::rowvec gt = E_x.row(i)*G_i_firstterm.slice(C(i)) + G_i_secondterm.row(C(i));
            GHG_sum = GHG_sum + as_scalar(gt*H_i.slice(C(i))*gt.t());
            mu_sigma_mu_sum = mu_sigma_mu_sum + mu_sigma_mu(C(i));
            det_H = det_H + std::log(det(H_i.slice(C(i)))) + log(det(Sigma_augmented_list.slice(C(i))));
          }					
          double lratio  = - 0.5*det_H + 0.5*(GHG_sum - mu_sigma_mu_sum) + ap;
          
          int kfinal = 0;  
          if(std::exp(lratio) > 1){
            kfinal = knew;
          }   
          if(kfinal != 0){  
            int newK = activeK + kfinal;
            recordF.rows(activeK, newK - 1) = FAugmented;
            recordV.rows(activeK, newK - 1) = VAugmented;
            recordlambdaU.rows(activeK, newK - 1) = lambdaUAugmented;		
            
            arma::mat Q_alter = arma::zeros<arma::mat>(n, knew);
            arma::mat term_11 = psi_x_inv*lambdaUAugmented.t();
            arma::mat D_inv_firstterm = lambdaUAugmented*term_11;
            
            for(int i = 0; i < n; ++i){ 
              arma::mat inv_Sigma_augmented_i = inv(Sigma_augmented_list.slice(C(i)));
              arma::mat D_inv = inv(D_inv_firstterm + inv_Sigma_augmented_i);
              arma::mat term_12 = E_x.row(i)*term_11 + Mu_augmented_list.row(C(i))*inv_Sigma_augmented_i;
              arma::rowvec W_i = term_12*D_inv;
              Q_alter.row(i) = mvrnormArma(W_i, D_inv);
            }
            
            recordQ.cols(activeK, newK - 1) = Q_alter;
            Mu.cols(activeK, newK - 1) = Mu_augmented_list;               
            activeK = newK;
          }
        }         
      }
      //Rcpp::Rcout << "testf" << std::endl;
      arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
      arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
      arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
      arma::mat lambdaUnew = subsetMat(recordlambdaU.rows(0, activeK - 1), dseq, TRUE); 
      arma::mat Qnew = subsetMat(recordQ.cols(0, activeK - 1), dseq, FALSE);
      
      activeK = Fnew.n_rows;
      
      arma::mat Munew = arma::zeros<arma::mat>(TruncateL, activeK);
      arma::cube Sigmanew = arma::zeros<arma::cube>(activeK, activeK, TruncateL); 
      
      //Rcpp::Rcout << "test!" << std::endl; 
      // Sample Mu_l, Sigma_l  
      int nl = 0;
      arma::vec currentL(TruncateL);
      currentL.zeros();
      for(int l = 0; l < TruncateL; ++l) {
        arma::mat sumTerm = arma::zeros<arma::mat>(activeK, activeK); 
        arma::rowvec muBar = arma::zeros<arma::rowvec>(activeK);
        if(any(C == l)){
          currentL(l) = 1;
          arma::uvec IDs = find(C == l);
          muBar = colmean(Qnew.rows(IDs));
          nl = IDs.n_elem;
          for(int ii = 0; ii < nl; ii++){
            arma::rowvec ll = Qnew.row(IDs(ii)) - muBar;
            sumTerm = sumTerm + ll.t()*ll;
          }  
        } 
        Sigmanew.slice(l) = riwishart(nu0 + nl, arma::eye(activeK, activeK) + sumTerm + kappa0*nl*muBar.t()*muBar/(kappa0 + nl));  
        Munew.row(l) = mvrnormArma(nl*muBar/(kappa0 + nl), Sigmanew.slice(l)/(kappa0 + nl)); //mu0 = 0
      }	
      
      // Sample V_l
      
      for(int l = 0; l < TruncateL; ++l) {
        arma::uvec ID1s = find(C == l);
        arma::uvec ID2s = find(C > l);
        if(any(C == l)){
          v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
        } else {
          v(l) = R::rbeta(1, alpha);
        }	
        if(l == 0){
          vProd_minus(l) = std::log(1-v(l)); 
          vProd(l) = std::log(v(l)); 
        } else {
          vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
          vProd(l) = vProd_minus(l-1) + std::log(v(l));
        }
      }
      
      // Sample Q_i 
      psi_x_inv.diag() = 1/(psi_x_vector + 0.000001);
      
      //Rcpp::Rcout << activeK << std::endl;
      A = X - Z*lambdaX;
      
      arma::mat term_13 = psi_x_inv*lambdaUnew.t();
      for(int i = 0; i < n; ++i) {
        arma::mat inv_Sigmanew_i = inv(Sigmanew.slice(C(i))); 
        arma::mat Dinv = inv(lambdaUnew*term_13 + inv_Sigmanew_i);
        arma::mat term_14 = A.row(i)*term_13 + Munew.row(C(i))*inv_Sigmanew_i;
        arma::rowvec Wi = term_14*Dinv;
        Qnew.row(i) = mvrnormArma(Wi, Dinv);
      }
      
      // Sample C_i
      for(int i = 0; i < n; ++i) {
        arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
        for(int l = 0; l < TruncateL; ++l) {
          likelihood(l) = vProd(l) + dmvnrmRowArma(Qnew.row(i), Munew.row(l), Sigmanew.slice(l), TRUE);
        }
        likelihood = likelihood + abs(max(likelihood));
        likelihood = exp(likelihood);
        prob = likelihood/sum(likelihood);
        C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
        
      }
      
      //  Sample Z 
      W = join_rows(Y, X - Qnew*lambdaUnew);
      arma::mat sharedterm2 = fastInverse(psi, lambda, TRUE);
      arma::mat term_15 =lambda*sharedterm2 ;
      arma::mat E_Z_w = (term_15*W.t()).t();
      arma::mat Var_Z_w = k1I - term_15*lambda.t();
      for(int i = 0; i < n; ++i){
        Z.row(i) = mvrnormArma(E_Z_w.row(i), Var_Z_w);
      }
      
      arma::mat sharedHterm = fastInverse(H, Z, TRUE);
      
      // Sample lambdaY 
      arma::mat term_16 = Y.t()*Z;
      arma::mat term_16_2 = (term_16*sharedHterm).t();
      lambdaY = sampleFromMND(k1, p, term_16_2, sharedHterm, psi_y);
      
      // Sample lambda_x
      arma::mat term_17 = (X - Qnew*lambdaUnew).t();
      arma::mat term_18 = Z*sharedHterm; 
      arma::mat term_19 = (term_17*term_18).t();
      lambdaX = sampleFromMND(k1, q, term_19, sharedHterm, psi_x);
      
      lambda = join_rows(lambdaY, lambdaX);
      
      // Sample psi_y, psi_x 
      E_y = Y - Z*lambdaY;
      for(int j = 0; j < p; ++j){
        psi_y_vector(j) = 1/(R::rgamma(g + n/2, 1/(h1 + sum(vectorise(E_y.col(j)%E_y.col(j)))/2)));
        if(psi_y_vector(j) < 0.000001){
          psi_y_vector(j) = 0.000001;
        }
      }
      psi_y.diag() = psi_y_vector;
      h1 = R::rgamma(1 + g*p, 1/(1 + sum(1/psi_y_vector)));
      
      E_x = X - Z*lambdaX - Qnew*lambdaUnew;
      
      for(int j = 0; j < q; ++j){
        psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h2 + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
        if(psi_x_vector(j) < 0.000001){
          psi_x_vector(j) = 0.000001;
        }
      }
      psi_x.diag() = psi_x_vector;
      h2 = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
      
      psi.diag() = join_cols(psi_y_vector, psi_x_vector);
      
      //Sample sigma
      sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(lambdaUnew(arma::find(Fnew == 1)), 2))/2));	
      sigma = sqrt(1/sigma2inv);
      d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));			
      
      double vProd_minus_add = 0;
      double counter = 0;
      for(int l = 0; l < TruncateL; ++l) {
        if(currentL(l) == 1){
          vProd_minus_add = vProd_minus_add + std::log(1-v(l));
          counter = counter + 1;
        }
      }
      
      // Sample alpha
      alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
      
      // Sample alpha2 
      alpha2 = R::rgamma(1 + activeK, 1/(1 + har));
      
      //Rcpp::Rcout << "test12" << std::endl; 
      
      recordF.rows(0, activeK-1) = Fnew;
      recordV.rows(0, activeK-1) = Vnew;
      recordlambdaU.rows(0, activeK-1) = lambdaUnew;   
      recordQ.cols(0, activeK-1) = Qnew; 
      Mu.cols(0, activeK-1) = Munew; 
      
      for(int l = 0; l < TruncateL; ++l){
        Sigma.slice(l).submat(arma::span(0, activeK-1), arma::span(0, activeK-1)) = Sigmanew.slice(l);
      } 
      
      if(iter - kk <= iter_to_average){
        for(int i = 0; i < n; ++i) {
          AveC(C(i), i) = AveC(C(i), i) + 1;
        }
      }
    }    
    
    AveC = AveC/iter_to_average;      
    
    return Rcpp::List::create(
      Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
      Rcpp::Named("F") = recordF.rows(0, activeK-1),
      Rcpp::Named("C") = C,
      Rcpp::Named("AveC") = AveC,
      Rcpp::Named("Z") = Z,
      Rcpp::Named("lambdaX") = lambdaX,
      Rcpp::Named("lambdaY") = lambdaY,
      Rcpp::Named("lambdaU") = recordlambdaU.rows(0, activeK-1), 
      Rcpp::Named("Mu") = Mu,
      Rcpp::Named("Sigma") = Sigma,
      Rcpp::Named("sigma") = sigma,
      Rcpp::Named("psi") = join_cols(psi_y_vector, psi_x_vector)); 
  }


// Using a partial least square framework to model ERCC spike-in
// A mixture factor model using a dirichlet process prior on the number of clusters (mixture components)
// and an IBP prior on the number of latent factors. Adding the modeling of dropout.

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List DirichletIBPPLSModelZero(arma::mat Y, arma::mat X, arma::mat Hind, arma::mat Gind, arma::vec tau1, arma::vec tau2, 
                                    double kappa_int, arma::vec kappagrid, int k1, int K, int iniL, int TruncateL, int iter, 
                                int maxK, int nu0, double sigma, double r, double s, double alpha, double alpha2, 
                                arma::rowvec mu0, arma::mat Sigma0, int kappa0, double m, double g, double c, double d, 
                                double kappa_ibp, double diagH, double h1, double h2, double s1, double s2, int iter_to_average){
  // initialization 
  arma::mat U1;
  arma::vec S1;
  arma::mat V1;
  
  bool SVD1success = false;
  while(SVD1success == false) {
    SVD1success = svd(U1, S1, V1, Y);
    if(SVD1success == false){
      Y += 1e-4;
    }
  }
  
  int p = Y.n_cols;
  int n = X.n_rows; 
  int q = X.n_cols;
  
  double kappa = kappa_int;
  int kappan = kappagrid.n_elem;
  
  double har = harmonic(q);
  
  arma::mat matd1 = arma::zeros<arma::mat>(k1, k1);
  matd1.diag() = sqrt(S1(arma::span(0, k1-1)));
  arma::mat Z_Y = U1.cols(0, k1-1)*matd1;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambdaY = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambdaX = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambdaX;
  
  arma::mat U;
  arma::vec S;
  arma::mat P;
  
  bool SVD2success = false;
  while(SVD2success == false) {
    SVD2success = svd(U, S, P, res_x);
    if(SVD2success == false){
      res_x += 1e-4;
    }
  }
  
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(S(arma::span(0, K-1)));
  arma::mat Z_X = U.cols(0, K-1)*matd;
  arma::mat iniQ = (Z_X - mean(vectorise(Z_X)))/stddev(vectorise(Z_X));
  arma::vec C = callKmeans(iniQ, iniL) - 1; 
  arma::mat lamU = initialize_lambda(iniQ, res_x); 
  
  arma::umat a = abs(lamU) > Cquantile(vectorise(abs(lamU)), 0.3);   
  arma::mat iniF = arma::conv_to<arma::mat>::from(a);
  arma::mat inilambdaU = iniF%lamU;
  arma::mat iniV = inilambdaU;
  
  arma::mat AveC = arma::zeros<arma::mat>(TruncateL, n);
  arma::mat Mu = arma::zeros<arma::mat>(TruncateL, maxK);  
  arma::cube Sigma = arma::zeros<arma::cube>(maxK, maxK, TruncateL); 
  arma::vec seqx = arma::zeros<arma::vec>(TruncateL);
  
  for(int l = 0; l < iniL; ++l){
    arma::uvec IDs = find(C == l);
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = colmean(iniQ.rows(IDs));
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = cov(iniQ.rows(IDs));
    seqx(l) = l;
  }
  
  arma::mat diagK = arma::eye(K, K);
  arma::rowvec vecK = arma::zeros<arma::rowvec>(K);
  
  for(int l = iniL; l < TruncateL; ++l){
    Mu.submat(arma::span(l, l), arma::span(0, K-1)) = mvrnormArma(vecK, diagK*2);
    Sigma.slice(l).submat(arma::span(0, K-1), arma::span(0, K-1)) = riwishart(nu0 + K, diagK);
    seqx(l) = l;
  }
  
  arma::vec probF(2);
  probF.zeros();
  arma::vec seqF = arma::zeros<arma::vec>(2);
  seqF(1) = 1;
  
  double sigma2inv = 1/(sigma*sigma);
  
  arma::mat k1I = arma::eye(k1, k1);
  
  arma::mat psi_x = arma::eye(q, q);
  arma::mat psi_x_inv = psi_x;
  arma::vec psi_x_vector = psi_x.diag();
  
  arma::mat psi_y = arma::eye(p, p);
  arma::vec psi_y_vector = psi_y.diag();
  
  arma::mat psi = arma::eye(p+q, p+q);
  
  arma::mat H = arma::zeros<arma::mat>(k1, k1);
  H.diag().fill(diagH);
  
  arma::vec prob(TruncateL);
  prob.zeros();
  
  arma::mat Z = Z_Y_norm;
  arma::mat lambda = join_rows(lambdaY, lambdaX);
  arma::mat E_y = Y - Z*lambdaY;
  arma::mat E_x = X - Z*lambdaX - iniQ*inilambdaU;
  
  arma::mat W = join_rows(Y, X - iniQ*inilambdaU);
  arma::mat A = X - Z*lambdaX;
  
  arma::vec v = arma::zeros<arma::vec>(TruncateL);
  arma::vec vProd = v;
  arma::vec vProd_minus = v;
  
  arma::mat recordF = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordV = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordlambdaU = arma::zeros<arma::mat>(maxK, q);
  arma::mat recordQ = arma::zeros<arma::mat>(n, maxK);
  
  int activeK = K;
  recordF.rows(0, activeK-1) = iniF;
  recordV.rows(0, activeK-1) = iniV;
  recordlambdaU.rows(0, activeK-1) = inilambdaU;   
  recordQ.cols(0, activeK - 1) = iniQ;
  
  for (int kk = 0; kk < iter; ++kk) { 
    
    // Sample F    
    for(int j = 0; j < q; ++j){    
      
      arma::mat F = recordF.rows(0, activeK - 1);
      arma::mat V = recordV.rows(0, activeK - 1);
      arma::mat lambdaU = recordlambdaU.rows(0, activeK - 1);
      arma::mat Q = recordQ.cols(0, activeK - 1);
      arma::mat B = arma::zeros<arma::mat>(n, activeK);
      
      A = X - Z*lambdaX; 
      
      for(int k = 0; k < activeK; ++k){           
        
        for(int i = 0; i < n; ++i){          
          double QFVsum = 0;
          for(int l = 0; l < activeK; ++l){
            if(l != k){
              QFVsum = QFVsum + Q(i, l)*F(l, j)*V(l, j);
            } 
          }  
          B(i, k) = A(i, j) - QFVsum;
        }
        
        double sumQB = 0;
        for(int i = 0; i < n; ++i){
          sumQB = sumQB + Q(i, k)*B(i, k);
        }
        
        double sumQ2 = as_scalar(Q.col(k).t()*Q.col(k));  
        double sigmakj = 1/(sumQ2/psi_x_vector(j) + sigma2inv);
        double mukj = sumQB*sigmakj/psi_x_vector(j);
        arma::uvec allselect = find(F.row(k) == 1);
        int active = 0;
        if( F(k, j) == 1 ){
          active = allselect.n_elem - 1;
        } else {
          active = allselect.n_elem;
        }
        
        double ratio = 0.5*std::log(sigma2inv*sigmakj) + 0.5*mukj*mukj/sigmakj + logratio(active, q);
        
        probF(0) = 1/(std::exp(ratio)+1);
        probF(1) = 1 - probF(0); 
        F(k, j) = as_scalar(RcppArmadillo::sample(seqF, 1, FALSE, as<NumericVector>(wrap(probF)))); 
        
        if(F(k, j) == 1){
          V(k, j)  = R::rnorm(mukj, sqrt(sigmakj));
        } 
        
        lambdaU(k, j) = V(k, j)*F(k, j);
        recordF(k, j) = F(k, j);
        recordV(k, j) = V(k, j);
        recordlambdaU(k, j) = lambdaU(k, j); 
      }  
      // Generate new features
      int knew = R::rpois(alpha2/(q-1));
      
      if(knew != 0){
        if(knew + activeK > maxK) break;
        double ap = 0;
        if(kappa == 0){
          ap = - callpois(knew, alpha2/(q-1), TRUE);
        } else {
          ap = callpois(knew, alpha2/(q-1), TRUE) - callpois(knew, kappa*alpha2/(q-1), TRUE); 
        }
        arma::mat FAugmented = arma::zeros<arma::mat>(knew, q);
        FAugmented.col(j).fill(1);
        arma::mat VAugmented = arma::randn<arma::mat>(knew, q)*sigma;
        arma::mat lambdaUAugmented = FAugmented%VAugmented;
        
        psi_x_inv.diag() = 1/(psi_x_vector + 0.000001);
        
        arma::mat term_10 = psi_x_inv*lambdaUAugmented.t();
        arma::cube H_i = arma::zeros<arma::cube>(knew, knew, TruncateL); 
        arma::cube G_i_firstterm = arma::zeros<arma::cube>(q, knew, TruncateL); 
        arma::mat G_i_secondterm = arma::zeros<arma::mat>(TruncateL, knew);
        
        arma::mat H_i_firstterm = lambdaUAugmented*term_10;
        arma::vec mu_sigma_mu = arma::zeros<arma::vec>(TruncateL);
        arma::cube Sigma_augmented_list = arma::zeros<arma::cube>(knew, knew, TruncateL);
        arma::mat Mu_augmented_list = arma::zeros<arma::mat>(TruncateL, knew);
        
        for(int l = 0; l < TruncateL; ++l){
          arma::mat knewI = arma::eye(knew, knew);
          Sigma_augmented_list.slice(l) = riwishart(knew + nu0, knewI); 
          Mu_augmented_list.row(l) = mvrnormArma(arma::zeros<arma::rowvec>(knew), Sigma_augmented_list.slice(l)/kappa0);
          if(any(C==l)){
            arma::rowvec Mu_augmented = Mu_augmented_list.row(l);
            arma::mat Sigma_augmented = Sigma_augmented_list.slice(l);
            arma::mat inv_Sigma_augmented = inv(Sigma_augmented);
            arma::mat term_9 = Mu_augmented*inv_Sigma_augmented;
            mu_sigma_mu(l) = as_scalar(term_9*Mu_augmented.t());
            H_i.slice(l) = H_i_firstterm + inv_Sigma_augmented;
            arma::mat inv_H_l = inv(H_i.slice(l));
            G_i_secondterm.row(l) = term_9*inv_H_l;
            G_i_firstterm.slice(l) = term_10*inv_H_l; 		     
          }	
        }
        
        E_x = X - Z*lambdaX - Q*lambdaU;
        double GHG_sum = 0; 
        double mu_sigma_mu_sum = 0;
        double det_H = 0;
        
        for(int i = 0; i < n; ++i){ 
          arma::rowvec gt = E_x.row(i)*G_i_firstterm.slice(C(i)) + G_i_secondterm.row(C(i));
          GHG_sum = GHG_sum + as_scalar(gt*H_i.slice(C(i))*gt.t());
          mu_sigma_mu_sum = mu_sigma_mu_sum + mu_sigma_mu(C(i));
          det_H = det_H + std::log(det(H_i.slice(C(i)))) + log(det(Sigma_augmented_list.slice(C(i))));
        }					
        double lratio  = - 0.5*det_H + 0.5*(GHG_sum - mu_sigma_mu_sum) + ap;
        
        int kfinal = 0;  
        if(std::exp(lratio) > 1){
          kfinal = knew;
        }   
        if(kfinal != 0){  
          int newK = activeK + kfinal;
          recordF.rows(activeK, newK - 1) = FAugmented;
          recordV.rows(activeK, newK - 1) = VAugmented;
          recordlambdaU.rows(activeK, newK - 1) = lambdaUAugmented;		
          
          arma::mat Q_alter = arma::zeros<arma::mat>(n, knew);
          arma::mat term_11 = psi_x_inv*lambdaUAugmented.t();
          arma::mat D_inv_firstterm = lambdaUAugmented*term_11;
          
          for(int i = 0; i < n; ++i){ 
            arma::mat inv_Sigma_augmented_i = inv(Sigma_augmented_list.slice(C(i)));
            arma::mat D_inv = inv(D_inv_firstterm + inv_Sigma_augmented_i);
            arma::mat term_12 = E_x.row(i)*term_11 + Mu_augmented_list.row(C(i))*inv_Sigma_augmented_i;
            arma::rowvec W_i = term_12*D_inv;
            Q_alter.row(i) = mvrnormArma(W_i, D_inv);
          }
          
          recordQ.cols(activeK, newK - 1) = Q_alter;
          Mu.cols(activeK, newK - 1) = Mu_augmented_list;               
          activeK = newK;
        }
      }         
    }
    //Rcpp::Rcout << "testf" << std::endl;
    arma::uvec dseq = reorderMatF(recordF.rows(0, activeK - 1));  
    arma::mat Fnew = subsetMat(recordF.rows(0, activeK - 1), dseq, TRUE);  
    arma::mat Vnew = subsetMat(recordV.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat lambdaUnew = subsetMat(recordlambdaU.rows(0, activeK - 1), dseq, TRUE); 
    arma::mat Qnew = subsetMat(recordQ.cols(0, activeK - 1), dseq, FALSE);
    
    activeK = Fnew.n_rows;
    
    arma::mat Munew = arma::zeros<arma::mat>(TruncateL, activeK);
    arma::cube Sigmanew = arma::zeros<arma::cube>(activeK, activeK, TruncateL); 
    
    //Rcpp::Rcout << "test!" << std::endl; 
    // Sample Mu_l, Sigma_l  
    int nl = 0;
    arma::vec currentL(TruncateL);
    currentL.zeros();
    for(int l = 0; l < TruncateL; ++l) {
      arma::mat sumTerm = arma::zeros<arma::mat>(activeK, activeK); 
      arma::rowvec muBar = arma::zeros<arma::rowvec>(activeK);
      if(any(C == l)){
        currentL(l) = 1;
        arma::uvec IDs = find(C == l);
        muBar = colmean(Qnew.rows(IDs));
        nl = IDs.n_elem;
        for(int ii = 0; ii < nl; ii++){
          arma::rowvec ll = Qnew.row(IDs(ii)) - muBar;
          sumTerm = sumTerm + ll.t()*ll;
        }  
      } 
      Sigmanew.slice(l) = riwishart(nu0 + nl, arma::eye(activeK, activeK) + sumTerm + kappa0*nl*muBar.t()*muBar/(kappa0 + nl));  
      Munew.row(l) = mvrnormArma(nl*muBar/(kappa0 + nl), Sigmanew.slice(l)/(kappa0 + nl)); //mu0 = 0
    }	
    
    // Sample V_l
    
    for(int l = 0; l < TruncateL; ++l) {
      arma::uvec ID1s = find(C == l);
      arma::uvec ID2s = find(C > l);
      if(any(C == l)){
        v(l) = R::rbeta(1 + ID1s.n_elem, alpha + ID2s.n_elem);
      } else {
        v(l) = R::rbeta(1, alpha);
      }	
      if(l == 0){
        vProd_minus(l) = std::log(1-v(l)); 
        vProd(l) = std::log(v(l)); 
      } else {
        vProd_minus(l) = vProd_minus(l-1) + std::log(1-v(l));
        vProd(l) = vProd_minus(l-1) + std::log(v(l));
      }
    }
    
    // Sample Q_i 
    psi_x_inv.diag() = 1/(psi_x_vector + 0.000001);
    
    //Rcpp::Rcout << activeK << std::endl;
    A = X - Z*lambdaX;
    
    arma::mat term_13 = psi_x_inv*lambdaUnew.t();
    for(int i = 0; i < n; ++i) {
      arma::mat inv_Sigmanew_i = inv(Sigmanew.slice(C(i))); 
      arma::mat Dinv = inv(lambdaUnew*term_13 + inv_Sigmanew_i);
      arma::mat term_14 = A.row(i)*term_13 + Munew.row(C(i))*inv_Sigmanew_i;
      arma::rowvec Wi = term_14*Dinv;
      Qnew.row(i) = mvrnormArma(Wi, Dinv);
    }
    
    // Sample C_i
    for(int i = 0; i < n; ++i) {
      arma::vec likelihood = arma::zeros<arma::vec>(TruncateL);
      for(int l = 0; l < TruncateL; ++l) {
        likelihood(l) = vProd(l) + dmvnrmRowArma(Qnew.row(i), Munew.row(l), Sigmanew.slice(l), TRUE);
      }
      likelihood = likelihood + abs(max(likelihood));
      likelihood = exp(likelihood);
      prob = likelihood/sum(likelihood);
      C(i) = as_scalar(RcppArmadillo::sample(seqx, 1, FALSE, as<NumericVector>(wrap(prob))));
      
    }
    
    //  Sample Z 
    W = join_rows(Y, X - Qnew*lambdaUnew);
    arma::mat sharedterm2 = fastInverse(psi, lambda, TRUE);
    arma::mat term_15 =lambda*sharedterm2 ;
    arma::mat E_Z_w = (term_15*W.t()).t();
    arma::mat Var_Z_w = k1I - term_15*lambda.t();
    for(int i = 0; i < n; ++i){
      Z.row(i) = mvrnormArma(E_Z_w.row(i), Var_Z_w);
    }
    
    arma::mat sharedHterm = fastInverse(H, Z, TRUE);
    
    // Sample lambdaY 
    arma::mat term_16 = Y.t()*Z;
    arma::mat term_16_2 = (term_16*sharedHterm).t();
    lambdaY = sampleFromMND(k1, p, term_16_2, sharedHterm, psi_y);
    
    // Sample lambda_x
    arma::mat term_17 = (X - Qnew*lambdaUnew).t();
    arma::mat term_18 = Z*sharedHterm; 
    arma::mat term_19 = (term_17*term_18).t();
    lambdaX = sampleFromMND(k1, q, term_19, sharedHterm, psi_x);
    
    lambda = join_rows(lambdaY, lambdaX);
    
    // Sample psi_y, psi_x 
    E_y = Y - Z*lambdaY;
    for(int j = 0; j < p; ++j){
      psi_y_vector(j) = 1/(R::rgamma(g + n/2, 1/(h1 + sum(vectorise(E_y.col(j)%E_y.col(j)))/2)));
      if(psi_y_vector(j) < 0.000001){
        psi_y_vector(j) = 0.000001;
      }
    }
    psi_y.diag() = psi_y_vector;
    h1 = R::rgamma(1 + g*p, 1/(1 + sum(1/psi_y_vector)));
    
    E_x = X - Z*lambdaX - Qnew*lambdaUnew;
    
    for(int j = 0; j < q; ++j){
      psi_x_vector(j) = 1/(R::rgamma(g + n/2, 1/(h2 + sum(vectorise(E_x.col(j)%E_x.col(j)))/2)));
      if(psi_x_vector(j) < 0.000001){
        psi_x_vector(j) = 0.000001;
      }
    }
    psi_x.diag() = psi_x_vector;
    h2 = R::rgamma(1 + g*q, 1/(1 + sum(1/psi_x_vector)));
    
    psi.diag() = join_cols(psi_y_vector, psi_x_vector);
    
    //Sample sigma
    sigma2inv = R::rgamma(c + sum(vectorise(Fnew))/2, 1/(d + sum(pow(lambdaUnew(arma::find(Fnew == 1)), 2))/2));	
    sigma = sqrt(1/sigma2inv);
    d = R::rgamma(1 + activeK, 1/(1 + sigma2inv*activeK));			
    
    double vProd_minus_add = 0;
    double counter = 0;
    for(int l = 0; l < TruncateL; ++l) {
      if(currentL(l) == 1){
        vProd_minus_add = vProd_minus_add + std::log(1-v(l));
        counter = counter + 1;
      }
    }
    
    // Sample alpha
    alpha = R::rgamma(s1 + counter, 1/(s2 - vProd_minus_add));	
    
    // Sample alpha2 
    alpha2 = R::rgamma(1 + activeK, 1/(1 + har));
    
    //Rcpp::Rcout << "test12" << std::endl; 
    
    recordF.rows(0, activeK-1) = Fnew;
    recordV.rows(0, activeK-1) = Vnew;
    recordlambdaU.rows(0, activeK-1) = lambdaUnew;   
    recordQ.cols(0, activeK-1) = Qnew; 
    Mu.cols(0, activeK-1) = Munew; 
    
    for(int l = 0; l < TruncateL; ++l){
      Sigma.slice(l).submat(arma::span(0, activeK-1), arma::span(0, activeK-1)) = Sigmanew.slice(l);
    } 
    
    arma::vec kappa_likelihood = arma::zeros<arma::vec>(kappan);
    
    // Sample X_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < q; ++j){
        if(Hind(i, j) == 1){
          double term_ij = 2*kappa*psi_x_vector(j) + 1;
          double mu_ij = (as_scalar(Qnew.row(i)*lambdaUnew.col(j)) + as_scalar(Z.row(i)*lambdaX.col(j)) 
                            - 2*kappa*psi_x_vector(j)*tau1(j))/term_ij ;
          double sigma_ij = psi_x_vector(j)/term_ij;
          X(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(X(i, j) + tau1(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(X(i, j) + tau1(j), 2)));
          }
        }
        
      }
    }
    
    // Sample Y_ij
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < p; ++j){
        if(Gind(i, j) == 1){
          double term_ij = 2*kappa*psi_y_vector(j) + 1;
          double mu_ij = (as_scalar(Z.row(i)*lambdaY.col(j)) 
                            - 2*kappa*psi_y_vector(j)*tau2(j))/term_ij ;
          double sigma_ij = psi_y_vector(j)/term_ij;
          Y(i, j) = R::rnorm(mu_ij, sigma_ij); 
          
          // Sample kappa
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) - kappagrid(s)*pow(Y(i, j) + tau2(j), 2);
          }
        } else {
          for(int s = 0; s < kappan; ++s){
            kappa_likelihood(s) = kappa_likelihood(s) + std::log(1 - exp(- kappagrid(s)*pow(Y(i, j) + tau2(j), 2)));
          }
        }
        
      }
    }
    
    kappa_likelihood = kappa_likelihood + abs(max(kappa_likelihood));
    kappa_likelihood = exp(kappa_likelihood);
    arma::vec kappa_prob = kappa_likelihood/sum(kappa_likelihood);
    kappa = as_scalar(RcppArmadillo::sample(kappagrid, 1, FALSE, as<NumericVector>(wrap(kappa_prob))));
    
    
    if(iter - kk <= iter_to_average){
      for(int i = 0; i < n; ++i) {
        AveC(C(i), i) = AveC(C(i), i) + 1;
      }
    }
  }    
  
  AveC = AveC/iter_to_average;      
  
  return Rcpp::List::create(
    Rcpp::Named("Q") = recordQ.cols(0, activeK-1),
    Rcpp::Named("F") = recordF.rows(0, activeK-1),
    Rcpp::Named("C") = C,
    Rcpp::Named("AveC") = AveC,
    Rcpp::Named("Z") = Z,
    Rcpp::Named("lambdaX") = lambdaX,
    Rcpp::Named("lambdaY") = lambdaY,
    Rcpp::Named("lambdaU") = recordlambdaU.rows(0, activeK-1), 
    Rcpp::Named("Mu") = Mu,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("psi") = join_cols(psi_y_vector, psi_x_vector),
    Rcpp::Named("kappa") = kappa); 
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List factorEM(arma::mat X, arma::mat Y, int k,  int iter) {
  int p = Y.n_cols;
  int q = X.n_cols;
  int n = X.n_rows;
  arma::mat W = join_rows(Y, X); 
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.diag().ones();
  
  // initialization 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Y);
  
  arma::mat matd = arma::zeros<arma::mat>(k, k);
  matd.diag() = sqrt(s(arma::span(0, k-1)));
  arma::mat Z_Y = U.cols(0, k-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  
  arma::mat lambda = join_rows(lambda_y, lambda_x);
  
  arma::vec psi_vec = psi.diag();
  arma::mat psi_inv = psi;
  psi_inv.diag() = 1/psi_vec;
  
  arma::mat shared_term = arma::zeros<arma::mat>(p+q, p+q);
  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.diag().ones();
  arma::mat E_Z_w = arma::zeros<arma::mat>(p, k);
  
  for (int i=0; i < iter; i++) {
    
    arma::mat psi_term = arma::zeros<arma::mat>(p+q, p+q);
    arma::mat E_ZZ_w_all = arma::zeros<arma::mat>(k, k);
    arma::mat w_E_z_w_all = arma::zeros<arma::mat>(p+q, k); 
    psi_vec = psi.diag();
    psi_inv.diag() = 1/(psi_vec+0.00001);
    
    shared_term = fastInverse(psi, lambda, TRUE);    
    arma::mat term_4 = lambda*shared_term;
    E_Z_w = (term_4*W.t()).t();
    arma::mat Var_Z_w = I - term_4*lambda.t();
    
    for (int j=0; j < n; j++) {
      E_ZZ_w_all = E_ZZ_w_all + Var_Z_w + E_Z_w.row(j).t()*E_Z_w.row(j);
      w_E_z_w_all = w_E_z_w_all + W.row(j).t()*E_Z_w.row(j);
    }
    arma::mat lambda_1_new = w_E_z_w_all*inv(E_ZZ_w_all);
    arma::mat lambda_new = lambda_1_new.t();
    
    for (int j = 0; j < n; j++) {
      arma::mat term_5 = W.row(j).t() - lambda_new.t()*E_Z_w.row(j).t();
      psi_term =	psi_term + term_5*W.row(j);
    }
    
    if(sum(square(vectorise(abs(lambda.t()*lambda - lambda_new.t()*lambda_new)))) < 0.1) break;
    psi.diag() = psi_term.diag()/n;
    psi_vec = psi.diag();
    lambda = lambda_new;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda, 
                            Rcpp::Named("Z") = E_Z_w, 
                            Rcpp::Named("psi") = psi_vec);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEM(arma::mat X, arma::mat Y, int k1, int k2, int iter, double epsilon) {
  int p = Y.n_cols;
  int q = X.n_cols;
  int n = X.n_rows;
  int k = k1 + k2;
  arma::mat W = join_rows(Y, X); 
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.diag().ones();
  
  // initialization 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Y);
  
  arma::mat matd = arma::zeros<arma::mat>(k1, k1);
  matd.diag() = sqrt(s(arma::span(0, k1-1)));
  arma::mat Z_Y = U.cols(0, k1-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambda_x;
  
  arma::mat U2;
  arma::vec s2;
  arma::mat V2;
  svd(U2, s2, V2, res_x);
  
  arma::mat matd2 = arma::zeros<arma::mat>(k2, k2);
  matd2.diag() = sqrt(s2(arma::span(0, k2-1)));
  arma::mat Z_U = U2.cols(0, k2-1)*matd2;
  arma::mat Z_U_norm = (Z_U - mean(vectorise(Z_U)))/stddev(vectorise(Z_U));
  arma::mat lambda_u = initialize_lambda(Z_U_norm, res_x); 
  
  arma::mat lambda_zero = arma::zeros<arma::mat>(k2, p);
  arma::mat lambda_1 = join_rows(lambda_y, lambda_x);
  arma::mat lambda_2 = join_rows(lambda_zero, lambda_u);
  arma::mat lambda = join_cols(lambda_1, lambda_2);
  
  arma::vec psi_vec = psi.diag();
  arma::mat psi_inv = psi;
  psi_inv.diag() = 1/(psi_vec+0.00001);
  
  arma::mat shared_term = arma::zeros<arma::mat>(p+q, p+q);
  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.diag().ones();
  arma::mat E_Z_w = arma::zeros<arma::mat>(p, k);
  
  for (int i = 0; i < iter; i++) {
    
    arma::mat psi_term = arma::zeros<arma::mat>(p+q, p+q);
    arma::mat E_ZZ_w_all = arma::zeros<arma::mat>(k, k);
    arma::mat y_E_z_s_all = arma::zeros<arma::mat>(p, k1);
    arma::mat x_E_z_s_all = arma::zeros<arma::mat>(q, k1);
    arma::mat u_E_z_x_all = arma::zeros<arma::mat>(q, k2);
    
    psi_vec = psi.diag();
    psi_inv.diag() = 1/(psi_vec+0.00001);
    // Rcpp::Rcout << "Mark 2 " << psi_vec << std::endl;
    shared_term = fastInverse(psi, lambda, TRUE);
    //  Rcpp::Rcout << "Mark 2 "<< psi_vec << std::endl;      
    arma::mat term_4 = lambda*shared_term;
    E_Z_w = (term_4*W.t()).t();
    arma::mat Var_Z_w = I - term_4*lambda.t();
    
    for (int j = 0; j < n; j++) {
      E_ZZ_w_all = E_ZZ_w_all + Var_Z_w + E_Z_w.row(j).t()*E_Z_w.row(j);
      y_E_z_s_all = y_E_z_s_all + W.row(j).subvec(0, p-1).t()*E_Z_w.row(j).subvec(0, k1-1);
      x_E_z_s_all = x_E_z_s_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(0, k1-1);
      u_E_z_x_all = u_E_z_x_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(k1, k-1);
    }
    
    arma::mat E_Z_s_Z_s_all = E_ZZ_w_all.submat(0, 0, k1-1, k1-1);
    arma::mat E_Z_x_Z_s_all = E_ZZ_w_all.submat(k1, 0, k-1, k1-1);
    arma::mat E_Z_s_Z_x_all = E_ZZ_w_all.submat(0, k1, k1-1, k-1);
    arma::mat E_Z_x_Z_x_all = E_ZZ_w_all.submat(k1, k1, k-1, k-1);
    
    arma::mat lambda_y_new = y_E_z_s_all*inv(E_Z_s_Z_s_all);
    arma::mat lambda_x_new = (x_E_z_s_all - lambda_u.t()*E_Z_x_Z_s_all)*inv(E_Z_s_Z_s_all);
    arma::mat lambda_u_new = (u_E_z_x_all - lambda_x_new*E_Z_s_Z_x_all)*inv(E_Z_x_Z_x_all);
    arma::mat lambda_1_new = join_rows(lambda_y_new.t(), lambda_x_new.t());
    arma::mat lambda_2_new = join_rows(lambda_zero, lambda_u_new.t());
    arma::mat lambda_new = join_cols(lambda_1_new, lambda_2_new);
    
    for (int j = 0; j < n; j++) {
      arma::mat term_5 = W.row(j).t() - lambda_new.t()*E_Z_w.row(j).t();
      psi_term =	psi_term + term_5*W.row(j);
    }
    if(sum(square(vectorise(abs(lambda.t()*lambda - lambda_new.t()*lambda_new)))) <= epsilon) break;
    psi.diag() = psi_term.diag()/n;
    lambda = lambda_new;
    lambda_u = lambda_u_new.t();
  }
  
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda, Rcpp::Named("Z") = E_Z_w, Rcpp::Named("psi") = psi_vec);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List NaivePCA(arma::mat X, arma::mat Y, int k1, int k2) {
  int p = Y.n_cols;
  int q = X.n_cols;
  arma::mat W = join_rows(Y, X); 
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.diag().ones();
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Y);
  
  arma::mat matd = arma::zeros<arma::mat>(k1, k1);
  matd.diag() = sqrt(s(arma::span(0, k1-1)));
  arma::mat Z_Y = U.cols(0, k1-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambda_x;
  
  arma::mat U2;
  arma::vec s2;
  arma::mat V2;
  svd(U2, s2, V2, res_x);
  
  arma::mat matd2 = arma::zeros<arma::mat>(k2, k2);
  matd2.diag() = sqrt(s2(arma::span(0, k2-1)));
  arma::mat Z_U = U2.cols(0, k2-1)*matd2;
  arma::mat Z_U_norm = (Z_U - mean(vectorise(Z_U)))/stddev(vectorise(Z_U));
  arma::mat lambda_u = inv(Z_U_norm.t()*Z_U_norm)*Z_U_norm.t()*res_x; 
  
  arma::mat lambda_zero = arma::zeros<arma::mat>(k2, p);
  arma::mat lambda_1 = join_rows(lambda_y, lambda_x);
  arma::mat lambda_2 = join_rows(lambda_zero, lambda_u);
  arma::mat lambda = join_cols(lambda_1, lambda_2);
  arma::mat term_1 = W*lambda.t();
  arma::mat term_2 = inv(lambda*lambda.t());
  arma::mat Z = term_1*term_2;
  
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda, Rcpp::Named("Z") = Z);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSinitilization(arma::mat Y, arma::mat X, int k1, int k2) {
  // initialization 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Y);
  
  arma::mat matd = arma::zeros<arma::mat>(k1, k1);
  matd.diag() = sqrt(s(arma::span(0, k1-1)));
  arma::mat Z_Y = U.cols(0, k1-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambda_x;
  
  arma::mat U2;
  arma::vec s2;
  arma::mat V2;
  svd(U2, s2, V2, res_x);
  
  arma::mat matd2 = arma::zeros<arma::mat>(k2, k2);
  matd2.diag() = sqrt(s2(arma::span(0, k2-1)));
  arma::mat Z_U = U2.cols(0, k2-1)*matd2;
  arma::mat Z_U_norm = (Z_U - mean(vectorise(Z_U)))/stddev(vectorise(Z_U));
  arma::mat lambda_u = initialize_lambda(Z_U_norm, res_x); 
  
  return Rcpp::List::create(Rcpp::Named("lambdaY") = lambda_y, Rcpp::Named("lambdaX") = lambda_x, Rcpp::Named("lambdaU") = lambda_u);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec callLasso(arma::mat Y, arma::mat X, arma::mat D, double lam){
  Environment myEnv("package:Citrus");
  Function myLasso = myEnv["LassoByGenlasso"];
  return as<NumericVector>(wrap(myLasso(Y, X, D, lam)));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMpenalty(arma::mat Y, arma::mat X, int k1, int k2, int iter, arma::mat lambda_y, arma::mat lambda_x, arma::mat lambda_u, double lam, double tol) {
  int p = Y.n_cols;
  int q = X.n_cols;
  int n = X.n_rows;
  int k = k1 + k2;
  arma::mat W = join_rows(Y, X); 
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.diag().ones();
  
  arma::mat lambda_zero = arma::zeros<arma::mat>(k2, p);
  arma::mat lambda_1 = join_rows(lambda_y, lambda_x);
  arma::mat lambda_2 = join_rows(lambda_zero, lambda_u);
  arma::mat lambda = join_cols(lambda_1, lambda_2);
  
  arma::vec psi_vec = psi.diag();
  arma::mat psi_inv = psi;
  psi_inv.diag() = 1/(psi_vec+0.00001);
  
  arma::mat shared_term = arma::zeros<arma::mat>(p+q, p+q);
  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.diag().ones();
  arma::mat E_Z_w = arma::zeros<arma::mat>(p, k);
  
  for (int i=0; i < iter; i++) {
    
    arma::mat psi_term = arma::zeros<arma::mat>(p+q, p+q);
    arma::mat E_ZZ_w_all = arma::zeros<arma::mat>(k, k);
    arma::mat y_E_z_s_all = arma::zeros<arma::mat>(p, k1);
    arma::mat x_E_z_s_all = arma::zeros<arma::mat>(q, k1);
    arma::mat u_E_z_x_all = arma::zeros<arma::mat>(q, k2);
    
    psi_vec = psi.diag();
    psi_inv.diag() = 1/(psi_vec+0.00001);
    
    shared_term = fastInverse(psi, lambda, TRUE);    
    arma::mat term_4 = lambda*shared_term;
    E_Z_w = (term_4*W.t()).t();
    arma::mat Var_Z_w = I - term_4*lambda.t();
    
    for (int j=0; j < n; j++) {
      E_ZZ_w_all = E_ZZ_w_all + Var_Z_w + E_Z_w.row(j).t()*E_Z_w.row(j);
      y_E_z_s_all = y_E_z_s_all + W.row(j).subvec(0, p-1).t()*E_Z_w.row(j).subvec(0, k1-1);
      x_E_z_s_all = x_E_z_s_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(0, k1-1);
      u_E_z_x_all = u_E_z_x_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(k1, k-1);
    }
    
    arma::mat E_Z_s_Z_s_all = E_ZZ_w_all.submat(0, 0, k1-1, k1-1);
    arma::mat E_Z_x_Z_s_all = E_ZZ_w_all.submat(k1, 0, k-1, k1-1);
    arma::mat E_Z_s_Z_x_all = E_ZZ_w_all.submat(0, k1, k1-1, k-1);
    arma::mat E_Z_x_Z_x_all = E_ZZ_w_all.submat(k1, k1, k-1, k-1);
    
    arma::mat lambda_y_new = y_E_z_s_all*inv(E_Z_s_Z_s_all);
    arma::mat lambda_x_new = (x_E_z_s_all - lambda_u.t()*E_Z_x_Z_s_all)*inv(E_Z_s_Z_s_all);
    arma::mat lambda_u_new = arma::zeros<arma::mat>(k2, q);
    
    for(int j=0; j < q; j++){
      arma::mat matA = (u_E_z_x_all.row(j) - lambda_x_new.row(j)*E_Z_s_Z_x_all);
      arma::mat matC = E_Z_x_Z_x_all;
      //arma::mat matA = (u_E_z_x_all.row(j) - lambda_x_new.row(j)*E_Z_s_Z_x_all)/n;
      //arma::mat matC = E_Z_x_Z_x_all/n;
      arma::mat Xtilde = arma::chol(matC); //SRoot
      arma::mat Ytilde = (matA*inv(Xtilde)).t();
      arma::mat D = arma::zeros<arma::mat>(k2, k2);
      if(psi_vec(j) < 0.00001){
        psi_vec(j) = 0.00001;
      }
      D.diag().fill(psi_vec(j)/2);
      //Rcpp::Rcout << "Mark 2 "<< psi_vec(j) << std::endl;
      //D.diag().fill(psi_vec(j)/(2*n));
      
      lambda_u_new.col(j) = callLasso(Ytilde, Xtilde, D, lam);
    }
    
    arma::mat lambda_1_new = join_rows(lambda_y_new.t(), lambda_x_new.t());
    arma::mat lambda_2_new = join_rows(lambda_zero, lambda_u_new);
    arma::mat lambda_new = join_cols(lambda_1_new, lambda_2_new);
    for (int j=0; j < n; j++) {
      arma::mat term_5 = W.row(j).t() - lambda_new.t()*E_Z_w.row(j).t();
      psi_term =	psi_term + term_5*W.row(j);
    }
    
    if(sum(square(vectorise(abs(lambda.t()*lambda - lambda_new.t()*lambda_new)))) < tol) break;
    psi.diag() = psi_term.diag()/n;
    lambda = lambda_new;
    lambda_u = lambda_u_new;
    lambda_y = lambda_y_new.t();
    lambda_x = lambda_x_new.t();
  }
  return Rcpp::List::create(Rcpp::Named("lambdaY") = lambda_y, 
                            Rcpp::Named("lambdaX") = lambda_x, 
                            Rcpp::Named("lambdaU") = lambda_u, 
                            Rcpp::Named("lambda") = lambda, 
                            Rcpp::Named("Z") = E_Z_w, 
                            Rcpp::Named("psi") = psi_vec,
                            Rcpp::Named("psiMat") = psi, 
                            Rcpp::Named("pen") = lam);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

float floss(arma::mat X, arma::mat Y){
  return sum(pow(vectorise(X - Y), 2)); 
}
    
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List InNotIn(arma::vec x, int n) {
  arma::vec out = arma::zeros<arma::vec>(n);
  for (int i = 0; i < n; ++i) {
    arma::uvec found = find(x == i);
    if (found.n_elem == 0) {
      out(i) = 1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("NotIn") = find(out == 1), 
                             Rcpp::Named("In") = find(out == 0)); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMpenaltyTrain(arma::mat Y, arma::mat X, int k1, int k2, int iter, arma::vec penalty, double tol) {
  
  int n = Y.n_rows;
  int trainingsize = floor(2*n/3);
  
  arma::vec seqx = arma::zeros<arma::vec>(n);
  for (int i = 0; i < n; ++i) {
    seqx(i) = i;
  }
  
  arma::vec training = RcppArmadillo::sample(seqx, trainingsize, FALSE); 
  Rcpp::List allID = InNotIn(training, n);
  arma::uvec trainingID = allID["In"];
  arma::uvec tuneID = allID["NotIn"];
	arma::mat trainingX = X.rows(trainingID);
	arma::mat trainingY = Y.rows(trainingID);
	arma::mat tuneX = X.rows(tuneID);
	arma::mat tuneY = Y.rows(tuneID); 
   
	Rcpp::List res = PLSinitilization(trainingY, trainingX, k1, k2);
	arma::mat lambda_y = res["lambdaY"];
	arma::mat lambda_x = res["lambdaX"];
	arma::mat lambda_u = res["lambdaU"];
	  
  arma::mat tuneW = join_rows(tuneY, tuneX); 
	arma::mat covtuneW = cov(tuneW);
	
  arma::vec flossList = arma::zeros<arma::vec>(penalty.n_elem);
  
  for(int j=0; j < penalty.n_elem; j++){
    double pen = penalty(j);
    Rcpp::List est = PLSfactorEMpenalty(trainingY, trainingX, k1, k2, iter, lambda_y, lambda_x, lambda_u, pen, tol);
    arma::mat matL = est["lambda"];
    arma::mat matpsi = est["psiMat"];
		arma::mat estCov = matL.t()*matL + matpsi;
		flossList(j) = floss(covtuneW, estCov);  
	}
  arma::vec pen = penalty.elem(arma::find(flossList == flossList.min()));
  Rcpp::List finalEst = PLSfactorEMpenalty(Y, X, k1, k2, iter, lambda_y, lambda_x, lambda_u, pen(0), tol);
	return finalEst;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec CountSparsity(arma::mat X){
  int n = X.n_rows;
  arma::vec sparisity = arma::zeros<arma::vec>(n); 
  for(int j = 0; j < n; j++){
    sparisity(j) = sum(X.row(j));
  }
  return sparisity;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMpenaltyCVonefold(arma::mat trainingY, arma::mat trainingX, arma::mat tuneY, arma::mat tuneX, int k1, int k2, int iter, arma::vec penalty, double tol) {
  
	Rcpp::List res = PLSinitilization(trainingY, trainingX, k1, k2);
	arma::mat lambda_y = res["lambdaY"];
	arma::mat lambda_x = res["lambdaX"];
	arma::mat lambda_u = res["lambdaU"];
	
  arma::mat tuneW = join_rows(tuneY, tuneX); 
	arma::mat covtuneW = cov(tuneW);
	
  arma::vec flossList = arma::zeros<arma::vec>(penalty.n_elem);
  arma::vec sparsityList = arma::zeros<arma::vec>(penalty.n_elem);
  
  for(int j = 0; j < penalty.n_elem; j++){
    double pen = penalty(j);
    Rcpp::List est = PLSfactorEMpenalty(trainingY, trainingX, k1, k2, iter, lambda_y, lambda_x, lambda_u, pen, tol);
    arma::mat matL = est["lambda"];
    arma::mat matpsi = est["psiMat"];
		arma::mat estCov = matL.t()*matL + matpsi;
    flossList(j) = floss(covtuneW, estCov);
    arma::mat estLambdaU = est["lambdaU"];
    arma::uvec sparseRows = arma::find(CountSparsity(estLambdaU) < 0.000001);
		sparsityList(j) = sparseRows.n_elem;  
	}
  return Rcpp::List::create(Rcpp::Named("Floss") = flossList, 
                             Rcpp::Named("Sparsity") = sparsityList); 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMpenaltyNfoldCV(arma::mat Y, arma::mat X, int k1, int k2, int iter, int kfold, arma::vec penalty, double tol){
  
	int n = Y.n_rows;
	int onePart = floor(n/kfold);
  
  arma::vec seqx = arma::zeros<arma::vec>(n);
  for (int i = 0; i < n; ++i) {
    seqx(i) = i;
  }
  
	arma::vec norder = RcppArmadillo::sample(seqx, n, FALSE);
  arma::vec bestFloss = arma::zeros<arma::vec>(kfold);
  arma::vec bestPen = arma::zeros<arma::vec>(kfold);
  arma::vec bestSparsity = arma::zeros<arma::vec>(kfold);
   
  for (int i = 0; i < kfold; ++i) {
    arma::vec sets = norder.subvec(0 + i*onePart, onePart - 1 + i*onePart);
    Rcpp::List allID = InNotIn(sets, n);
    arma::uvec trainingID = allID["NotIn"];
    arma::uvec tuneID = allID["In"];
    arma::mat trainingX = X.rows(trainingID);
	  arma::mat trainingY = Y.rows(trainingID);
	  arma::mat tuneX = X.rows(tuneID);
	  arma::mat tuneY = Y.rows(tuneID); 
    Rcpp::List est = PLSfactorEMpenaltyCVonefold(trainingY, trainingX, tuneY, tuneX, k1, k2, iter, penalty, tol);
    arma::vec Floss = est["Floss"];
  	arma::vec Sparsity = est["Sparsity"];
    bestSparsity(i) = Sparsity.min();
    bestFloss(i) = Floss.elem(arma::find(Sparsity == bestSparsity(i))).min();
    arma::vec penLoss = penalty.elem(arma::find(Floss == bestFloss(i)));
    bestPen(i) = penLoss(0);
  }    
    
  arma::vec pen = arma::zeros<arma::vec>(1);  
  arma::uvec ID = arma::find(bestSparsity == bestSparsity.min());
  if(ID.n_elem == 1){
    pen = bestPen.elem(ID);
  } else {
    pen = bestPen.elem(arma::find(bestFloss == bestFloss(ID).min()));
  }
  
  Rcpp::List res = PLSinitilization(Y, X, k1, k2);
	arma::mat lambda_y = res["lambdaY"];
	arma::mat lambda_x = res["lambdaX"];
	arma::mat lambda_u = res["lambdaU"];
	
	Rcpp::List final = PLSfactorEMpenalty(Y, X, k1, k2, iter, lambda_y, lambda_x, lambda_u, pen(0), tol);
	return final;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMpenaltyGivenPen(arma::mat Y, arma::mat X, int k1, int k2, int iter, double pen, double tol){
  
  Rcpp::List res = PLSinitilization(Y, X, k1, k2);
  arma::mat lambda_y = res["lambdaY"];
	arma::mat lambda_x = res["lambdaX"];
	arma::mat lambda_u = res["lambdaU"];
	
	Rcpp::List final = PLSfactorEMpenalty(Y, X, k1, k2, iter, lambda_y, lambda_x, lambda_u, pen, tol);
	return final;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List EM_for_one_chunk(arma::mat W, arma::mat lambda, arma::mat psi, int k1, int k2, int p, 
                            arma::mat E_Z_w, arma::mat Var_Z_w){
  
  int n = W.n_rows;
  int q = W.n_cols - p;
  int k = k1 + k2;
  arma::mat E_ZZ_w_all = arma::zeros<arma::mat>(k, k);
  arma::mat y_E_z_s_all = arma::zeros<arma::mat>(p, k1);
  arma::mat x_E_z_s_all = arma::zeros<arma::mat>(q, k1);
  arma::mat u_E_z_x_all = arma::zeros<arma::mat>(q, k2);
  
  arma::mat lambda_u = lambda.submat(k1, p, k-1, p+q-1);
  arma::mat lambda_zero = arma::zeros<arma::mat>(k2, p);
  
  for (int j = 0; j < n; j++) {
    E_ZZ_w_all = E_ZZ_w_all + Var_Z_w + E_Z_w.row(j).t()*E_Z_w.row(j);
    y_E_z_s_all = y_E_z_s_all + W.row(j).subvec(0, p-1).t()*E_Z_w.row(j).subvec(0, k1-1);
    x_E_z_s_all = x_E_z_s_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(0, k1-1);
    u_E_z_x_all = u_E_z_x_all + W.row(j).subvec(p, p+q-1).t()*E_Z_w.row(j).subvec(k1, k-1);
  }
  arma::mat E_Z_s_Z_s_all = E_ZZ_w_all.submat(0, 0, k1-1, k1-1);
  arma::mat E_Z_x_Z_s_all = E_ZZ_w_all.submat(k1, 0, k-1, k1-1);
  arma::mat E_Z_s_Z_x_all = E_ZZ_w_all.submat(0, k1, k1-1, k-1);
  arma::mat E_Z_x_Z_x_all = E_ZZ_w_all.submat(k1, k1, k-1, k-1);
  
  arma::mat lambda_y_new = y_E_z_s_all*inv(E_Z_s_Z_s_all);
  arma::mat lambda_x_new = (x_E_z_s_all - lambda_u.t()*E_Z_x_Z_s_all)*inv(E_Z_s_Z_s_all);
  arma::mat lambda_u_new = (u_E_z_x_all - lambda_x_new*E_Z_s_Z_x_all)*inv(E_Z_x_Z_x_all);
  arma::mat lambda_1_new = join_rows(lambda_y_new.t(), lambda_x_new.t());
  arma::mat lambda_2_new = join_rows(lambda_zero, lambda_u_new.t());
  arma::mat lambda_new = join_cols(lambda_1_new, lambda_2_new);
  
  arma::vec psi_vec = psi.diag();
  arma::mat psi_inv = psi;
  psi_inv.diag() = 1/(psi_vec+0.000001);
  
  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.diag().ones();
  
  arma::mat shared_term = fastInverse(psi, lambda_new, TRUE);
  arma::mat term_4 = lambda_new*shared_term;
  E_Z_w = (term_4*W.t()).t();
  Var_Z_w = I - term_4*lambda_new.t();
  arma::mat psi_term = arma::zeros<arma::mat>(p+q, p+q);
  for (int j = 0; j < n; j++) {
    arma::mat term_5 = W.row(j).t() - lambda_new.t()*E_Z_w.row(j).t();
    psi_term =	psi_term + term_5*W.row(j);
  }
  psi.diag() = psi_term.diag()/n;
  
  return Rcpp::List::create(Rcpp::Named("E_Z") = E_Z_w,
                            Rcpp::Named("Var_z") = Var_Z_w,
                            Rcpp::Named("lambda") = lambda_new,
                            Rcpp::Named("psi") = psi);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PLSfactorEMchunk(arma::mat X, arma::mat Y, int k1, int k2, int iter, double chunk) {
  
  int p = Y.n_cols;
  int q = X.n_cols;
  int n = X.n_rows;
  int k = k1 + k2;
  arma::mat W = join_rows(Y, X); 
  arma::mat psi = arma::zeros<arma::mat>(p+q, p+q);
  psi.diag().ones();
  
  // initialization 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  bool SVDsuccess = false;
  while(SVDsuccess == false) {
    SVDsuccess =  svd(U, s, V, Y);
    if(SVDsuccess == false){
      Y.diag() = Y.diag() + 0.1;
    }
  }

  arma::mat matd = arma::zeros<arma::mat>(k1, k1);
  matd.diag() = sqrt(s(arma::span(0, k1-1)));
  arma::mat Z_Y = U.cols(0, k1-1)*matd;
  arma::mat Z_Y_norm = (Z_Y - mean(vectorise(Z_Y)))/stddev(vectorise(Z_Y));
  arma::mat lambda_y = initialize_lambda(Z_Y_norm, Y);
  arma::mat lambda_x = initialize_lambda(Z_Y_norm, X); 
  arma::mat res_x = X - Z_Y_norm*lambda_x;
  
  arma::mat U2;
  arma::vec s2;
  arma::mat V2;
  
  bool SVDsuccess2 = false;
  while(SVDsuccess2 == false) {
    SVDsuccess2 = svd(U2, s2, V2, res_x);
    if(SVDsuccess2 == false){
      res_x.diag() = res_x.diag() + 0.1;
    }
  }
   
  arma::mat matd2 = arma::zeros<arma::mat>(k2, k2);
  matd2.diag() = sqrt(s2(arma::span(0, k2-1)));
  arma::mat Z_U = U2.cols(0, k2-1)*matd2;
  arma::mat Z_U_norm = (Z_U - mean(vectorise(Z_U)))/stddev(vectorise(Z_U));
  arma::mat lambda_u = initialize_lambda(Z_U_norm, res_x); 
  
  arma::mat lambda_zero = arma::zeros<arma::mat>(k2, p);
  arma::mat lambda_1 = join_rows(lambda_y, lambda_x);
  arma::mat lambda_2 = join_rows(lambda_zero, lambda_u);
  arma::mat lambda = join_cols(lambda_1, lambda_2);

  arma::mat I = arma::zeros<arma::mat>(k, k);
  I.diag().ones();
  
  arma::mat shared_term = arma::zeros<arma::mat>(p+q, p+q);
  shared_term = fastInverse(psi, lambda, TRUE);
  arma::mat term_1 = lambda*shared_term;
  arma::mat E_Z_w = (term_1*W.t()).t();
  arma::mat Var_Z_w = I - term_1*lambda.t();
  
  arma::vec seqx = arma::zeros<arma::vec>(q);
  for (int i = 0; i < q; ++i) {
    seqx(i) = i;
  }
  int onepiece = ceil(q/chunk);
  arma::vec random_sample = RcppArmadillo::sample(seqx, q, FALSE); 
  
  arma::cube W_chunk_list = arma::zeros<arma::cube>(n, onepiece + p, chunk);
  arma::cube psi_chunk_list = arma::zeros<arma::cube>(onepiece + p, onepiece + p, chunk);
  arma::cube lambda_chunk_list = arma::zeros<arma::cube>(k, onepiece + p, chunk);
  
  int len = 0;
  for (int f = 0; f < chunk; f++) {
    
    if(f != chunk - 1) {
      arma::uvec IDs = arma::zeros<arma::uvec>(onepiece);
      for (int m = 0; m < onepiece; m++) {
        IDs(m) = random_sample(m + onepiece*f);
      }
      arma::mat subX = X.cols(IDs);
      W_chunk_list.slice(f) = join_rows(Y, subX);
      lambda_chunk_list.slice(f) = join_rows(lambda.cols(arma::span(0, p-1)), lambda.cols(IDs+p-1));
    } else {
      arma::uvec IDs = arma::zeros<arma::uvec>(q - (onepiece)*f);
      for (int m = 0; m < q - (onepiece)*f; m++) {
        IDs(m) = random_sample(m + onepiece*f);
      }
      len = IDs.n_elem;
      arma::mat subX = X.cols(IDs);
      W_chunk_list.slice(f).cols(arma::span(0, p + len-1)) = join_rows(Y, subX);
      lambda_chunk_list.slice(f).cols(arma::span(0, p + len-1)) = join_rows(lambda.cols(arma::span(0, p-1)), lambda.cols(IDs+p-1));
    }
    psi_chunk_list.slice(f).diag().ones();
  }
  
  for (int i = 0; i < iter; i++) {
    arma::mat E_Z_w_k = arma::zeros<arma::mat>(n, k);
    arma::mat Var_Z_w_k = arma::zeros<arma::mat>(k, k);
    
    for (int f = 0; f < chunk; f++) {
      if(f != chunk - 1) {
        Rcpp::List res = EM_for_one_chunk(W_chunk_list.slice(f), lambda_chunk_list.slice(f), psi_chunk_list.slice(f), k1, k2, p, E_Z_w, Var_Z_w);
        arma::mat mat1 = res(0);
        E_Z_w_k = E_Z_w_k + mat1;
        arma::mat mat2 = res(1);
        Var_Z_w_k = Var_Z_w_k + mat2;
        arma::mat mat3 = res(2);
        lambda_chunk_list.slice(f) = mat3;
        arma::mat mat4 = res(3);
        psi_chunk_list.slice(f) = mat4;
      } else {
        Rcpp::List res = EM_for_one_chunk(W_chunk_list.slice(f).cols(arma::span(0, p + len-1)), lambda_chunk_list.slice(f).cols(arma::span(0, p + len-1)), 
                                          psi_chunk_list.slice(f).submat(0,  0, p + len-1, p + len-1), k1, k2, p, E_Z_w, Var_Z_w);
        arma::mat mat1 = res(0);
        E_Z_w_k = E_Z_w_k + mat1;
        arma::mat mat2 = res(1);
        Var_Z_w_k = Var_Z_w_k + mat2;
        arma::mat mat3 = res(2);
        lambda_chunk_list.slice(f).cols(arma::span(0, p + len-1)) = mat3;
        arma::mat mat4 = res(3);
        psi_chunk_list.slice(f).submat(0, 0, p + len-1, p + len-1) = mat4;
      }  
    }
    E_Z_w = E_Z_w_k/chunk;
    Var_Z_w = Var_Z_w_k/chunk; 
    if(i == iter - 1) {
      lambda.cols(arma::span(0, p-1)) = lambda_chunk_list.slice(0).cols(0, p-1); 
      psi.submat(0, 0, p-1, p-1) = psi_chunk_list.slice(0).submat(0, 0, p-1, p-1); 
      for (int f = 0; f < chunk; f++) {
        if(f != chunk - 1) {
          arma::uvec IDs = arma::zeros<arma::uvec>(onepiece);
          for (int m = 0; m < onepiece; m++) {
            IDs(m) = random_sample(m + onepiece*f);
          }
          lambda.cols(IDs+p-1) = lambda_chunk_list.slice(f).cols(arma::span(p, p + onepiece - 1));
          psi.submat(IDs+p-1, IDs+p-1) = psi_chunk_list.slice(f).submat(p, p, p + onepiece - 1, p + onepiece - 1);
        } else {
          arma::uvec IDs = arma::zeros<arma::uvec>(q - (onepiece)*f);
          for (int m = 0; m < q - (onepiece)*f; m++) {
            IDs(m) = random_sample(m + onepiece*f);
          }
          lambda.cols(IDs+p-1) = lambda_chunk_list.slice(f).cols(arma::span(p, p + len - 1));
          psi.submat(IDs+p-1, IDs+p-1) = psi_chunk_list.slice(f).submat(p, p, p + len - 1, p + len - 1);
        }
      }
      Rcpp::List final = EM_for_one_chunk(W, lambda, psi, k1, k2, p, E_Z_w, Var_Z_w);
      arma::mat final_lambda = final(2);
      arma::mat final_psi = final(3);
      arma::mat final_E_Z_w = final(0);
      lambda = final_lambda;
      psi = final_psi;
      E_Z_w = final_E_Z_w;
    }  
  }  
  arma::vec psi_vec = psi.diag();
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda, Rcpp::Named("Z") = E_Z_w, Rcpp::Named("psi") = psi_vec);
}

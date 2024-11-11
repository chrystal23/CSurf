# include <RcppArmadillo.h>
# include <cmath>
# include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// ## single change with lasso from Lee et.al. The lasso for high dimensional regression with a
// # possible change point. JRSSB, 2016.
// [[Rcpp::export]]
Rcpp::List tlasso1(const arma::mat & x, const arma::vec & y, const arma::vec & q, 
                   double s, const arma::vec & grid_tau) {
  
  
  // Obtain environment containing function
  Rcpp::Environment package_env("package:ILAMM"); 
  // Make function callable from C++
  Rcpp::Function ncvxReg = package_env["ncvxReg"];  
  
  
  
  int ngrid = grid_tau.size();
  int nobs = y.size();
  int nreg = x.n_cols * 2.0;
  
  
  arma::mat delta_hat_grid(ngrid, nreg);
  arma::mat obj_v(ngrid, 1.0);
  arma::mat norm_x(1.0, nreg);
  arma::mat delta_hat(1.0, nreg);
  arma::vec ssr(ngrid);
  for(int i = 0; i < ngrid; i++){
    
    arma::mat x_reg_ind(nobs, x.n_cols);
    for (int k = 0; k < x.n_cols; k++){
      x_reg_ind.col(k) = x.col(k) % (q > grid_tau(i));
    }
    arma::mat x_reg = join_rows(x, x_reg_ind);
    
    Rcpp::List sm_fit = ncvxReg(x_reg, y, Rcpp::_["lambda"] = s);
    arma::vec beta = Rcpp::as<arma::vec>(sm_fit["beta"]);
    arma::vec yhat = x_reg * beta.subvec(1, x_reg.n_cols); 
    arma::vec uhat = y - yhat;
    ssr(i) = arma::as_scalar(uhat.t() * uhat) / nobs;
    for(int j = 0; j < nreg; j++){
      norm_x.col(j) = arma::sqrt( x_reg.col(j).t() * x_reg.col(j) / (nobs) );
    }
    
    arma::vec pen = norm_x * arma::abs(delta_hat_grid.row(i).t());
    obj_v.row(i) = ssr(i) + s * pen;
  }
  
  int opt = obj_v.index_min();
  delta_hat = delta_hat_grid.row(opt);
  double tau_hat = grid_tau(opt);
  ssr = ssr(opt);
  
  
  return Rcpp::List::create(Rcpp::Named("s") = s, Rcpp::Named("tau.hat") = tau_hat,
                            Rcpp::Named("delta.hat") = delta_hat, Rcpp::Named("ssr") = ssr);
  
}



arma::mat colCen(const arma::mat & X){
  int nRows = X.n_rows;
  int nCols = X.n_cols;
  arma::mat out(nRows, nCols);
  for(int i = 0; i < nCols; i++){
    out.col(i) = X.col(i) - mean(X.col(i));
  }
  return(out);
}


// [[Rcpp::export]]
arma::mat ConstructBasis(const arma::vec& x, const int m = 3){
  
  int n = x.size();
  int nk = n - m; // number of knots
  double x_shift = arma::mean(x); // shift used to enhance stability
  
  
  arma::vec xord = arma::sort(x);
  arma::vec knots = arma::zeros(n-m);
  bool odd = (m % 2 == 0);
  if(odd) knots = xord.subvec((m/2), (n - m/2 - 1)) - x_shift;
  else knots = xord.subvec(((m-1)/2), (n - (m-1)/2 - 1)) - x_shift;
  arma::vec x_c = x - x_shift;
  
  arma::mat X(n, n-1, arma::fill::zeros);
  for (int i = 0; i < (m-1); i++) {
    X.col(i) = pow(x_c, i+1);
  }
  for (int i = 0; i < nk; i++) {
    X.col(i+m-1) = pow(x_c - knots(i), m - 1) % (x_c > knots(i));
  }
  
  arma::mat X_cen = colCen(X);
  
  return X_cen;
  
}

arma::vec cmptLambda(const arma::vec& beta, const double lambda, 
                     const int p, const arma::vec& group, const int G,
                     const std::string penalty) {
  
  arma::vec subNorm = arma::zeros(G);
  for (int i = 1; i <= p; i++) {
    subNorm(group(i)) += beta(i) * beta(i);
  }
  subNorm = arma::sqrt(subNorm);
  
  arma::vec rst = arma::zeros(G);
  if (penalty == "gLasso") {
    rst = lambda * arma::ones(G);
  } else if (penalty == "gSCAD") {
    double a = 3.7;
    double abBeta;
    for (int i = 0; i < G; i++) {
      abBeta = subNorm(i);
      if (abBeta <= lambda) {
        rst(i) = lambda;
      } else if (abBeta <= a * lambda) {
        rst(i) = (a * lambda - abBeta) / (a - 1);
      }
    }
  } else if (penalty == "gMCP") {
    double a = 3;
    double abBeta;
    for (int i = 0; i < G; i++) {
      abBeta = subNorm(i);
      if (abBeta <= a * lambda) {
        rst(i) = lambda - abBeta / a;
      }
    }
  }
  return rst;
}



double loss(const arma::vec& Y, const arma::vec& Ynew) {
  double rst = 0;
  rst = arma::mean(arma::square(Y - Ynew)) / 2;
  
  return rst;
}

arma::vec gradLoss(const arma::mat& X, const arma::vec& Y, const arma::vec& beta, const bool interecept) {
  arma::vec res = Y - X * beta;
  arma::vec rst = arma::zeros(beta.size());
  rst = -1 * (res.t() * X).t();
  
  if (!interecept) {
    rst(0) = 0;
  }
  return rst / Y.size();
}

arma::vec updateBeta(const arma::mat& X, const arma::vec& Y, arma::vec beta, 
                     const int p, const arma::vec& group, const int G, const double phi,
                     const arma::vec& Lambda, const bool intercept) {
  
  
  arma::vec subNorm = arma::zeros(G);
  arma::vec betaNew = beta - gradLoss(X, Y, beta, intercept) / phi;
  for (int i = 1; i <= p; i++) {
    subNorm(group(i)) += betaNew(i) * betaNew(i);
  }
  subNorm = arma::max(1.0 - Lambda / (phi * arma::sqrt(subNorm)), arma::zeros(G));
  for (int i = 1; i <= p; i++) {
    betaNew(i) *= subNorm(group(i));
  }
  
  
  return betaNew;
}

double cmptPsi(const arma::mat& X, const arma::vec& Y, const arma::vec& betaNew,
               const arma::vec& beta, const double phi, const bool intercept) {
  arma::vec diff = betaNew - beta;
  double rst = loss(Y, X * beta)
    + arma::as_scalar((gradLoss(X, Y, beta, intercept)).t() * diff)
    + phi * arma::as_scalar(diff.t() * diff) / 2;
    return rst;
}


Rcpp::List LAMM(const arma::mat& X, const arma::vec& Y, const arma::vec& Lambda, 
                arma::vec beta, const int p, const arma::vec& group, const int G,
                const double phi, const double gamma, const bool interecept) {
  double phiNew = phi;
  arma::vec betaNew = arma::vec();
  double FVal;
  double PsiVal;
  while (true) {
    betaNew = updateBeta(X, Y, beta, p, group, G, phiNew, Lambda, interecept);
    FVal = loss(Y, X * betaNew);
    PsiVal = cmptPsi(X, Y, betaNew, beta, phiNew, interecept);
    if (FVal <= PsiVal) {
      break;
    }
    phiNew *= gamma;
  }
  return Rcpp::List::create(Rcpp::Named("beta") = betaNew, Rcpp::Named("phi") = phiNew);
}


// [[Rcpp::export]]
Rcpp::List ncvxGroupReg(arma::mat X, const arma::vec& Y, const arma::vec& group, const int G,
                        double lambda = -1, std::string penalty = "gSCAD", const double phi0 = 0.001,
                        const double gamma = 1.5, const double epsilon_c = 0.001,
                        const double epsilon_t = 0.001, const int iteMax = 500,
                        const bool intercept = false, const bool itcpIncluded = false) {
  if (!itcpIncluded) {
    arma::mat XX(X.n_rows, X.n_cols + 1);
    XX.cols(1, X.n_cols) = X;
    XX.col(0) = arma::ones(X.n_rows);
    X = XX;
  }
  int n = Y.size();
  int p = X.n_cols - 1;
  if (lambda <= 0) {
    double lambdaMax = arma::max(arma::abs(Y.t() * X)) / n;
    double lambdaMin = 0.01 * lambdaMax;
    lambda = std::exp((long double)(0.7 * std::log((long double)lambdaMax)
                                      + 0.3 * std::log((long double)lambdaMin)));
  }
  arma::vec beta = arma::zeros(p + 1);
  arma::vec betaNew = arma::zeros(p + 1);
  // Contraction
  arma::vec Lambda = cmptLambda(beta, lambda, p, group, G, penalty);
  double phi = phi0;
  int ite = 0;
  Rcpp::List listLAMM;
  while (ite <= iteMax) {
    ite++;
    listLAMM = LAMM(X, Y, Lambda, beta, p, group, G, phi, gamma, intercept);
    betaNew = Rcpp::as<arma::vec>(listLAMM["beta"]);
    phi = listLAMM["phi"];
    phi = std::max(phi0, phi / gamma);
    if (arma::norm(betaNew - beta, "inf") <= epsilon_c) {
      break;
    }
    beta = betaNew;
  }
  int iteT = 0;
  // Tightening
  if (penalty != "gLasso") {
    arma::vec beta0 = arma::zeros(p + 1);
    while (iteT <= iteMax) {
      iteT++;
      beta = betaNew;
      beta0 = betaNew;
      Lambda = cmptLambda(beta, lambda, p, group, G, penalty);
      phi = phi0;
      ite = 0;
      while (ite <= iteMax) {
        ite++;
        listLAMM  = LAMM(X, Y, Lambda, beta, p, group, G, phi, gamma, intercept);
        betaNew = Rcpp::as<arma::vec>(listLAMM["beta"]);
        phi = listLAMM["phi"];
        phi = std::max(phi0, phi / gamma);
        if (arma::norm(betaNew - beta, "inf") <= epsilon_t) {
          break;
        }
        beta = betaNew;
      }
      if (arma::norm(betaNew - beta0, "inf") <= epsilon_t) {
        break;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("beta") = betaNew, Rcpp::Named("phi") = phi,
                            Rcpp::Named("penalty") = penalty, Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("IteTightening") = iteT);
}



#' @title CSurf
#'
#' @description Change Surface (CSurf) regression method (additive model based) for nonlinear subgroup identification
#'
#'
#' @param y A vector (length n) of the response for subgroup identification.
#' @param X A matrix (n x p) of the covariates for subgroup identification. The covariates are also used as potential thresholding variables for subgroup identification.
#' @param f_hat0 A vector (length n) of the initialization for the threshold function (change surface).
#' @param Tighten A boolean indicating whether to include the tightening stage in optimization.
#' @param tol_out The tolerance value for optimization (stopping criterion).
#' @param maxiter_out The maximal number of iterations allowed for optimization.
#' @param m An integer indicating the power of the threshold function.
#' @param lam2rate The ratio between the two tunning parameters, lambda1 and lambda2.
#'
#' @return A list containing the following fields: 
#' \item{f_all.hat}{The estimated threshold functions. Each column represent one component function.}
#' \item{tau.hat}{The estimated change points over the threshold functions.}
#' \item{alp.hat}{The estimated regression coefficients.}
#' \item{coverged_out}{An boolean indicating whether the algorithm successfully converged.}
#' \item{iter_out}{The number of iterations it took to converge.}
#' \item{...}{Other estimation results.}
#' 
#' @examples
#' data(syn_data)
#' results <- CSurf(y = y, X = X, f_hat0 = f_hat0, h = .5, maxiter_out = 10, tol_out = 0.0001)

CSurf <- function(y, X, f_hat0, h, Tighten = FALSE, tol_out = 1e-4, maxiter_out = 10, m = 3, lam2rate = 1){
  
  n = nrow(X)
  p = ncol(X)
  ##
  tau.trace = f.trace = alp.trace = X.trace = vector("list", maxiter_out)
  SIC = loss = IteContraction = rep(NA, maxiter_out)
  iter_out = 1
  converged_out <- FALSE
  while (iter_out <= maxiter_out){
    
    
    ####### step 1: update number of change 
    cat('-----------------------------------------', '\n')
    cat("Step 1: Start Splitting", '\n')
    cat('-----------------------------------------', '\n')
    ## order data
    ord = order(f_hat0)
    X = X[ord,]
    y = y[ord]
    f_hat0 = f_hat0[ord]
    
    ## choose the best split c
    c <- c(round(n*0.02*c(0.25, 0.5, 0.75, 1)))
    split.fit <- NULL
    for(j in 1:length(c)) {
      split.fit[[j]] <- split_Fun(X = X, Y = y, c[j], penalty = "gSCAD")
    }
    BIC <- sapply(split.fit, function(x){x$BIC})
    split.fit <- split.fit[[which.min(BIC)]]
    mcp.cp1 <- split.fit$mcp.cp1
    Numtau <- length(mcp.cp1)
    cat("choose no.", which.min(BIC), '\n')
    
    ## choose c = 10
    # split.fit = split_Fun(X = X, Y = y, c = 10, penalty = "gSCAD")
    # mcp.cp1 = split.fit$mcp.cp1
    # Numtau = length(mcp.cp1)
    # cat("Possiable changes in groups: ", mcp.cp1, " selected by group SCAD using LAMM.", '\n')
    
    
    if(Numtau > 0) {
      
      #### step 2: update tau and alp
      cat('-----------------------------------------', '\n')
      cat("Step 2: Estimate thresholds", '\n')
      cat('-----------------------------------------', '\n')
      ## get tau
      Itil = split.fit$Itil
      tau.hat = rep(Inf, Numtau)
      for (k in 1:Numtau) {
        x = X[Itil[[k]],]
        yy = y[Itil[[k]]]
        grid.tau = f_hat0[Itil[[k]]]
        ## low dimension
        # id = css(y = yy, X = x)
        # tau.hat[k] = ifelse(id > 0, grid.tau[id], Inf)
        
        ## high dimension: projection
        # tau.hat[k] = tau_hd_Fun(x = x, yy = yy, grid.tau = grid.tau)
        ## high dimension: threshold lasso
        tau.hat[k] = tau_Fun(x = x, yy = yy, z = grid.tau)
      }
      tau.hat = tau.hat[tau.hat< Inf]
      ## get alp
      fit_alp = alp_Fun(tau = tau.hat, X = X, y = y, f = f_hat0, penalty = "SCAD")
      alp.hat = fit_alp$alp
      tau.hat = fit_alp$tau
      tau.hat = tau.hat[tau.hat< Inf]
      cat("Thresholds at ", tau.hat, '\n')
      # cat("Coeff: ", alp.hat, '\n')
      if (length(tau.hat) == 0) {
        cat("NO Change.", "\n")
        tau.hat = NULL
        fit = ILAMM::cvNcvxReg(X, y, penalty = "SCAD")
        alp.hat = as.vector(fit$beta[-1,]) 
        loss[iter_out] = mean((y - X%*%alp.hat)^2)
        
        return(list(alp.hat = alp.hat, tau.hat = tau.hat, loss = loss))
        
        break
      }
      
      ### step 3: update f
      cat('-----------------------------------------', '\n')
      cat("Step 3: Start smoothed LAMM", '\n')
      cat('-----------------------------------------', '\n')
      
      ## choose the lambda1
      
      m = m  # 2,4
      if (iter_out <= 1) {
        f_all.hat0 = NULL
      } else{
        f_all.hat0 = f_all.hat[ord,]
      }
      print(utils::head(f_all.hat0))
      
      
      lam_max <- grad_Fun(f_hat = matrix(0, nrow = n, ncol = p), y=y, X=X, alp=alp.hat, tau=tau.hat, h=h)
      lam_max <- norm_f(lam_max)
      
      c <- seq(lam_max*0.02, lam_max*0.8, length.out=10)
      # c <- c(7,8,9,10)
      # c <- c(0.1, .5, 1, 1.5, 2)
      
      BIC_cm = f_all.hat = Niter = NULL
      for(j in 1:length(c)) {
        
        lambda1 = c[j]#*sqrt(log(p)/(n*h)) 
        lambda2 = lambda1^2 * lam2rate  # for minimax rate
        cat("lambda1 = ", lambda1, '\n')
        fit <- LAMMad(f_all.hat0 = f_all.hat0, phi0 = .01, y = y, X = X, Z = X, m = m, 
                      lambda1 = lambda1, lambda2 = lambda2, max_iter = 20, tol = tol_out, 
                      Tighten = Tighten, alp = alp.hat, tau = tau.hat, h = h)
        ##
        f_all.hat[[j]] = fit$fhat
        Niter = c(Niter, fit$IteContraction)
        if (sum(abs(fit$fhat))==0) BIC_cm = c(BIC_cm, Inf)
        else BIC_cm = c(BIC_cm, 
                        n * log(obj_f_Fun(f_hat = fit$fhat, y = y, X = X, 
                                          alp = alp.hat, tau = tau.hat, h = h)) 
                        + sum(colMeans(abs(fit$fhat)) != 0)*2 ) #(log(n)*log(p)/2) #log(n)
        
      }
      f_all.hat <- f_all.hat[[which.min(BIC_cm)]]
      IteContraction[iter_out] <- Niter[which.min(BIC_cm)]
      cat("choose no.", which.min(BIC_cm), '\n')
      loss[iter_out] = obj_f_Fun(f_hat = f_all.hat, y = y, X = X, alp = alp.hat, tau = tau.hat, h = h) * 2

      f_hat0 = rowSums(f_all.hat)
      
    } else {
      
      cat("NO Change.", "\n")
      tau.hat = NULL
      fit = ILAMM::cvNcvxReg(X, y, penalty = "SCAD")
      alp.hat = as.vector(fit$beta[-1,]) 
      loss[iter_out] = mean((y - X%*%alp.hat)^2)
      
      return(list(alp.hat = alp.hat, tau.hat = tau.hat, loss = loss))
      
      break
    }
    
    
    ##
    tau.trace[[iter_out]] <- tau.hat
    f.trace[[iter_out]] <- f_all.hat
    alp.trace[[iter_out]] <- alp.hat
    X.trace[[iter_out]] <- X
    SIC[iter_out] <- n * log(loss[iter_out])+ length(tau.trace[[iter_out]]) * log(n)
    
    
    ### epsilon-optimal stopping criterion
    if ( iter_out > 1) {

      temp_res1 <- mean((f_all.hat - f_all.hat0)^2)
      temp_res2 <- mean((f_all.hat0)^2)
      print(temp_res1/(temp_res2+1e-30))
      if(temp_res1/(temp_res2+1e-30) <= tol_out) {
        cat("Finish!", "\n")
        converged_out <- TRUE
        break
      }

    }
    
    cat("iteration = ", iter_out <- iter_out + 1, "\n")
    
    
  }
  
  
  res = list(f_all.hat = f_all.hat, 
             X.trace = X.trace, yord = y, Xord = X,
             alp.hat = alp.hat, 
             tau.hat = tau.hat, 
             tau.trace = tau.trace,
             alp.trace = alp.trace,
             f_all.trace = f.trace, 
             loss = loss, SIC = SIC, 
             converged_out = converged_out,
             iter_out = iter_out, IteContraction = IteContraction[1:iter_out])
  
  return(res)
  
}



###############################################
## step 1: split
split_Fun <- function(X, Y, c, trace = TRUE, penalty=c("gLasso", "gSCAD")) {
  
  
  penalty <- match.arg(penalty)
  
  n = length(Y)
  p = dim(X)[2]
  m = ceiling(sqrt(n * c))
  q = floor(n / m)
  
  ###########transform###################
  K_temp=matrix(0, nrow = q, ncol=q, byrow=TRUE)
  for(i in 1:q) K_temp[i,1:i]=rep(1,i)
  
  X_temp= X #cbind(1,X)
  Y_temp=c(Y)
  
  x=NULL
  y=NULL
  x[[1]]=as.matrix(X_temp[1:((n-(q-1)*m)),])
  y[[1]]=Y_temp[1:((n-(q-1)*m))]
  
  for(i in 2:q) {
    x[[i]] = as.matrix(X_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m)),])
    y[[i]] = Y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
  }
  X_temp1 <- lapply(1:length(x), function(j, mat, list)
    kronecker(K_temp[j, , drop = FALSE], x[[j]]), mat = K_temp, list = x)
  Xn=do.call("rbind",X_temp1)
  rm(X_temp1)
  #
  Yn=NULL
  for(i in 1:q) {
    Yn=c(Yn,y[[i]])
  }
  #
  groupindex <- kronecker(c(1:q),rep(1,p))
  G = length(unique(groupindex))
  group = c(0, groupindex - 1)
  ## group SCAD by I-LAMM
  nLambda = 10
  lambda = lambdaseq(X = Xn, y = Yn, groupindex = groupindex, lambdaRatio = 0.001, nLambda = nLambda)$lambda
  beta.all = matrix(0, ncol(Xn), nLambda)
  loss = df = rep(NA, nLambda)
  for (i in 1:nLambda) {
    fit = ncvxGroupReg(X = Xn, Y = Yn, group = group, G = G, lambda = lambda[i], penalty = penalty, 
                       intercept = FALSE, itcpIncluded = FALSE)
    beta.all[,i] = fit$beta[-1]
    loss[i] = mean((Yn - Xn%*%beta.all[,i])^2)/2
    df[i] = sum(beta.all[,i]!=0)/p
  }
  BIC = n*log(loss) + df*log(n)*log(p*q)
  mcp.coef = beta.all[,which.min(BIC)]
  BIC = BIC[which.min(BIC)]
  loss = loss[which.min(BIC)]
  
  ##
  mcp.coef.v.m <- abs(matrix(c(mcp.coef), q, p, byrow = TRUE))
  mcp.coef.m <- c(apply(mcp.coef.v.m, 1, norm))
  mcp.cp <- which(mcp.coef.m!=0)
  if (length(mcp.cp) > 1) {
    for (i in 2:length(mcp.cp))
    {
      if (mcp.cp[i] - mcp.cp[i - 1] == 1)
        mcp.cp[i] = 0
    }
  }
  mcp.cp1 <- mcp.cp[mcp.cp > 1 & mcp.cp < q]
  
  ## consistent change location interval
  Khat <- length(mcp.cp1)
  if(Khat > 0){
    Itil = NULL
    for(k in 1:Khat) {
      Itil[[k]] = c( (n-(q - mcp.cp1[k]+2)*m) : (n-(q - mcp.cp1[k])*m) )
    }
    
  } else {
    
    Itil = 1:n
  }
  
  return(list(BIC = BIC, loss = loss, mcp.cp1 = mcp.cp1, Itil = Itil))
  
}


## step 1:  single change point detection in high dimensional regression model
## given estimate tau for each interval \tilde I_j
tau_Fun <- function (x, yy, z) {
  
  nx = ncol(x)
  nobs <- length(yy)
  ##
  trim.factor <- round(length(z)*.10,0)
  grid.tau <- z[trim.factor:(length(z)-trim.factor)]
  t_0 <- grid.tau[length(grid.tau)]
  
  lambda_max = max(abs(t(x)%*%yy/nobs))*stats::sd(yy)
  nLambda = 10
  lambda_all=exp(seq(log(lambda_max),log(lambda_max*0.001),length=nLambda))
  
  ##
  # y.odd= yy[seq(1, (nobs-1), 2)]
  # X.odd= x[seq(1, (nobs-1), 2),]
  # z.odd = z[seq(1,(nobs-1),2)]
  y.even =  yy[seq(2, nobs, 2)]
  X.even= x[seq(2, nobs, 2),]
  z.even =  z[seq(2, nobs, 2)]
  ##
  delta.hat = matrix(NA, nrow = 2*nx, ncol = nLambda)
  tau.hat = res.test = HBIC = NULL
  for (j in 1:nLambda) {
    # fit = tlasso1(x = x, y = yy, q = z, s = lambda_all[j], grid_tau = grid.tau)
    fit = tlasso(x = x, y = yy, q = z, s = lambda_all[j], grid.tau = grid.tau)
    delta.hat[,j] = fit$delta.hat
    # HBIC = c(HBIC, length(yy)*log(fit$ssr) 
    #         + sum(delta.hat[,j]!=0)*log(length(yy)) )
    
    x.reg <- cbind(X.even, X.even*(z.even > fit$tau.hat))
    res.test = c(res.test, mean((y.even - x.reg%*%delta.hat[,j])^2))
    tau.hat = c(tau.hat, fit$tau.hat)
  }
  tau.hat = tau.hat[which.min(res.test)]
  delta.hat = delta.hat[, which.min(res.test)]
  lambda.min = lambda_all[which.min(res.test)]
  
  cat("choose no.", which.min(res.test), " among ", lambda_all,  '\n')
  
  ## if delta = 0 than no change
  if(sum(abs(delta.hat[(nx+1):(2*nx)])) == 0) tau.hat = Inf
  ## boundary 
  if(tau.hat == t_0) tau.hat = Inf
  
  return(tau.hat = tau.hat)
  
}


## single change with lasso from Lee et.al. The lasso for high dimensional regression with a
# possible change point. JRSSB, 2016.
tlasso <- function (x, y, q, s, grid.tau) {
  
  ngrid <- length(grid.tau)
  nobs <- length(y)
  nreg <- ncol(x)*2
  delta.hat.grid <- matrix(rep(NA), nrow=ngrid, ncol=nreg)
  obj.v <- matrix(rep(NA), nrow=ngrid, ncol=1)
  norm.x <- matrix(rep(NA), nrow=1, ncol=nreg)
  delta.hat <- matrix(rep(NA), nrow=1, ncol=nreg)
  ssr <- rep(NA, ngrid)
  for(i in 1:ngrid) {
    ind <- ( q > grid.tau[i] )
    x.reg <- cbind(x,x*ind)
    for(j in 1:nreg){
      norm.x[1,j] <- sqrt( t((x.reg[,j]))%*% (x.reg[,j]) / (nobs) )
    }
    # m = wncvxReg(X = x.reg, Y = y, weight = c(0, Matrix::drop(norm.x)), lambda = s, penalty = "SCAD", intercept = FALSE)
    # delta.hat.grid[i,] = Matrix::drop(m$beta[-1])
    # yhat = x.reg%*%delta.hat.grid[i,]
    m = glmnet::glmnet(x = x.reg, y = y, lambda = s, standardize = FALSE, penalty.factor = Matrix::drop(norm.x), intercept = FALSE)
    delta.hat.grid[i,] <- Matrix::drop(m$beta) #coef(m, s=s)[-1]
    yhat <- stats::predict(m, newx = x.reg, s = s)
    uhat <- y - yhat
    ssr[i] <- t(uhat)%*%uhat / nobs
    p <- SCAD(x = norm.x %*% abs(delta.hat.grid[i,]), lambda = s)
    obj.v[i,1] <- ssr[i] + p # s*p    	
  }
  opt<-which.min(obj.v)
  delta.hat <- delta.hat.grid[opt,]
  tau.hat <- grid.tau[opt]
  ssr <- ssr[opt]
  
  
  # rval <- list(opt= opt, est=c(s, tau.hat, delta.hat))
  rval <- list(s = s, tau.hat = tau.hat, delta.hat = delta.hat, ssr = ssr)
  
  return(rval)
  
}



#An  quick hypothesis test of a singe structral break  
#by  the likelihood ratio method proposed by Davis et al.
#(1995, Testing for a change in
#the parameter values and order of an autoregressive model,
#Annals of Statistics 23(1), 282-304).
#lm(formula=y~X-1). Confidence level is 0.05.
#Return the estimate of chang-point if there is one and 0 otherwise.
css <- function(y, X) {
  #
  n <- dim(X)[1]
  p <- dim(X)[2]
  bn <- (2 * log(log(n)) + (p + 1)/2 * log(log(log(n))) - log(gamma((p +
                                                                       1)/2)))^2/(2 * log(log(n)))
  an <- sqrt(bn/(2 * log(log(n))))
  
  # Recursive Calculation of the likelihood ratio statistic
  
  S <- t(y[1:(n)]) %*% X[1:(n), ]
  C1 <- solve(t(X[1:(n), ]) %*% X[1:(n), ], t(S))
  Q <- S %*% C1
  sigma <- (sum(y^2) - Q)/n
  Dic <- c(rep(0, n))
  C0 <- solve(t(X[1:(2 * p), ]) %*% X[1:(2 * p), ])
  C1 <- solve(t(X[(2 * p + 1):n, ]) %*% X[(2 * p + 1):n, ])
  S0 <- t(y[1:(2 * p)]) %*% X[1:(2 * p), ]
  S1 <- t(y[(2 * p + 1):(n)]) %*% X[(2 * p + 1):(n), ]
  for (i in (2 * p + 1):(n - (2 * p + 1))) {
    xi <- matrix(X[i, ], 1, p)
    C0 <- C0 - C0 %*% t(xi) %*% xi %*% C0/c(1 + xi %*% C0 %*% t(xi))
    C1 <- C1 + C1 %*% t(xi) %*% xi %*% C1/c(1 - xi %*% C1 %*% t(xi))
    S0 <- S0 + y[i] * xi
    S1 <- S1 - y[i] * xi
    Q1 <- S0 %*% C0 %*% t(S0)
    Q2 <- S1 %*% C1 %*% t(S1)
    Dic[i] <- Q1 + Q2
  }
  
  cp <- which(Dic == max(Dic))
  d <- (Dic[cp] - Q)/sigma
  if (((d - bn)/an) >= 2 * log(-2/(log(1 - 0.05))))
    return(cp) else return(0)
  
}



## given gam and tau, estimate alpha = c(beta, delta) 
alp_Fun <- function(tau, X, y, f, penalty = "SCAD") {
  
  n = nrow(X)
  p = ncol(X)
  
  
  ts = c(-Inf, tau, Inf)
  Numtau = length(ts)-1
  bb = matrix(NA, p, Numtau)
  for (k in 1:Numtau) {
    id = which(f > ts[k] & f <= ts[k+1])
    X.reg <- X[id,]
    y.reg <- y[id]
    
    BIC = NULL
    # lambda_all = lambdaseq(X.reg, y.reg, groupindex = 1:ncol(X.reg))$lambda
    # nLambda = length(lambda_all)
    lambda_max = max(abs(t(X.reg)%*%y.reg/length(y.reg)))
    nLambda = 10
    lambda_all=exp(seq(log(lambda_max),log(lambda_max*0.0001),length=nLambda))
    alp = matrix(NA, nrow = ncol(X.reg), ncol = nLambda)
    for (i in 1:nLambda) {
      if (penalty == "SCAD"){
        fit = ILAMM::ncvxReg(X = X.reg, Y = y.reg, lambda = lambda_all[i], penalty = "SCAD",
                      intercept = FALSE)
        alp[,i] = as.vector(fit$beta[-1,])
      } else {
        fit = glmnet::glmnet(x = X.reg, y = y.reg, lambda = lambda_all[i], intercept = FALSE, standardize = FALSE)
        alp[,i] = as.vector(fit$beta)
      }
      
      
      BIC = c(BIC, length(y.reg)*log(mean((y.reg - X.reg%*%alp[,i])^2)) 
              + sum(alp[,i]!=0)*(log(length(y.reg))) )
    }
    bb[,k] = alp[, which.min(BIC)]
    
  }
  delta = matrixStats::rowDiffs(bb)
  idx = apply(matrixStats::rowDiffs(bb), 2, function(x) sum(abs(x))==0 )
  if(sum(idx) > 0) {
    tau[idx] = Inf
    delta = matrixStats::rowDiffs(bb)[,!idx]
  } 
  beta = bb[,1]
  alp = c(beta, delta)
  
  return(list(alp = alp, tau = tau, beta = beta, delta = delta))
  
}





###########################################################
########################## LAMM algorithm for additive model
LAMMad <- function(f_all.hat0 = NULL, y, X, Z, m = 3, phi0 = .1, nu = 2, lambda1, lambda2, max_iter = 100, tol = 1e-4, 
                 Tighten = TRUE, alp, tau, h) {
  
  # Initialize parameters using warm starts.
  n <- nrow(Z)
  d <- ncol(Z)
  if(is.null(f_all.hat0)) {
    f_all.hat0 <- matrix(0, ncol = d, nrow = n)
  }
  old_ans <- f_all.hat0
  ## Contraction
  cat("Contraction stage", "\n")
  Lambda1 = apply(old_ans, 2, SCAD_f.dev, lambda = lambda1)
  phi = phi0
  # Initialize some objects.
  iter <- 1
  converged <- FALSE
  
  while(iter < max_iter & !converged) {
    
    
    linfit = LineSearch(nu = nu, phi = phi, y = y, f_hat = old_ans, X = X, Z = Z, m = m, 
                        Lambda1 = Lambda1, lambda2 = lambda2, alp = alp, tau = tau, h = h)
    new_ans = linfit$new_ans
    phi = linfit$phi
    phi = max(phi0, phi / nu);
    print(phi)
    
    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    temp_res1 <- mean((new_ans - old_ans)^2) 
    temp_res2 <- mean((old_ans)^2) 
    if(temp_res1/(temp_res2+1e-30) <= tol) {
      converged <- TRUE
    } else {
      iter <- iter + 1;
      old_ans <- new_ans;
    }
  }
  if(converged == FALSE) {
    expr <- paste0("Algorithm did not converge for Lambda1 = ", signif(lambda1, 2),
                   " and Lambda2 = ", signif(lambda2,2));
    warning(expr)
  }
  
  
  ## Tightening
  iteT = 0;
  if (Tighten) {
    cat("Tightening stage", "\n")
    
    while (iteT <= max_iter) {
      iteT = iteT + 1
      old_ans = new_ans
      old_ans0 = new_ans
      Lambda1 = apply(old_ans, 2, SCAD_f.dev, lambda = lambda1)
      phi = phi0
      ite = 0;
      while (ite <= max_iter) {
        ite = ite + 1
        linfit = LineSearch(nu = nu, phi = phi, y = y, f_hat = old_ans, X = X, Z = Z, m = m, 
                            Lambda1 = Lambda1, lambda2 = lambda2, alp = alp, tau = tau, h = h)
        new_ans = linfit$new_ans
        phi = linfit$phi
        phi = max(phi0, phi / nu)
        print(phi)
        
        temp_res1 <- mean((new_ans - old_ans)^2) 
        temp_res2 <- mean((old_ans)^2) 
        if(temp_res1/(temp_res2+1e-30) <= tol){
          break
        }
        old_ans <- new_ans;
      }
      if (mean((new_ans - old_ans0)^2)/(mean((old_ans0)^2)+1e-30) <= tol) {
        break;
      }
    }
    
  }
  
  
  list("fhat" = new_ans,
       "conv" = converged,
       "IteContraction" = iter,
       "IteTightening" = iteT)
}



##
LineSearch <- function(nu = 2, phi, y, f_hat, X, Z, m, Lambda1, lambda2, tau, alp, h) {
  
  n = nrow(X)
  phiNew = phi
  old_ans = f_hat
  vector_r <- grad_Fun(f_hat = old_ans, y = y, X = X, alp = alp, tau = tau, h = h)
  F_old <- obj_f_Fun(f_hat = old_ans, y = y, X = X, alp = alp, tau = tau, h = h)
  
  while (TRUE) {
    
    new_ans <- update_f(f_hat = old_ans, X = X, Z = Z, phi = phiNew, m = m, y = y,
                        lambda1 = Lambda1, lambda2 = lambda2, alp = alp, tau = tau, h = h)
    FVal = obj_f_Fun(f_hat = new_ans, y = y, X = X, alp = alp, tau = tau, h = h)
    PsiVal = F_old + sum(vector_r %*% (new_ans - old_ans)/n) + phiNew * sum((new_ans - old_ans)^2)/(n*2) 
    if (FVal <= PsiVal) {
      break;
    }
    phiNew = nu*phiNew;
  }
  
  # Return the next iterate based on the selected step size.
  return(list(new_ans = new_ans, phi = phiNew))
}


##
solve_prox <- function(r, Zj, lbd1, lbd2, m = 3) {

  # m = 3
  n = length(r)
  psi = ConstructBasis(x = Zj, m = m)

  # library(glmnet)
  sm_fit = glmnet::glmnet(x = psi, y = r, lambda = lbd2, intercept = FALSE, standardize = FALSE,
                  penalty.factor = c(rep(0, m-1),rep(1,ncol(psi)-(m-1))) )
  fhat_j = glmnet::predict.glmnet(object = sm_fit, newx = psi, s = lbd2)  # lambda2 = sqrt(log(n-1)/n)

  ## soft-threshold
  fhat_j = max(1 - lbd1/norm_f(fhat_j),0)* fhat_j

  return(fhat_j)

}


##
update_f <- function(f_hat, X, Z, phi, m = 3, y, lambda1, lambda2, alp, tau, h, ind.to.be.positive = 1) {
  
  n <- nrow(f_hat)
  d <- ncol(f_hat)
  
  
  vector_r <- grad_Fun(f_hat = f_hat, y = y, X = X, alp = alp, tau = tau, h = h)
  inter_step <- apply(f_hat, 2, "-", vector_r/phi)
  ans <- matrix(0, nrow = n, ncol = d)
  
  for(j in 1:d) {
    ans[, j] <- solve_prox(r = inter_step[,j] - mean(inter_step[,j]), Zj = Z[, j],
                           lbd1 = lambda1[j]/phi, lbd2 = lambda2/phi, m = m)
  }
  
  # ## identifiable constrain
  fvar = sqrt(sum(apply(ans, 2, stats::var))) + 1e-30
  ans = scale(ans, center = TRUE, scale = rep(fvar,d))
  if(!is.null(ind.to.be.positive)){
    tmp = (Z[, ind.to.be.positive] - mean(Z[, ind.to.be.positive]))%*%ans[, ind.to.be.positive]
    if(tmp < 0) ans <- -1*ans      # for the (sign) identifiability
  }
  
  return(ans)
}
##


#' @title f0_Fun
#'
#' @description The function to generate threshold function initialization using change plane model.
#'
#'
#' @param y A vector (length n) of the true response.
#' @param X A matrix (n x p) of the covariates. The covariates are also used as potential thresholding variables for subgroup identification.
#'
#' @return A list containing the following fields: 
#' \item{f0}{A vector (length n) of the generated threshold function initialization.}
#' \item{omeg.hat}{The estimated change plane parameters.}


## function to initialize threshold function estimation
## adapted from Wei and Kosorok, 2018
f0_Fun <- function(X, y) {
  
  p=ncol(X)
  n=nrow(X)
  
  ## get \Omega_0
  km = stats::kmeans(x = X, centers = n/10)
  Omeg0 = NULL
  for (i in which(km$size > 4)) {
    
    if( length(which(km$cluster == i)) > 10 ){
      id = sample(which(km$cluster == i), 10)
    } else {
      id = which(km$cluster == i)
    }
    all.part = partitions::setparts(partitions::restrictedparts(length(id),2)[,-1])
    Omeg0 = cbind(Omeg0, apply(all.part, 2, function(x){
      m = sda::centroids(x = X[id,], L = x, lambda.var=0, var.groups = FALSE, verbose = FALSE)$mean[,1:2]
      vec = matrixStats::rowDiffs(m)
      unit.dir = vec / sqrt(sum(vec^2)) 
      if(unit.dir[1] < 0) unit.dir <- -1*unit.dir 
      
      return(unit.dir)
      
    }))
    
  }
  
  ## get gam
  gam_Fun = function(omeg, gam.seq, X, y){
    
    nobs = nrow(X)
    ngrid <- length(gam.seq)
    ssr <- rep(NA, ngrid)
    for (i in 1:ngrid) {
      
      ind <- Matrix::drop( X%*%omeg >= gam.seq[i] )
      x.reg <- cbind(X, X*ind)
      m = stats::glm.fit(x = x.reg, y = y, intercept = FALSE)
      # bb = m$coefficients
      uhat <- y - m$fitted.values
      ssr[i] <- t(uhat)%*%uhat / (2*nobs)
      
    }
    opt<-which.min(ssr)
    gam.hat <- gam.seq[opt]
    ssr <- ssr[opt]
    
    return(list(gam = gam.hat, ssr = ssr))
    
  }
  
  ## slice invese regression
  discretize = function(y, h) {
    n = length(y)
    m = round(n / h)
    y = y + .00001 * mean(y) * stats::rnorm(n)
    yord = y[order(y)]
    divpt = numeric()
    for (i in 1:(h - 1))
      divpt = c(divpt, yord[i * m + 1])
    y1 = rep(0, n)
    y1[y < divpt[1]] = 1
    y1[y >= divpt[h - 1]] = h
    
    for (i in 2:(h - 1)){
      y1[(y >= divpt[i - 1]) & (y < divpt[i])] = i
    }
    
    return(y1)
  }
  # ##
  # matpower = function(a, alpha) {
  #   a = round((a + t(a)) / 2, 7)
  #   tmp = eigen(a)
  #   return(tmp$vectors %*% diag((tmp$values) ^ alpha) %*% t(tmp$vectors))
  # }
  ##
  Omega_sir <- function(X, y, omeg, gam.seq, nslices = n/10){
    
    p=ncol(X)
    n=nrow(X)
    ##
    gam.hat = gam_Fun(omeg = omeg, gam.seq = gam.seq, X = X, y = y)$gam
    # changeplane.id = (X%*%omeg >= gam.hat)
    ##
    y0 = y[(X%*%omeg < gam.hat)]
    y1 = y[(X%*%omeg >= gam.hat)]
    
    ydis0=discretize(y0, nslices)
    ylabel0 = unique(ydis0)
    ydis1=discretize(y1, nslices)
    ylabel1 = unique(ydis1)
    
    prob1 = prob0 = rep(NA, nslices)
    exy1 = exy0 = matrix(NA, nrow = nslices, ncol = p)
    for (i in 1:nslices){
      prob1[i] = length(ydis1[ydis1 == ylabel1[i]]) / length(y1)
      exy1[i,] = colMeans(X[(ydis1 == ylabel1[i]), ,drop = FALSE])/ prob1[i] - colMeans(X)
      prob0[i] = length(ydis0[(ydis0 == ylabel0[i])]) / length(y0)
      exy0[i,] = colMeans(X[(ydis0 == ylabel0[i]), ,drop = FALSE])/ prob0[i] - colMeans(X)
    }
    sirmat=t(exy1)%*%diag(prob1)%*%exy1 + t(exy0)%*%diag(prob0)%*%exy0;  
    # print(sirmat)
    signrt = solve(stats::var(X))%*%eigen(sirmat)$vectors[,1];   # matpower_cpp(stats::var(X)+stats::diag(lambda,p,p),-1)
    unit.dir = signrt / sqrt(sum(signrt^2)) 
    if(unit.dir[1] < 0) unit.dir <- -1*unit.dir 
    
    return(unit.dir)
  }
  
  ##
  gam.seq = seq(-1,1, length.out = 21)
  cat('\n', 'sieve size:', ncol(Omeg0), '\n')
  
  Omega.hat = sapply(1:ncol(Omeg0), function(k){
    Omega_sir(X = X, y = y, omeg = Omeg0[,k], gam.seq = gam.seq)
  })

  ##
  ssr = NULL
  for (i in 1:ncol(Omega.hat)) {
    # print(i)
    ssr = c(ssr, gam_Fun(omeg = Omega.hat[,i], gam.seq = gam.seq, X = X, y = y)$ssr)
  }
  opt = which.min(ssr)
  omeg.hat = Omega.hat[,opt]
  gam.hat = gam_Fun(omeg = omeg.hat, gam.seq = gam.seq, X = X, y = y)$gam
  
  
  return(list(f0=X%*%omeg.hat, omeg.hat=omeg.hat))
  
}



#' @title pred_Fun
#'
#' @description The function for prediction based on the Change Surface (CSurf) regression results
#'
#'
#' @param y A vector (length n) of the true response.
#' @param X A matrix (n x p) of the covariates. The covariates are also used as potential thresholding variables for subgroup identification.
#' @param f A vector (length n) the threshold function (change surface) obtained by CSurf.
#' @param tau The change points over the threshold functions obtained by CSurf.
#' @param alp The estimated regression coefficients obtained by CSurf.
#'
#' @return A list containing the following fields: 
#' \item{Yhat}{A vector of the predicted response.}
#' \item{mse}{The mean squared error (MSE) of response prediction.}


pred_Fun <- function(X, y, f, tau, alp){
  
  n <- nrow(X)
  p <- ncol(X)
  s <- length(tau)
  
  if (s>=1){
    
    tmp <- NULL
    
    for (j in 1:s) {
      
      delta_j = alp[(j*p+1) : ((j+1)*p)]
      ind = which(f>tau[j])
      ktemp = matrix(0, nrow=n, ncol=1)
      ktemp[ind,] = 1
      tmp[[j]] <- X %*% delta_j * ktemp 
      
    }

    sm <- Reduce("+", tmp)
    Yhat <- X%*%alp[1:p]+sm
    
  } else {
    
    Yhat <- X%*%alp
    
  }
  
  mse = mean((y-Yhat)^2)
  
  return(list(Yhat=Yhat, mse=mse))
  
}



#################################
#################################
### Define auxiliary functions###
#################################
#################################
hausdorff.distance=function(v1,v2,n=1){
  p1=length(v1)
  p2=length(v2)
  if ( p1*p2==0){ dis=n} else{
    distance.mat=matrix(0,nrow=p2,ncol=p1)
    for( i in 1: p2){
      for( j in 1: p1){
        distance.mat[i,j]=abs(v2[i]-v1[j])
        
      }
    }
    # dis=max(max( apply(distance.mat, 1, min)),max( apply(distance.mat, 2, min)))
    dis= max( apply(distance.mat, 1, min))
  }
  return(dis)
  
}

##
lambdaseq <- function(X, y, groupindex, lambdaRatio = 0.01, nLambda = 20, addZeroLambda = FALSE){
  
  n = nrow(X)
  ng = length(unique(groupindex))
  group = lambda = NULL
  for (tt in 1:ng) {
    group[[tt]] = which(groupindex == tt)
    lambda = c(lambda, norm(t(X[,group[[tt]]])%*%y/n))
  }
  
  lambdaMax <-  max(lambda)
  lambdaMin <- lambdaMax * lambdaRatio
  
  loghi <- log(lambdaMax)
  loglo <- log(lambdaMin)
  logrange <- loghi - loglo
  interval <- -logrange/(nLambda-1)
  
  lambda <-  exp(seq.int(from = loghi, to = loglo, by = interval))
  if(addZeroLambda == T){
    lambda[length(lambda)] <-  0
  }else{
    lambda[length(lambda)] <-  lambdaMin;
  }
  
  return(list("lambda" = lambda, "lambdaMin" = lambdaMin,
              "lambdaMax" = lambdaMax))
}

##
norm <- function(beta){
  temp <- sqrt(sum(beta^2))
  return(temp)
}
##
norm_f <- function(f){
  temp <- sqrt(mean(f^2))
  return(temp)
}
## SCAD penalty function
SCAD <- function(x, lambda) {
  
  x <- abs(x)
  # x <- norm_f(x)
  a <- 3.7
  
  u <- (x<=lambda)
  penalty <- lambda*x*u + (-0.5*pmax(a*lambda-x,0)^2/(a-1)  + 0.5*(a+1)*lambda^2)*(1-u)
  return(penalty)
}
## first order derivative of SCAD_f penalty
SCAD_f.dev <- function(x, lambda) {
  
  # x <- abs(x)
  x <- norm_f(x)
  a <- 3.7
  
  u <- (x<=lambda)
  penalty.derivative <- lambda * u + (pmax(a*lambda-x,0)/(a-1)) * (1-u)
  return(penalty.derivative)
}
##
obj_f_Fun <- function(f_hat, y, X, alp, tau, h){
  
  
  f = rowSums(f_hat)
  K = length(tau)
  n = dim(X)[1]
  p = dim(X)[2]
  
  beta = alp[1:p]
  delta = alp[-c(1:p)]
  
  tmp = NULL
  for (k in 1:K) {
    delta_k = alp[(k*p+1) : ((k+1)*p)]
    tmp[[k]] = c(X%*%delta_k)* stats::pnorm(c(f - tau[k])/h)
  }
  tmp = Reduce("+", tmp)
  ##
  obj = mean((y - X%*%beta - tmp)^2)/2
  
  return(obj)
  
}

##
grad_Fun <- function(f_hat, y, X, alp, tau, h){
  
  ## cubic spline basis 
  f = rowSums(f_hat)
  
  K = length(tau)
  n = dim(X)[1]
  p = dim(X)[2]
  
  beta = alp[1:p]
  delta = alp[-c(1:p)]
  
  tmp1 = tmp2 = NULL
  for (k in 1:K) {
    delta_k = alp[(k*p+1) : ((k+1)*p)]
    tmp1[[k]] = c(X%*%delta_k)* stats::dnorm(c(f - tau[k])/h)
    tmp2[[k]] = c(X%*%delta_k)* stats::pnorm(c(f - tau[k])/h)
  }
  tmp1 = Reduce("+", tmp1)
  tmp2 = Reduce("+", tmp2)
  
  ##
  grad_gam = (-1/h)*(y - X%*%beta - tmp2)*tmp1
  
  return(c(grad_gam))
  
}


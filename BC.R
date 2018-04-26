#=====================================================================
#estVAR(X, p, demean, intercept)
#=====================================================================
est_VAR <- function(X, p=1, demean=FALSE, intercept=TRUE){
  
  # function [Gamma_hat, alpha_hat, Omega_hat] = regressVAR(X)
  # inputs:
  # X: T*N data matrix
  # p: order of the VAR (default 1)
  # demean: demean the data before estimation? (default false)
  # intercept: include an intercept in the regression? (default true)
  #
  # output:
  # Gamma_hat : N*N
  # alpha_hat : N*1
  # Omega_hat : N*N
  #
  # X(t+1) = alpha + Gamma*X(t) + eps(t+1), cov(eps(t+1)) = Omega
  #
  # Compute the maximum likelihood estimates of Gamma, alpha, Omega
  #
  # NOTE: The MLE estimates of Gamma, alpha do not depend on Omega.
  # That is, the argmax_{Gamma,alpha} [L(X|Gamma,alpha,Omega)] = f(X)
  # So this function provides MLE of Gamma, alpha for a fixed Omega.
  
  # parameters/flags
  T <- dim(X)[1]
  k <- dim(X)[2]
  
  if(demean==TRUE){
    X <- X - matrix(1,T,k)*mean(X)
  }
  
  if(p==1){
    
    Yt <- X[1:nrow(X)-1,]  #(T-1)*k
    Ytp1 <- X[2:nrow(X),]  #(T-1)*k
    Y <- t(Ytp1)   #(T-1)*k
    
    if(intercept==TRUE){
      
      Z <- t(cbind(matrix(1,T-1,1), Yt))
      
    } else {
      Z <- t(Yt)     #(T-1)*k
    }
    
    A <- Y%*%t(Z)%*%solve(Z%*%t(Z))  #k*(k+1) / k*k
    
    if(intercept==TRUE){
      
      alpha_hat <- A[,1]
      Gamma_hat <- A[,2:ncol(A)]
      
    } else {
      
      alpha_hat <- NaN
      Gamma_hat <- A
      
    }
    
    residuals = Ytp1 - t(A%*%Z)         #(T-1)*N
    Omega_hat = 1/(T-1)*t(residuals)%*%residuals
    
  } else if(p==2){
    
    Ytm2 <- X[1:nrow(X)-2,]  #(T-2)*k
    Ytm1 <- X[2:nrow(X)-1,]  #(T-2)*k
    Y <- t(X[3,nrow(X),])    #k*(T-2)
    
    if(intercept==TRUE){
      
      Z <- t(cbind(matrix(1,T-2,1), Ytm1, Ytm2))     #(2*k+1)*(T-2)
      
    } else {
      Z <- t(cbind(Ytm1, Ytm2))     #2k*(T-2)
    }
    
    A <- Y%*%t(Z)%*%solve(Z%*%t(Z))  #k*(2k+1) / 2k*2k
    
    if(intercept==TRUE){
      
      alpha_hat <- A[,1]
      Gamma_hat <- A[,2:ncol(A)]
      
    } else {
      
      alpha_hat <- NaN
      Gamma_hat <- A
      
    }
    
    Gamma_hat <- cbind(Gamma_hat, diag(k), matrix(0,k,k))
    residuals = Y - t(A%*%Z)         #(T-2)*N
    Omega_hat = 1/(T-2)*t(residuals)%*%residuals
    
  } else {stop('not implemented for order>2')}
  
  return(list(Gamma_hat=Gamma_hat, alpha_hat=alpha_hat, Omega_hat=Omega_hat))
}


#=====================================================================
#genVAR(Phi, M, Y, p)
#=====================================================================
genVAR <- function(Phi, M, Y, p=1){
  
  #generate M data sets from VAR(p) model
  T <- dim(Y)[1]
  k <- dim(Y)[2]
  
  Y_mean <- colMeans(Y)
  Y_mean_rep <- Y_mean%*%matrix(1,1,M)
  X_sim <- array(0, c(k, T, M))
  
  if(p==1){
    
    # VAR(1)
    # 1. obtain residuals
    # use bootstrapped residuals
    resid <- matrix(0, k, T-1)
    
    for(t in 2:T){
      temp1 <- as.matrix(t(Y[t,]) - Y_mean)
      temp2 <- as.matrix(Y[t-1,] - Y_mean)
      resid[,t-1] <- t(temp1) - (Phi%*%temp2)
    }
    
    # 2. generate series
    # randomly select initial values from the data Y
    ind_start <- sample(T,M,replace = TRUE)
    X_sim[,1,] <- t(Y[ind_start,])
    
    for (t in 2:T) {
      u_sim <- resid[,sample(T-1,M,replace = TRUE)]
      X_sim[,t,] <- Y_mean_rep + Phi %*% (aa=as.matrix(X_sim[,t-1,],k,M) - Y_mean_rep) + u_sim
      
    }
    
  } else if(p==2){
    
    # VAR(2)
    Phi1 <- Phi[1:k,1:k]
    Phi2 <- Phi[1:k, k+1:2*k]
    
    # randomly select initial values from the data Y
    ind <- sample(T-1,M,TRUE)
    X_sim[,1,] <- t(Y[ind,])
    X_sim[,2,] <- t(Y[ind+1,])
    
    # use bootstrapped residuals
    # 1. obtain residuals
    resid <- matrix(0, k, T-2)
    
    for(t in 3:T){
      temp1 <- as.matrix(t(Y[t,]) - Y_mean)
      temp2 <- as.matrix(Y[t-1,] - Y_mean)
      resid[,t-2] <- temp1 - (Phi1%*%temp2) - (Phi2%*%temp2) 
    }
    
    # 2. generate series
    for (t in 3:T) {
      X_sim[,t,] <- Y_mean_rep + Phi1%*%(as.matrix(X_sim[,t-1,],k,M) - Y_mean_rep) 
                          + Phi2%*%(as.matrix(X_sim[,t-2,],k,M) - Y_mean_rep) + resid[,sample(T-2,M,replace = TRUE)]
    }
  } else {stop('not implemented for p>2')}
  
  return(X_sim)
}

#=====================================================================
#shrink_phi(Phi_tilde, Phi_hat, ev_restr, trace)
#=====================================================================
shrink_phi <- function(Phi_tilde, Phi_hat, ev_restr, trace=1){
  
  # stationarity adjustment
  # see Kilian (1998, REStat) "Small-sample confidence intervals for impulse response functions"
  # impose eigenvalue restriction
  ev_bc <- max(abs(eigen(Phi_tilde)$values))     # get scalar largest absolute eigenvalues
  ev_ols <- max(abs(eigen(Phi_tilde)$values))
  
  if(ev_bc > ev_restr){                # impose restriction:  if ev_bc is larger than ev_restr
    if(trace>0){
      fprintf('*** eigenvalue constraint binding ***\n')
      fprintf('*** constraint: %8.6f \n', ev_restr)
    }
    
    if(ev_ols<ev_restr){                     # and if Phi_hat actually satisfies eigenvalue restr
                                             # then shrink Phi_tilde to Phi_hat
      Phi_diff <- Phi_hat - Phi_tilde        # difference/bias adjustment -- this is what we shrink to zero
      delta <- 1 
      
      # as long as restriction is not satisfied, go through the following loop
      while (ev_bc > ev_restr){
        
      delta <- delta - .01                        # reduce scaling factor
      Phi_diff <- delta*Phi_diff                  # shrink Phi_diff to zero
      Phi_tilde <- Phi_hat - Phi_diff             # try this value of Phi_diff
      ev_bc <- max(abs(eig(Phi_tilde)))           # calculate largest eigenvalue
        
      }
    } else {
      
      # OLS estimates do not satisfy restriction
      # warning('ev_ols > ev_restr -- setting Phi_bc equal to Phi_ols')
      Phi_tilde <- Phi_hat
      
    }
    
    if(trace > 0){
      fprintf('largest absolute eigenvalue after imposing eigenvalue restriction:  %8.6f \n', ev_bc)
    }
  }
  
  return(Phi_tilde)
}

#=====================================================================
# m_var(theta, M, p, Y, flag_mean)
#=====================================================================
m_var <- function(theta, M, p, Y, flag_mean=TRUE){
  
  # find mean/median of OLS when DGP is VAR(p)
  # works for VAR(1) and VAR(2) only
  # inputs:
  #  theta - vector of parameters determining Phi_1 (and Phi_2)
  #  M - number of Monte Carlo replications
  #  p - order of the VAR
  #  Y - data, T x k matrix
  #  flag_mean - flag whether mean (TRUE) or median (FALSE) is to be returned
  
  TT <- dim(Y)[1]
  k <- dim(Y)[2]
  
  if(p==1){
    
    # get coefficient matrix
    Phi_tilde <- matrix(theta, k,k, byrow = FALSE)
    
    # simulate M datasets
    X_sim <- genVAR(Phi_tilde, M, Y, p)
    
    # estimate VAR(1) on each series
    theta_new_i <- matrix(0, M, k^2)
    
    for (m in 1:M) {
      Phi_new <- est_VAR(t(X_sim[,,m]), 1, TRUE, FALSE)$Gamma_hat
      theta_new_i[m,] <- Phi_new
    }
    
    if(flag_mean==TRUE){
      # find mean of OLS estimates
      Phi_new <- matrix(colMeans(theta_new_i), k, k, byrow = FALSE)
    } else {
      # find median of OLS estimates
      Phi_new <- matrix(apply(theta_new_i, 1, median), k, k, byrow = FALSE)
    }
    
    Phi_sample <- theta_new_i
    
  } else if(p==2){
    
    # get coefficient matrix
    Phi_tilde <- matrix(theta, k,2*k, byrow = FALSE)
    
    # simulate M datasets
    X_sim <- genVAR(Phi_tilde, M, Y, p)
    
    # estimate VAR(2) on each series
    theta_new_i <- matrix(0, M, 2*k^2)
    
    for (m in 1:M) {
      F_new <- est_VAR(t(X_sim[,,m],2, TRUE, FALSE))
      temp <-F_new[1:k,]
      theta_new_i[m,] <- temp
    }
    
    if(flag_mean==TRUE){
      # find mean of OLS estimates
      Phi_new <- matrix(colMeans(theta_new_i), k, 2*k, byrow = FALSE)
    } else {
      # find median of OLS estimates
      Phi_new <- matrix(apply(theta_new_i, 1, median), k, 2*k, byrow = FALSE)
    }
    
  } else {
    stop("not yet implemented for order>2")
  }
  
  return(list(Phi_new=Phi_new, Phi_sample=Phi_sample))
}


#=====================================================================
# est_bc_var(X, ev_restr, p, flag_mean, N, N_burn, B, check, B_check)
#=====================================================================
est_bc_var <- function(X, ev_restr=1, p=1, flag_mean=TRUE, N=5000, N_burn=1000, B=10, check=TRUE, B_check=100000){
  
  # est_bc_var - bias-corrected VAR estimation using stochastic approximation (SA)
  # see Bauer, Rudebusch, Wu (2012, JBES) Correcting Estimation Bias in Dynamic
  # Term Structure Models
  # inputs:
  #  X       REQUIRED. data matrix, T x k
  #  ev_restr  largest eigenvalue not to be exceeded by Phi_tilde
  #  p       order of the VAR, default VAR(1), p=2 possible, p>2 not implemented
  #  flag_mean     flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired
  #           default TRUE (mean-unbiased)
  #  N       number of iterations of the SA algorithm after burn-in (default 10,000)
  #  N_burn  number of burn-in iterations (default 100)
  #  B       number of bootstrap samples to calculate noisy measure of mean/median
  #          of the OLS estimator (default 50)
  #  check   flag whether closeness check is to be performed in the end (default TRUE)
  #  B_check number of bootstrap samples for closeness check (default 100,000)
  #
  # outputs:
  #  Phi_tilde  mean-reversion matrix
  #  mu_tilde   intercept
  #  V_tilde    variance-covariance matrix of residuals
  gamma_i <- 0.2
  
  if(flag_mean==TRUE){
    fprintf('Mean-')
  } else {
    fprintf('Median-')
  }
  
  fprintf('Bias-adjustment for VAR estimates\n')
  cat("N = ", N, "N_burn = ", N_burn, "B = ", B, "B_check = ", B_check, "p = ", p, "\n")
  
  TT <- dim(X)[1]
  k <- dim(X)[2]
  
  if(p==1){
    
    # first-order VAR

    # OLS
    Phi_hat <- est_VAR(X,1,TRUE,FALSE)$Gamma_hat
    cat("largest absolute eigenvalue OLS: ", max(abs(eigen(Phi_hat)$values)), "\n")
    
    # initialization for SA algorithm
    theta <- matrix(0, k^2, N_burn+N)
    theta_hat <- Phi_hat
    theta[,1] <- theta_hat              # starting value
    
    # SA algorithm
    for (j in 1:(N_burn+N-1)) {     
      
      Phi_new <- m_var(theta[,j],B,1,X,flag_mean)$Phi_new
      theta_new <- Phi_new
      d <- theta_hat - theta_new
      theta[,j+1] <- theta[,j] + gamma_i*d
      
    }
 
    theta_tilde <- apply(theta[,(N_burn+1):(N_burn+N)],1,mean)
    
    if(check==TRUE){
      
      # check whether mean/median of OLS---given that DGP is VAR with Phi_tilde---is close to actual OLS estimates
      cat("... checking closeness of mean/median to actual estimates ...")
      
      Phi_new <- m_var(theta_tilde, B_check, 1, X, flag_mean)$Phi_new
      dist <- sqrt(sum((theta_new - theta_hat)^2)/length(theta_new))
      fprintf("root mean square distance: %8.6f \n", dist)
      
    }
    
    # bias_adjusted estimates
    Phi_tilde <- matrix(theta_tilde,k,k, byrow = FALSE)
    ev_bc <- max(abs(eigen(Phi_tilde)$values))                       # get scalar largest absolute eigenvalues
    fprintf("largest eigenvalue after BC:  %8.6f \n", ev_bc)
    
    # impose restriction on eigenvalue
    Phi_tilde <- shrink_phi(Phi_tilde, Phi_hat, ev_restr)
    
    # choose intercept to match sample mean
    mu_tilde <- (diag(k) - Phi_tilde)%*%matrix(colMeans(X),5,1)
    
    # residuals and their variance-covariance matrix
    Xdem <- X - matrix(1,TT,1)%*%colMeans(X)
    resid_tilde <- t(Xdem[2:TT,]) - Phi_tilde %*% t(Xdem[1:(TT-1),])
    V_tilde <- resid_tilde %*% t(resid_tilde)/(TT-1)
  } else {
    
    stop("bias-corrected VAR estimation only implemented for first-order VAR")
  }
  
  return(list(Phi_tilde=Phi_tilde, mu_tilde=mu_tilde, V_tilde=V_tilde))
}



#=====================================================================
# Estimation of Term Structure Model
#=====================================================================
DTSM_est <- function(score, W, m, m_m, initial=NULL, ns, specId){
  
  score <- score[,1:3]             # PCs
  bigt <- dim(score)[1]            # sample size
  
  # estimate short rate equation by OLS
  x <- cbind(matrix(1,bigt,1), score)  
  delta <- ols(x,m[,1]/4)$bhat        # regress short yield on PCs
  vcov3 <- ols(x,m[,1]/4)$vcov
  delta0 <- delta[1]
  delta1 <- delta[2:4]
  mats <- ns/4                        # mats in years
  score_m_m <- cbind(score,m_m)
  
  if(specId==1){ #OLS
    
    cat("OLS               \n")
    temp0 <- varest(score_m_m,1,TRUE)
    bhat <- temp0$bhat
    sigma <- temp0$sigma
    phi <- bhat[,2:ncol(bhat)]
    mu.temp <- colMeans(score_m_m)
    mu <- matrix(mu.temp,length(mu.temp),1) - phi%*%matrix(mu.temp,length(mu.temp),1)
    cat("largest P-eigenvalue: ", num2str(max(abs(eig(phi)))), "\n")
    
  } else if(specId==2){ # BC2 -- indirect inference bias correction
    
    cat("            \n")
    cat("BC2 -- indirect inference bias correction               \n")
    temp1 <- est_bc_var(score_m_m,1,1,1,5000,1000,5,0)
    phi <- temp1$Phi_tilde
    mu <- temp1$mu_tilde
    sigma <- temp1$V_tilde
    
  } else {   # BC3 -- bootstrap bias correction
    
    cat("            \n")
    cat("BC3 -- bootstrap bias correction             \n")
    bhat <- varest(score_m_m,1,TRUE)$bhat
    phi_ols <- bhat[,2:ncol(bhat)]
    cat('largest eigenvalue OLS: ', num2str(max(abs(eigen(phi_ols)$values))), "\n")
    
    # bootstrap bias correction
    phi <- 2*phi_ols - m_var(phi_ols, 5000, 1, cbind(score, m_m),TRUE)$Phi_new
    cat('largest eigenvalue BC: ', num2str(max(abs(eigen(phi)$values))))
    # impose eigenvalue restriction
    phi <- shrink_phi(phi, phi_ols, 1, 1)
    mu <- (diag(5)-phi)%*%matrix(colMeans(score_m_m),ncol(score_m_m),1)
    # residual covariance matrix
    Xdem <- score_m_m - matrix(1,bigt,1)%*%colMeans(score_m_m)
    resids <- t(Xdem[2:bigt,]) - phi%*%t(Xdem[1:(bigt-1),])
    sigma <- resids %*% t(resids)/(bigt-1)
    
  }
  
  sigma <- t(chol(sigma))      # Cholesky decomposition of cov matrix (lower triangular)
  
  # numerical optimization for cross-sectional parameters
  if(specId==1){
    
    # OLS: starting values as Wright's code
    temp2 <- c(-0.01,-0.04,-0.08,0.08,sigma[1:5,1],sigma[2:5,2],sigma[3:5,3],sigma[4:5,4],sigma[5,5])
    paramstart <- matrix(temp2,length(temp2),1)
    
  } else {
    
    # BC: use OLS estimates as starting values
    paramstart <- initial
      
  }
  
  # parameter restrictions imposed within jpslikel
  # Wright's restrictions:
  # 1. lamQ must have difference of at least .01 and sorted descending
  # 2. lamQ must be between -.002 and -.6
  # 3. rinfQ must be between 0 and .4
  # here slightly weaker to get better fit:
  # - lamQ's arbitrarily close to each other
  # - lamQ's maximum -.0001
    
  result <- nlm(jpslikel, paramstart, m=m, W=W, mats=mats, Sigma_cP_OLS=sigma%*%t(sigma))
  param <- result$estimate
  code <- result$code
  
  if(code==1){
    cat("relative gradient is close to zero, current iterate is probably solution.\n")
  } else if(code==2){
    cat("successive iterates within tolerance, current iterate is probably solution. \n")
  } else if(code==3){
    cat("last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.\n")
  } else if(code==4){
    cat("iteration limit exceeded.\n")
  } else {
    cat("maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.\n")
  }
  
  ch <- matrix(c(param[5],0,0,param[6],param[10],0,param[7],param[11],param[14]),3,3,byrow = TRUE)
  temp3 <- jszLoadings(W, diag(param[1:3]), param[4], ch%*%t(ch), mats, 1/4)
  
  BcP <- temp3$BcP
  AcP <- temp3$AcP
  AX <- temp3$AX
  BX <- temp3$BX
  
  temp4 <- jszRotation(W, diag(param[1:3]), param[4], 1/4, NULL, NULL, BX, AX)
  K0Q_cP <- temp4$K0Q_cP
  K1Q_cP <- temp4$K1Q_cP
  rho0_cP <- temp4$rho0_cP
  rho1_cP <- temp4$rho1_cP
  
  cP <- m%*%t(W)            # risk factors -- only yields
  mfit <- (matrix(1,bigt,1)%*%AcP)+(cP%*%BcP)
    
  # fitting errors
  u <- mfit-m
  cat("Root Mean square fitting error (percentage points): ", sqrt(mean(colMeans((100*u)^2))), "\n")

  # save parameters
  K1Q_X <- diag(param[1:3])
  rinfQ <- param[4]
  
  # average expected short rate
  avexp = get_exp(score_m_m, delta0, delta1, mu, phi)
  
  return(list(param=param, mu=mu, phi=phi, sigma=sigma, K1Q_X=K1Q_X, rinfQ=rinfQ, ch=t(ch)))
}



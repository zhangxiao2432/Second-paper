#===============================================================================================
#gaussianDiscreteYieldLoadingsDiagonal(maturities, K0d, K1d_diag, H0d, rho0d, rho1d, timestep)
#===============================================================================================
gaussianDiscreteYieldLoadingsDiagonal <- function(maturities, K0d, K1d_diag, H0d, rho0d, rho1d, timestep){
 
  #function [By,Ay, dAyH0d] = gaussianDiscreteYieldLoadingsDiagonal(maturities, K0d, K1d_diag, H0d, rho0d, rho1d, timestep)
  #
  # DOESN'T HANDLE UNIT ROOTS!!
  #
  # THIS FUNCTION ASSUMES K1d is diagonal
  # K0d      : N*1
  # K1d_diag : N*1
  # H0d      : N*N
  # rho0d    : scalar  
  # rho1d    : N*1
  # timestep : optional argument.
  #
  # By : N*M
  # Ay : 1*M  (faster to not compute with only one output argument)
  #
  # r(t) = rho0d + rho1d'Xt
  #      = 1 period discount rate
  # P(t) =  price of  t-period zero coupon bond
  #      = EQ0[exp(-r0 - r1 - ... - r(t-1)]
  #      = exp(A+B'X0)
  # yields = Ay + By'*X0
  #   yield is express on a per period basis unless timestep is provided.
  #   --For example, if the price of a two-year zero is exp(-2*.06)=exp(-24*.005),
  #   --and we have a monthly model, the function will return Ay+By*X0=.005
  #   --unless timestep=1/12 is provided in which case it returns Ay+By*X0=.06
  #
  # Where under Q:
  #   X(t+1) - X(t) = K0d + K1d*X(t) + eps(t+1),  cov(eps(t+1)) = H0d
  #
  # We can compute the loadings by recurrence relations:
  #   A1 = -rho0d
  #   B1 = -rho1d
  #   At = A(t-1) + K0d'*B(t-1) .5*B(t-1)'*H0d*B(t-1) - rho0d
  #   Bt = B(t-1) + K1d'*B(t-1) - rho1d
  #
  # Or in closed form by noting that 
  #    r0+r1+..+r(t-1) = c.X(0) + alpha0 + alpha1*eps1 + ... + alpha(t-1)*eps(t-1)
  #                    ~ N(c.X(0) + alpha0, alpha1'H0d*alph1 + ... + alpha(t-1)'*H0d*alpha(t-1))
  #
  # And then use the MGF of Y~N(mu,Sigma) is E[exp(a.Y)] = a'*mu + .5*a'*Sigma*a
  # (or similarly use the partial geometric sum formulas repeatedly)
  #
  # Let G = K1+I
  # X(0)
  # X(1) = K0 + G*X(0) + eps1
  # X(2) = K0 + G*K0 + G^2*X(0) + G*eps1 + eps2
  # X(3) = K0 + G*K0 + G^2*K0 + G^3*X(0) + G^2*eps1 + G*eps2 + eps3
  # X(n) = sum(I+G+..+G^(n-1))*K0 + G^n*X(0) + sum(i=1..n,G^(n-i)*epsi)
  #      = (I-G\(I-G^n)*K0 + G^n*X(0) + sum(i=1..n,G^(n-i)*epsi)
  #
  # cov(G^n*eps) = G^n*cov(eps)*(G^n)'
  # vec(cov(G^n*eps) = kron(G^n,G^n)*vec(eps)
  #                   = (kron(G,G)^n)*vec(eps)
  #
  # sum(X(i),i=1..n) = mu0 + mu1*X0 + u
  #    mu0 = (I-G)\(I - (I-G)\(G-G^(n+1))*K0 
  #    mu1 = (I-G)\(G-G^(n+1))
  #    vec(cov(u)) = see below.  
  #  u = (I-G)\(I-G^n)*eps1 + 
  #      (I-G)\(I-G^(n-1))*eps2 + ..
  #      (I-G)\(I-G)*epsn
  #  cov(u) = (I-G)\Sig/(I-G)'
  # Sig = sum(cov(eps)) + sum(i=1..n,G^i*cov(eps)) + 
  #       sum(i=1..n,cov(eps)G^i') + sum(i=1..n,G^i*cov(eps)*G^i')
  # compute the last one using vec's.  see below.
  
  # K1d_diag is N*1 -- the diagonal of K1d
  M <- length(maturities)
  N <- length(K0d)
  Ay <- matrix(0,1,M)
  By <- matrix(0,N,M)
  dAyH0d <- array(rep(0, N*N*M),c(N,N,M))
  
  I <- matrix(1,N,1)
  G <- K1d_diag+matrix(1,N,1)                            # N*1
  
  GG <-  G %*% t(G)                                         # N*N
  GpG <- G %*% matrix(1,1,N) + matrix(1,N,1)%*%t(G)         # N*N (i,j) entry is G(i)+G(j)

  for (m in 1:M) {
    
    mat <- maturities[m]
    
    if(mat==1){
      
      By[,m] <- rho1d
      Ay[,m] <- rho0d
      mu0 <- matrix(0,N,1)
      mu1 <- matrix(1,N,1)
      Sigma0 <- matrix(0,N,1)
      
    }
    
    i <- mat-1                                               # # of random innovations X(0) + X(1) + ... + X(mat)=X(0)+...+X(i)
                                                             # X(0) + ... + X(i)~N(mu0 + mu1*X(0), Sigma0)
    mu1 <- I -  (G - G^(i+1))/K1d_diag             # N*1
    mu0 <- -((i+1)*I - mu1)*K0d/K1d_diag
    By[,m] <- rho1d * mu1 /mat
    
    Sigma_term1 <- i*H0d                                   # N*N
    Sigma_term2 <- ((mu1 - 1)%*%matrix(1,1,N))%*%H0d         # N*N
    Sigma_term3 <- t(Sigma_term2)
    Sigma_term4 <- (GG - GG^(i+1))*H0d/(1 - GG) 
    Sigma_term4 <- matrix(Sigma_term4,N,N, byrow = TRUE)
    Sigma0 <- (Sigma_term1 - Sigma_term2 - Sigma_term3 + Sigma_term4)/K1d_diag%*%t(K1d_diag)

    Ay[,m] <- rho0d + (t(rho1d)%*%mu0 - .5%*%t(rho1d)%*%Sigma0%*%rho1d)/mat
    
    dAyH0d[,,m] <- -.5*i*rho1d%*%t(rho1d)+.5*((mu1-1)%*%matrix(1,1,N))*(rho1d%*%t(rho1d))+.5*(matrix(1,N,1)%*%t(mu1-1))*(rho1d%*%t(rho1d))-
           .5*matrix((GG - GG^(i+1))/(1 - GG),N,N,byrow = TRUE)*(rho1d%*%t(rho1d))
    dAyH0d[,,m] <- dAyH0d[,,m]/mat/(K1d_diag%*%t(K1d_diag))
    
  }
  
  By <- By/timestep
  Ay <- Ay/timestep
  dAyH0d <- dAyH0d/timestep
  
  return(list(By=By,Ay=Ay,dAyH0d=dAyH0d))
}

#===========================================================================================================
#gaussianDiscreteYieldLoadingsRecurrence(mats_periods, K0Q_X, K1Q_X, matrix(0,N,N), rho0d*dt, rho1d*dt, dt)
#===========================================================================================================
gaussianDiscreteYieldLoadingsRecurrence <- function(maturities, K0d, K1d, H0d, rho0d, rho1d, timestep){
  
  # function [By, Ay] = gaussianDiscreteYieldLoadingsRecurrence(maturities, K0d, K1d, H0d, rho0d, rho1d, timestep)
  #
  # K0d      : N*1
  # K1d      : N*1
  # H0d      : N*N
  # rho0d    : scalar  
  # rho1d    : N*1
  # timestep : optional argument.
  #
  # By : N*M
  # Ay : 1*M  (faster to not compute with only one output argument)
  #
  # r(t)   = rho0d + rho1d'Xt
  #        = 1 period discount rate
  # P(t)   =  price of  t-period zero coupon bond
  #        = EQ0[exp(-r0 - r1 - ... - r(t-1)]
  #        = exp(A+B'X0)
  # yields = Ay + By'*X0
  #   yield is express on a per period basis unless timestep is provided.
  # --For example, if the price of a two-year zero is exp(-2*.06)=exp(-24*.005),
  # --and we have a monthly model, the function will return Ay+By*X0=.005
  # --unless timestep=1/12 is provided in which case it returns Ay+By*X0=.06
  #
  # Where under Q:
  #   X(t+1) - X(t) = K0d + K1d*X(t) + eps(t+1),  cov(eps(t+1)) = H0d
  #
  # A1 = -rho0d
  # B1 = -rho1d
  # At = A(t-1) + K0d'*B(t-1) + .5*B(t-1)'*H0d*B(t-1) - rho0d
  # Bt = B(t-1) + K1d'*B(t-1) - rho1d
  #
  # mautirities: 1*M # of periods

  M <- length(maturities)
  N <- length(K0d)
  Atemp <- 0
  Btemp <- matrix(0,N,1)
  Ay <- matrix(NaN,1,M)
  By <- matrix(NaN,N,M)

  curr_mat <- 1
  
  for (i in 1:maturities[M]){
    
    Atemp <- Atemp + t(K0d)%*%Btemp +.5*t(Btemp)%*%H0d%*%Btemp - rho0d
    Btemp <- Btemp + matrix(t(K1d)%*%Btemp,N,1) - rho1d
  
    if (i==maturities[curr_mat]){
      
      Ay[1,curr_mat] <- -Atemp/maturities[curr_mat]
      By[,curr_mat] <- -Btemp/maturities[curr_mat]
      curr_mat <- curr_mat + 1
      
    }
  }
  
  Ay = Ay/timestep
  By = By/timestep
  
  return(list(By=By,Ay=Ay))
}


#=====================================================================
#jszLoadings(W, K1Q_X, rinfQ, Sigma_cP, mats, dt, Sigma_X=NULL)
#=====================================================================
jszLoadings <- function(W, K1Q_X, rinfQ, Sigma_cP, mats, dt, Sigma_X=NULL){
  
  #=============================================================================
  # Inputs:
  #   mats       : 1*J,      maturities in years
  #   dt         : scalar,   length of period in years
  #   W          : N*J,      vector of portfolio weights to fit without error.
  #   K1Q_X      : N*N
  #   rinfQ      : scalar,   the long run mean under Q of the annualized short rate
  #   Sigma_cP, Sigma_X : N*N  covariance of innovations. PROVIDE ONE OR THE OTHER
  #
  # Returns:
  #   AcP : 1*J
  #   BcP : N*J
  #   AX  : 1*J
  #   BX  : N*J
  #   Sigma_X : N*N
  #
  #
  # This function:
  # 1. Compute the loadings for the normalized model:
  #     X(t+1) - X(t) = K1Q_X*X(t) + eps_X(t+1), cov(eps_X)=Sigma_X
  #     and r(t) = rinfQ + 1.X(t)  
  #     where r(t) is the annualized short rate, (i.e. price of 1-period zero coupon bond at time t is exp(-r(t)*dt))
  #    If Sigma_X is not provided, it is solved for so that Sigma_cP (below) is matched.
  #    yt = AX' + BX'*Xt
  #  
  # 2. For cPt = W*yt and the model above for Xt, find AcP, BcP so that
  #    yt = AcP' + BcP'*cPt
  #=============================================================================
  J <- length(mats)
  N <- dim(K1Q_X)[1]
  K0Q_X <- matrix(0,N,1)
  rho0d <- rinfQ
  rho1d <- matrix(1,N,1)
  mats_periods <- round(mats/dt)
  M <- max(mats_periods)
  
  K1Q_X <- jszAdjustK1QX(K1Q_X)$K1Q_X
  isTypicalDiagonal <- jszAdjustK1QX(K1Q_X)$isTypicalDiagonal
  
  #=============================================================================
  # If Sigma_cP is provided, we need to compute Sigma_X by 
  # first computing BX
  #
  
  if (is.null(Sigma_X)==TRUE){
  # First compute the loadings ignoring the convexity term -- BX will be correct
  # yt = AX' + BX'*Xt  
  # yt is J*1
  # AX is 1*J
  # BX is N*J
  # Xt is N*1
  #
  # cPt = W*yt  (cPt N*1, W is N*J) 
  #     = W*AX' + W*BX'*Xt
  #     = WAXp + WBXp*Xt
  #
  # Substituting:
  # yt = AX' + BX'*(WBXp\(cPt - WAXp))
  #    = (I - BX'*(WBXp\WAXp))*AX' + BX'*WBXp\cPt
  #    = AcP' + BcP'*cPt
  # where AcP = AX*(I - BX'*(WBXp\WAXp))'
  #       BcP = (WBXp)'\BX
  #
  # Sigma_cP = W*BX'*Sigma_X*(W*BX')'
  # Sigma_X = (W*BX')\Sigma_cP/(W*BX')'
  #
    
  # If K1d isn't diagonal, we should use the Recurrence solver:.
    if (isTypicalDiagonal==TRUE){
      
      temp1 <- gaussianDiscreteYieldLoadingsDiagonal(mats_periods, K0Q_X, diag(K1Q_X), matrix(0,N,N), rho0d*dt, rho1d*dt, dt)  # N*J
      AX <- temp1$Ay
      BX <- temp1$By
      dAyH0d <- temp1$dAyH0d
      
    } else {
      
      temp2 <- gaussianDiscreteYieldLoadingsRecurrence(mats_periods, K0Q_X, K1Q_X, matrix(0,N,N), rho0d*dt, rho1d*dt, dt)  # N*J 
      AX <- temp2$Ay
      BX <- temp2$By
      dAyH0d <- temp2$dAyH0d
      
    }

    WBXp <- W%*%t(BX)    # N*N
    Sigma_X <- qr.solve((W%*%t(BX)),1e-3*Sigma_cP)%*%solve(BX%*%t(W))
                              
  }
  #=============================================================================
  # Now with Sigma_X in hand, compute loadings for AX
  if (isTypicalDiagonal==TRUE){
    
    temp1 <- gaussianDiscreteYieldLoadingsDiagonal(mats_periods, K0Q_X, diag(K1Q_X), Sigma_X, rho0d*dt, rho1d*dt, dt)
    AX <- temp1$Ay
    BX <- temp1$By
    dAyH0d <- temp1$dAyH0d
    
  } else {
    
    temp2 <- gaussianDiscreteYieldLoadingsRecurrence(mats_periods, K0Q_X, K1Q_X, Sigma_X, rho0d*dt, rho1d*dt, dt)
    AX <- temp2$Ay
    BX <- temp2$By
    dAyH0d <- temp2$dAyH0d
    
  }
  #=============================================================================
    
  #=============================================================================
  # Finally, rotate the model to obtain the AcP, BcP loadings.
  # (See above for calculation)
  BcP <- qr.solve(t(W%*%t(BX)), BX)
  AcP <- AX%*%t(diag(J) - t(BX)%*%qr.solve(W%*%t(BX),W))     # 1*J
  #=============================================================================
  
  return(list(BcP=BcP,AcP=AcP,AX=AX,BX=BX,Sigma_X=Sigma_X))
  
}


#=====================================================================
#jszAdjustK1QX(K1Q_X, eps0=1e-5)
#=====================================================================
jszAdjustK1QX <- function(K1Q_X, eps0=1e-5){
  
  # function [K1Q_X, isTypicalDiagonal] = jszAdjustK1QX(K1Q_X, eps0);
  #
  # This function adjusts diagonal K1Q_X to give a non-diagonal but more
  # computationally tractable K1Q_X.
  #
  #
  # K1Q_X can fall into a few cases:
  #   0. diagonal
  #   1. not diagonal
  #   2. zero eigenvalue
  #   3. near repeated roots
  # In cases 1-3, the diagonal closed form solver doesn't work, so compute differently.
  #
  # For case 2/3, we use the method in Joslin, Le, Singleton.  
  #
  # We order the diagonal of diagonal K1Q.
  #
  #
  
  diag_K1Q_X <- diag(K1Q_X)
  isDiagonal <- all(K1Q_X==diag(diag_K1Q_X))
  
  if (isDiagonal == TRUE){
    
    diag_K1Q_X <- sort(diag_K1Q_X)
    K1Q_X <- diag(diag_K1Q_X)
    
    hasNearUnitRoot <- !all(abs(diag_K1Q_X)>eps0)                                 # Only applicable for diagonal
    hasNearRepeatedRoot <- !all(abs(diff(diag_K1Q_X)) > eps0)                     # Only applicable for diagonal
    isTypicalDiagonal <- all(isDiagonal,!hasNearUnitRoot, !hasNearRepeatedRoot)
    
  } else {
    
    isTypicalDiagonal <- FALSE
    
  }
  
  
  if (isDiagonal==TRUE & isTypicalDiagonal==FALSE) {
    
    inds <- (abs(diff(diag_K1Q_X)) < eps0 | abs(diag_K1Q_X[1:length(diag_K1Q_X)-1]) <= eps0)
    if (abs(diag_K1Q_X[length(diag_K1Q_X)]) <= eps0){
      inds[length(inds)] <- TRUE
    }
    
    K1Q_X[1:nrow(K1Q_X)-1,2:ncol(K1Q_X)] = K1Q_X[1:nrow(K1Q_X)-1,2:ncol(K1Q_X)] + diag(inds)
  
  }

  return(list(K1Q_X=K1Q_X[c(3,2,1),c(3,2,1)], isTypicalDiagonal=isTypicalDiagonal))
}

#=====================================================================
#jszRotation(W, K1Q_X, rinfQ, dt, Sigma_cP, mats, BX, AX)
#=====================================================================
jszRotation <- function(W, K1Q_X, rinfQ, dt, Sigma_cP=NULL, mats=NULL, BX, AX){
  
  #function [K0Q_cP, K1Q_cP, rho0_cP, rho1_cP] = jszRotation(W, K1Q_X, rinfQ, dt, Sigma_cP, mats, BX, AX)
  #
  # Either provide (Sigma_cP, mats) or (BX, AX)           
  #
  # Inputs:
  #   W          : N*J,      vector of portfolio weights to fit without error.
  #   K1Q_X      : N*N
  #   rinfQ      : scalar,   the long run mean under Q of the annualized short rate
  #   dt         : scalar,   length of period in years
  #   Sigma_cP   : N*N  covariance of innovations
  #   mats       : 1*J,      maturities in years
  #   BX         : N*J  (BX, AX) are optional (saves time)
  #   AX         : 1*J
  #
  # Returns:
  #   K0Q_cP : N*1
  #   K1Q_cP : N*N
  #   rho0_cP : scalar
  #   rho1_cP : N*1
  #
  #
  # r(t) = rho0_cP + rho1_cP'*cPt
  #      = rinfQ + 1'*Xt
  #      = 1 period discount rate (annualized)
  #
  # Under Q:
  #   X(t+1) - X(t)   =          K1Q_X*X(t)  + eps_X(t+1),   cov(eps_X(t+1)) = Sigma_X
  #   cP(t+1) - cP(t) = K0Q_cP + K1Q_cP*X(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Sigma_cP
  #
  # Where Sigma_X is chosen to match Sigma_cP 
  #
  # cPt = W*yt  (cPt N*1, W is N*J)
  #     = W*AX' + W*BX'*Xt
  #     = WAXp + WBXp*Xt
  #
  # Delta(cP) = WBXp*Delta(Xt)
  #           = WBXp*(K1Q_X*Xt + sqrt(Sigma_X)*eps(t+1))
  #           = WBXp*(K1Q_X)*(WBXp\(cPt - WAXp)) + sqrt(Sigma_cP)*eps(t+1)
  #           = WBXp*(K1Q_X)/WBXp*cPt - WBXp*(K1Q_X)/WBXp*WAXp] + sqrt(Sigma_cP)*eps(t+1)
  #
  # rt = rinfQ + 1'*Xt  [annualized 1-period rate]
  #    = rinfQ + 1'*(WBXp\(cPt - WAXp))
  #    = [rinfQ - 1'*(WBXp\WAXp)] + ((WBXp)'1)'*cPt

  if (is.null(Sigma_cP) == FALSE){
    
    # Adjust K1Q_X in case we have near repeated root/zero eigenvalue
    K1Q_X = jszAdjustK1QX(K1Q_X);
    temp <- jszLoadings(W, K1Q_X, rinfQ, Sigma_cP, mats, dt)
    BcP <- temp$BcP
    AcP <- temp$AcP
    AX <- temp$AX
    BX <- temp$BX
    
  }

  N <- dim(K1Q_X)[1]
  WBXp <- W%*%t(BX)
  WAXp <- W%*%t(AX)
  
  K1Q_cP <- WBXp%*%K1Q_X%*%solve(WBXp)
  K0Q_cP <- -K1Q_cP%*%WAXp
  
  rho0_cP <- rinfQ - matrix(1,1,N)%*%qr.solve(WBXp,WAXp)
  rho1_cP <- qr.solve(t(WBXp),matrix(1,N,1))

  return(list(K0Q_cP=K0Q_cP, K1Q_cP=K1Q_cP, rho0_cP=rho0_cP, rho1_cP=rho1_cP))
}




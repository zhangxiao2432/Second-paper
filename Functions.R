#=====================================================================
#OLS(X, Y)
#=====================================================================
ols <- function(x,y){

  bhat <- solve(t(x)%*%x)%*%t(x)%*%y
  n <- dim(x)[1]
  p <- dim(x)[2]
  res <- y - (x%*%bhat)
  sigs <- sum(res^2)/(n-p)
  vcov <- sigs*solve(t(x)%*%x)
  vcov <- sqrt(diag(vcov))
  
  return(list(bhat=bhat, vcov=vcov))
}

#=====================================================================
#get_exp(X, delta0, delta1, mu, phi)
#=====================================================================
get_exp <- function(X, delta0, delta1, mu, phi){
  
  # Calculate expectations component (average expected short-term interest
  # rate) of long-term interest rate
  # Bauer, Rudebusch, Wu (2014, AER) Comment on ``Term Premia and
  # Inflation Uncertainty: Empirical Evidence from an International Panel Dataset''
  #
  # Arguments:
  #  X - TxN matrix with risk factors
  #  delta0, delta1 - parameters of short-rate equation
  #  mu, phi - VAR parameters
  #
  # Returns:
  #  avexp - expectations component of forward rate from hmin to hmax
  #  quarters
  
  hmin <- 21
  hmax <- 40
  T <- dim(X)[1]
  N <- dim(X)[2]
  avexp <- matrix(0,T,1)
  
  for (t in 1:T) {
    
    g <- matrix(0,N,hmax)
    g[,1] <- X[t,]
    
    for (j in 2:hmax) {
      
      g[,j] <- mu + phi%*%g[,j-1]
      
    }
    g <- g[1:3, hmin:hmax]
  
    temp <- delta0+t(delta1)%*%g
    temp[which(temp< 0)]  <- 0
    avexp[t] <- mean(temp)
  }
  
  return(avexp)
}

#=====================================================================
#irf_var1(Phi, maxlag)
#=====================================================================
irf_var1 <- function(Phi, maxlag=500){
  
  # calculate impulse response function for a VAR(1)
  # for the first variable in response to shocks to the first variable
  if(length(Phi) > 1){
    
    irfvec <- matrix(0, maxlag, 1)
    Psi <- Phi
    irfvec[1] <- Phi[1,1]
    
    for (i in 2:maxlag){
      Psi <- Phi %*% Psi
      irfvec[i] <- Psi[1,1]
    }
  } else {stop('scalar')}

  return(irfvec)
}

#=====================================================================
#Plotting figure 0
#=====================================================================
fig0 <- function(data,k, color){
  
  h0 <- ggplot(NULL,aes(x=data[,1])) + ggtitle(cvec[k])+
    geom_line(aes(y=data[,5], col="1 Yr"),size=1,linetype="solid")+
    geom_line(aes(y=data[,7], col="2 Yr"),size=1,linetype="solid")+
    geom_line(aes(y=data[,9], col="4 Yr"),size=1,linetype="solid")+
    geom_line(aes(y=data[,11], col="6 Yr"),size=1,linetype="solid")+
    geom_line(aes(y=data[,13], col="8 Yr"),size=1,linetype="solid")+ 
    geom_line(aes(y=data[,15], col="10 Yr"),size=1,linetype="solid")+
    labs(x="", y="")+
    scale_color_manual(values=c("1 Yr"=color[1], "2 Yr"=color[2],"4 Yr"=color[3],"6 Yr"=color[4],"8 Yr"=color[5],"10 Yr"=color[6]))+
    theme_bw()+theme(legend.position = "none", plot.title = element_text(hjust=0.5))
  
  return(h0)
}


#=====================================================================
#Plotting figure 1
#=====================================================================
fig1 <- function(data,k){
  
    h1 <- ggplot(NULL,aes(x=data[,1])) + ggtitle(cvec[k])+
      geom_line(aes(y=data[,2], col="Data-generated"),size=1,linetype="solid")+ 
      geom_line(aes(y=data[,3], col="Long-bond"),size=1,linetype="longdash")+ 
      geom_line(aes(y=data[,4], col="Risk-neutral"),size=1,linetype="twodash")+
      geom_ribbon(aes(ymin=data[,5],ymax=data[,6]),alpha=0.2)+
      scale_y_continuous(limits = c(0,12))+
      labs(x="", y="Percent")+
      scale_color_manual(values=c("Data-generated"="#00AFBB", "Long-bond"="#f8766d","Risk-neutral"="#E7B800"))+
      theme_bw()+theme(legend.position = "none", plot.title = element_text(hjust=0.5))
    
    return(h1)
}

#=====================================================================
#Plotting figure 2
#=====================================================================
fig2 <- function(data,startRec,endRec,k){
  
  if(k != 4){
    
    recession <- list()
    recession$start <- matrix(startRec,length(startRec),1) 
    recession$end <- matrix(endRec,length(endRec),1) 
    recession <- as.data.frame(recession)
    
    h2 <- ggplot(NULL,aes(x=data[,1])) + ggtitle(cvec[k])+
      geom_line(aes(y=data[,2], col="OLS"),size=1,linetype="solid")+ 
      geom_line(aes(y=data[,3], col="BC"),size=1,linetype="twodash")+
      scale_y_continuous(limits = c(-2,6))+
      labs(x="", y="Percent")+
      scale_color_manual(values=c("OLS"="#f8766d","BC"="#00AFBB"))+
      theme_bw()+theme(legend.position = "none",plot.title = element_text(hjust=0.5))+
      geom_rect(data= recession, inherit.aes = FALSE,
                aes(xmin=start, xmax=end, ymin=-Inf, ymax=+Inf),fill='pink', alpha=0.3)    
    
  } else {
    
    h2 <- ggplot(NULL,aes(x=data[,1])) + ggtitle(cvec[k])+
      geom_line(aes(y=data[,2], col="OLS"),size=1,linetype="solid")+ 
      geom_line(aes(y=data[,3], col="BC"),size=1,linetype="twodash")+
      scale_y_continuous(limits = c(-2,6))+
      labs(x="", y="Percent")+
      scale_color_manual(values=c("OLS"="#f8766d","BC"="#00AFBB"))+
      theme_bw()+theme(legend.position = "none",plot.title = element_text(hjust=0.5))
    
  }
  
  return(h2)
}

#=====================================================================
# varest(z,p,isintercept)
#=====================================================================
varest <- function(z,p,isintercept){
  
  bigt <- dim(z)[1]
  n <- dim(z)[2]
  y <- z[(p+1):bigt,]
  x <- matrix(1,bigt-p,1)
  
  for(j in 1:p){
    x <- cbind(x,z[(p+1-j):(bigt-j),])
  }
  
  if(isintercept==FALSE){
    x <- x[,2:ncol(x)]
  }
  
  bhat <- t(y)%*%x%*%solve(t(x)%*%x)
  res <- y - (x%*%t(bhat))
  sigma <- (t(res)%*%res)/(bigt-(n*p)-1)
  vcov1 <- kron(solve(t(x)%*%x), sigma)
  d <- duplication(n)
  dplus <- solve(t(d)%*%d)%*%t(d)
  vcov2 <- 2*dplus%*%kron(sigma, sigma)%*%t(dplus)
  vcov1full <- vcov1
  vcov1 <- sqrt(diag(vcov1))
  vcov2 <- sqrt(diag(vcov2))
  
  return(list(bhat=bhat,sigma=sigma,vcov1=vcov1,vcov2=vcov2,vcov1full=vcov1full))
}

#=====================================================================
# duplication(n)
#=====================================================================
duplication <- function(n){
  
  # duplication(n)
  # Returns Magnus and Neudecker's duplication matrix of size n
  # Author: Thomas P Minka (tpminka@media.mit.edu)

  a <- matrix(0,n,n)
  k <- 1
  
  for(j in 1:n){
    for(i in 1:n){
      if(i >= j){
        a[i,j] <- k
        k <- k+1
      } else {
        a[i,j] <- a[j,i]
      }
    }
  }
  
  j <- matrix(a,(dim(a)[1])*(dim(a)[2]),1)
  m <- n*(n+1)/2
  d <- matrix(0,n*n,m)
  
  for(r in 1:dim(d)[1]){
    d[r, j[r]] <- 1
  }
  
  return(d)
}

#=====================================================================
# jpslikel(param,m,W,mats,Sigma_cP_OLS)
#=====================================================================
jpslikel <- function(param,m,W,mats,Sigma_cP_OLS){

  # Value of likelihood function for DTSM with unspanned macro risks
  # Bauer, Rudebusch, Wu (2014, AER) Comment on ``Term Premia and
  # Inflation Uncertainty: Empirical Evidence from an International Panel Dataset'' 
  #
  # This code was generously provided to us by Jonathan Wright.
  
  lb <- matrix(c(-0.6, -0.6, -0.6, 0, -Inf*matrix(1,15,1)), 19, 1) 
  ub <- matrix(c(-.0001, -.0001, -.0001, 0.4, Inf*matrix(1,15,1)), 19, 1) 
  A <- cbind(matrix(c(-1, 1, 0, 0, -1, 1),2,3,byrow=TRUE), matrix(0,2,16))
  b <- matrix(c(0,0),2,1)
  
  if(all(param<ub)==FALSE | all(param>lb)==FALSE | all(A%*%param<b)==TRUE){
    likel <- 100000
  } else {
    
    nf <- dim(Sigma_cP_OLS)[1]
    TT <- dim(m[2:nrow(m),])[1]
    J <- dim(m[2:nrow(m),])[2]
    N <- dim(W)[1]
    ch <- param[5:length(param)]
    
    ch.temp <- matrix(0,nf,nf)
    ch.temp[lower.tri(ch.temp,diag = TRUE)] <- ch
    
    ch <- ch.temp
    Sigma_cP <- ch%*%t(ch)
    K1Q_X <- diag(param[1:3])
    rinfQ <- param[4]
    cP <- m%*%t(W)                  # (T+1)*N, cP stands for math caligraphic P.
    
    temp1 <- jszLoadings(W,K1Q_X,rinfQ,Sigma_cP[1:3,1:3],mats,1/4)
    BcP <- temp1$BcP
    AcP <- temp1$AcP
    AX <- temp1$AX
    BX <- temp1$BX
    
    temp2 <- jszRotation(W,K1Q_X,rinfQ,1/4,NULL,NULL,BX,AX)
    K0Q_cP <- temp2$K0Q_cP
    K1Q_cP <- temp2$K1Q_cP
    rho0_cP <- temp2$rho0_cP
    rho1_cP <- temp2$rho1_cP
    
    yields_m <- matrix(1,TT+1,1)%*%AcP + cP%*%BcP                    # (T+1)*J, model-implied yields
    yield_errors <- m[2:nrow(m),] - yields_m[2:nrow(yields_m),]      # T*J
    square_orthogonal_yield_errors <- (yield_errors)^2             # T*J, but N-dimensional projection onto W is always 0, so effectively (J-N) dimensional
    sigma_e <- sqrt(sum(square_orthogonal_yield_errors)/(TT*(J-N)))
    llkQ <- 0.5*sum(apply(square_orthogonal_yield_errors,1,sum)/sigma_e^2 + (J-N)*log(2*pi) + (J-N)*log(sigma_e^2))
    llkP <- (0.5*N*TT*log(2*pi))+(0.5*TT*log(det(Sigma_cP)))+(0.5*TT*sum(diag(ginv(Sigma_cP)%*%Sigma_cP_OLS)))
    likel <- llkQ + llkP
    
  }
  
  return(likel)
}


#=====================================================================
# Helper function to add titles
#=====================================================================
# - sheet : sheet object to contain the title
# - rowIndex : numeric value indicating the row to 
#contain the title
# - title : the text to use as title
# - titleStyle : style object to use for title
xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}

#=====================================================================
# Long-Term recovery
#=====================================================================
LongTermRec <- function(m, W, K1Q_X, rinfQ, Sigma_cP, mats, dt, Sigma_X=NULL){

  #   K1Q_X      : N*N
  #   rinfQ      : scalar,   the long run mean under Q of the annualized short rate
  #   Sigma_cP, Sigma_X : N*N  covariance of innovations. PROVIDE ONE OR THE OTHER
  #
  # Returns:
  #   AX  : 1*1
  #   BX  : N*1
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
 
  while(mats_periods[length(mats_periods)] < 1000){
    
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

    WAXp <- W%*%t(AX)    # N*1
    WBXp <- W%*%t(BX)    # N*N
      
    # Rotate the model to obtain the AcP, BcP loadings.
    # (See above for calculation)
  
    BcPL <- ginv(BX%*%t(W))%*%BX
    AcPL <- AX%*%t(diag(J) - t(BX)%*%ginv(W%*%t(BX))%*%W)     # 1*J
    #=============================================================================
    
    # risk factors
    cP=m%*%t(W)
    
    # fitted yield of long bond
    fitted <- (matrix(1,bigt,1)%*%AcPL)+(cP%*%BcPL)
    
    # fitted bond price
    P.t0 <- exp(fitted[,ncol(fitted)-1])
    P.t1 <- exp(fitted[,ncol(fitted)])
    P_ratio <- P.t1/P.t0
    Mn <- max(P_ratio)
    mn <- min(P_ratio)

    if(abs(Mn - mn) < 1e-3){
      
      cat("Long term rate converged! \n")
      break
      
    } else {
      
      mats_periods <- mats_periods + 4
      
    }
  }
  
  NT <- round(mats_periods[1]/4)
  lamda.long <- -log((Mn+mn)/2)
  eigen.fun <- exp(lamda.long*NT)*P.t0
  yield.long <- fitted[,ncol(fitted)]

  
  return(list(NT=NT, AcPL=AcPL, BcPL=BcPL, eigen.fun=eigen.fun, yield.long=yield.long))

}

#=====================================================================
# Estimate 10 countries lamda
#=====================================================================
agg.countries.lamda <- function(m10, m_m10, ns){
  
  # Principle Component Analysis
  PCA <- princomp(m10)
  score <- PCA$scores
  score <- score[,1:3]
  score <- cbind(-score[,1:2],score[,3]) # keep consistent with matlab
  bigt <- dim(score)[1]
  coeff <- as.matrix(PCA$loadings)
  W <- t(coeff[,1:3])
  mats <- ns/4
  
  param_est <- DTSM_est(score, W, m10, m_m10, NULL, ns, 1)
  
  sigma <- param_est$sigma
  phi <- param_est$phi
  mu <- param_est$mu
  param <- param_est$param
  ch <- param_est$ch
  K1Q_X <- param_est$K1Q_X
  rinfQ <- param_est$rinfQ
  
  LB <- LongTermRec(m10, W, K1Q_X, rinfQ, ch%*%t(ch), mats, 1/4)
  BcPL <- LB$BcPL
  AcPL <- LB$AcPL
 
  return((list(BcPL=BcPL, AcPL=AcPL)))
}
  

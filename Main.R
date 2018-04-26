rm(list = ls())

library("R.matlab")
library("readxl")
library("amap")
library("expm")
library("MASS")
library("pracma")
library("ggplot2")
library("grid")
library("gridExtra")
library("xtable")
library("xlsx")

source("Functions.R")
source("JSZ.R")
source("BC.R")

set.seed(1030)

est_vec <- c('OLS', 'BC')

M <- 2000  # bootstrap sample
N <- 5

#Import zero coupon yield data(10 countries) from matlab
data1 <- readMat("10 countries/zcbis.mat")
#Australia
zcaus <- as.matrix(as.data.frame(data1[1]))
#Canada
zcca <- as.matrix(as.data.frame(data1[2]))
#Switzerland
zcch <- as.matrix(as.data.frame(data1[3]))
#Germany
zcde <- as.matrix(as.data.frame(data1[4]))
#Japan
zcjp <- as.matrix(as.data.frame(data1[5]))
#Norway
zeno <- as.matrix(as.data.frame(data1[6]))
#New Zealand
zenz <- as.matrix(as.data.frame(data1[7]))
#Sweden
zcse <- as.matrix(as.data.frame(data1[8]))
#U.K.
zcuk <- as.matrix(as.data.frame(data1[9]))
#U.S.
zcus <- as.matrix(as.data.frame(data1[10]))

remove(data1)

#Import macro data(10 countries) from matlab
data2 <- readMat("10 countries/macrovars.mat")
#Australia
aus <- as.matrix(as.data.frame(data2[c(1,11,23,39,49)]))
aus[is.nan(aus)==TRUE] <- NA
#Canada
ca <- as.matrix(as.data.frame(data2[c(2,12,24,35,45)]))
ca[is.nan(ca)==TRUE] <- NA
#Switzerland
chh <- as.matrix(as.data.frame(data2[c(8,18,25,38,48)]))
chh[is.nan(chh)==TRUE] <- NA
#Germany
de <- as.matrix(as.data.frame(data2[c(3,13,26,33,43)]))
de[is.nan(de)==TRUE] <- NA
#Japan
jp <- as.matrix(as.data.frame(data2[c(4,14,27,32,42)]))
jp[is.nan(jp)==TRUE] <- NA
#Norway
no <- as.matrix(as.data.frame(data2[c(6,16,28,36,46)]))
no[is.nan(no)==TRUE] <- NA
#New Zealand
nz <- as.matrix(as.data.frame(data2[c(5,15,30,40,50)]))
nz[is.nan(nz)==TRUE] <- NA
#Sweden
se <- as.matrix(as.data.frame(data2[c(7,17,29,37,47)]))
se[is.nan(se)==TRUE] <- NA
#U.K.
uk <- as.matrix(as.data.frame(data2[c(9,19,22,34,44)]))
uk[is.nan(uk)==TRUE] <- NA
#U.S.
us <- as.matrix(as.data.frame(data2[c(10,20,21,31,41)]))
us[is.nan(us)==TRUE] <- NA
remove(data2)

cvec <- c('Canada','Switzerland','Germany','Norway','Sweden','United Kingdom.','United States.','Australia','New Zealand','Japan')
ki <- c(7, 10, 3, 6, 1, 4, 5, 2, 8, 9)

dtplot <- as.matrix(zcus[219:nrow(zcus),1] + zcus[219:nrow(zcus),2]/12)
dtplot <- as.matrix(dtplot[seq(3,nrow(dtplot),3)])

# Read in ECRI dates
ecri <- read_excel("10 countries/Ecridates.xls", sheet=1)
ecri <- as.data.frame(ecri[362:594,c(1,10,9,4,5,3,6,8,2,7)])
ecri <- ecri[seq(2,nrow(ecri),3),]
ecri <- ecri[-nrow(ecri),]
ecri <- as.matrix(cbind(ecri[,2:6],0,ecri[,7:10]))
Time <- NROW(ecri[,1])

# difference between early and late part of the sample
d_forw <- matrix(0,10,1)                                      #hard coding
d_fhat <- matrix(0,10,2)
d_exp <- matrix(0,10,2)
d_ftp <- matrix(0,10,2)
ind_early <- c(3:7)
ind_late <- c(73:77)

figure0 <- figure1 <- figure2 <- list()
param_est.ols <- param_est.bc2 <- param_est.bc3 <- lamda.long <- list()
table1 <- table2 <- term_premium.ols <- term_premium.bc2 <- m10.panel <- m_m10.panel <- c()
color <- c("black","red","green3","blue","cyan","magenta","yellow")

for(i in 1:10){

  k <- ki[i]
  if(k==1){
    m <- zcca
    m_m <- ca[,1:2]
  } else if(k==2){
    m <- zcch
    m_m <- chh[,1:2]
  } else if(k==3){
    m <- zcde
    m_m <- de[,1:2]
  } else if(k==4){
    m <- zeno
    m_m <- no[,1:2]
  } else if(k==5){
    m <- zcse
    m_m <- se[,1:2]
  } else if(k==6){
    m <- zcuk
    m_m <- uk[,1:2]
  } else if(k==7){
    m <- zcus
    m_m <- us[,1:2]
  } else if(k==8){
    m <- zcaus
    m_m <- aus[,1:2]
  } else if(k==9){
    m <- zenz
    m_m <- nz[,1:2]
  } else {
    m <- zcjp
    m_m <- jp[,1:2]
  }
  
  if(nrow(m)>233){m <- tail.matrix(m,233)}
  forw <- as.matrix(2*m[,42]-m[,22])
  
  # Go to quarterly data
  j <- dim(m)[1]
  m <- m[seq(3,j,3),]
  forw <- forw[seq(3,j,round((3*((j/3)-floor(j/3)))+1)),]
  m_m <- tail.matrix(m_m,dim(m)[1])
  
  ns <- c(1,2,3,4,6,8,12,16,20,24,28,32,36,40)
  dt <- m[,1:2]
  m <- 0.01*m[,ns+2]
  colnames(m) <- NULL
  
  # Plot data
  data.temp0 <- cbind(dtplot[1:length(m[,1])],m)
  h0 <- fig0(data.temp0,k,color)
  print(h0)
  figure0[[i]] <- h0
  
  # Principle Component Analysis
  PCA <- princomp(m)
  score <- PCA$scores
  score <- score[,1:3]
  score <- cbind(-score[,1:2],score[,3]) # keep consistent with matlab
  bigt <- dim(score)[1]
  
  coeff <- as.matrix(PCA$loadings)
  
  fhat <- frn <- ftp <- fhat.LB <- matrix(0,bigt,2)
  rmse <- maxeig <- maxeigQ <- irf_5y <- irfQ_5y <- vol_f <- vol_frn <- vol_ftp <- matrix(NA,1,2)
  
  #==================================================================================
  # Estimate parameters using 1. OLS 
  #                           2. BC2 -- indirect inference bias correction
  #                           3. BC3 -- bootstrap bias correction
  #==================================================================================
  cat("\n")
  cat("=======================", cvec[k], "======================= \n")
  cat("\n")
  
  # for(specId in 1:3){
  #   
  #   if(specId==1){
  #     param_est.ols[[i]] <- DTSM_est(score, coeff, m_m, NULL, ns, specId)
  #     initial <- param_est.ols[[i]]$param
  #   } 
  #   else if(specId==2){
  #     param_est.bc2[[i]] <- DTSM_est(score, coeff, m_m, initial, ns, specId)
  #   }
  #   else {
  #     param_est.bc3[[i]] <- DTSM_est(score, coeff, m_m, initial, ns, specId)
  #   }
  #   
  # }
  #==================================================================================
  
  for (iEst in 1:2) {
    
# Read matlab data

    param_ols <- list.files("10 countries/param_ols", pattern = ".mat")
    pathname_ols <- file.path("10 countries/param_ols",param_ols)

    param_bc2 <- list.files("10 countries/param_bc2", pattern = ".mat")
    pathname_bc2 <- file.path("10 countries/param_bc2",param_bc2)

    if(iEst==1){

      param.ols <- readMat(pathname_ols[i])
      sigma <- param.ols$sigma
      phi <- param.ols$phi
      mu <- param.ols$mu
      param <- param.ols$param
      ch <- param.ols$ch
      K1Q_X <- param.ols$K1Q.X
      rinfQ <- param.ols$rinfQ

    } else {

      param.bc2 <- readMat(pathname_bc2[i])
      sigma <- param.bc2$sigma
      phi <- param.bc2$phi
      mu <- param.bc2$mu
      param <- param.bc2$param
      ch <- param.bc2$ch
      K1Q_X <- param.bc2$K1Q.X
      rinfQ <- param.bc2$rinfQ

    }

    
    # 
    # if(iEst==1){
    #   
    #   param.ols <- param_est.ols[[i]]
    #   sigma <- param.ols$sigma
    #   phi <- param.ols$phi
    #   mu <- param.ols$mu
    #   param <- param.ols$param
    #   ch <- param.ols$ch
    #   K1Q_X <- param.ols$K1Q_X
    #   rinfQ <- param.ols$rinfQ
    #   
    # } else {
    #   
    #   param.bc2 <- param_est.bc2[[i]]
    #   sigma <- param.bc2$sigma
    #   phi <- param.bc2$phi
    #   mu <- param.bc2$mu
    #   param <- param.bc2$param
    #   ch <- param.bc2$ch
    #   K1Q_X <- param.bc2$K1Q_X
    #   rinfQ <- param.bc2$rinfQ
    #   
    # }

    x <- cbind(1,score)
    delta <- as.matrix(ols(x,m[,1]/4)$bhat)
    vcov3 <- as.matrix(ols(x,m[,1]/4)$vcov)
    delta0 <- delta[1]
    delta1 <- c(delta[2],-delta[3],delta[4])  #different sign with Matlab
    
    W <- t(coeff[,1:3])
    mats <- ns/4
    jsz.temp1 <- jszLoadings(W,K1Q_X,rinfQ,ch%*%t(ch),mats,1/4)
    
    BcP <- jsz.temp1$BcP
    AcP <- jsz.temp1$AcP
    AX <- jsz.temp1$AX
    BX <- jsz.temp1$BX
    
    jsz.temp2 <- jszRotation(W,K1Q_X,rinfQ,1/4,NULL,NULL,BX,AX)
      
    K0Q_cP <- jsz.temp2$K0Q_cP
    K1Q_cP <- jsz.temp2$K1Q_cP
    rho0_cP <- jsz.temp2$rho0_cP
    rho1_cP <- jsz.temp2$rho1_cP
    phiQ <- rbind(cbind(K1Q_cP + diag(1,3,3),matrix(0,3,2)), matrix(0,2,5))
    
    X <- cbind(score, m_m)
    
    # risk factors
    cP=m%*%t(W)
    
    # fitted yields
    mfit <- (matrix(1,bigt,1)%*%AcP)+(cP%*%BcP)
    
    # Long term bond Yield
    LB <- LongTermRec(m, W, K1Q_X, rinfQ, ch%*%t(ch), mats, 1/4)
    fhat.LB[,iEst] <- 100*LB$yield.long 
    AcPL <- LB$AcPL
    BcPL <- LB$BcPL
    lamda.long[[i]] <- list(AcPL, BcPL)
    
    # fitted forward rates
    fhat[,iEst] <- 100*(2*mfit[, ncol(mfit)] - mfit[, ncol(mfit)-5])
    frn[,iEst] <- 400*get_exp(X, delta0, delta1, mu, phi)                           #a little bit different with matlab
    ftp[,iEst] <- fhat[,iEst] - frn[,iEst]
    
    if(iEst==2){
      
      # Aggregate 10 countries
      m10 <- cbind(i,1:nrow(m),m)
      m10.panel <- rbind(m10.panel,m10)

      m_m10 <- cbind(i,1:nrow(m_m),m_m)
      m_m10.panel <- rbind(m_m10.panel,m_m10)
      
      # confidence intervals for risk-neutral rates
      phi_ols <- est_VAR(X,1,TRUE,FALSE)$Gamma_hat
      frn_smpl <- matrix(0,bigt, M)
      ftp_smpl <- matrix(0,bigt, M)
      X_sim <- genVAR(phi, 2000, X, 1)    # bootstrap sample
      
      for(b in 1:M){
        
        # estimate VAR on each one
        phi_ols_i <- est_VAR(t(X_sim[,,b]),1,TRUE,FALSE)$Gamma_hat
        # bias correction
        phi_i <- phi_ols_i - (phi_ols - phi)                    #need to review
        # stationarity adjustment
        phi_i <- shrink_phi(phi_i, phi_ols_i, 1, 0)
        mu_i <- (diag(N) - phi_i) %*% as.matrix(colMeans(X))
        frn_smpl[,b] <- 400*get_exp(X, delta0, delta1, mu_i, phi_i)
        ftp_smpl[,b] <- fhat[,iEst] - frn_smpl[,b]
      }
      
      frn_lb <- matrix(0,bigt,1)
      frn_ub <- matrix(0,bigt,1)
      ftp_lb <- matrix(0,bigt,1)
      ftp_ub <- matrix(0,bigt,1)
      
      for(t in 1:bigt){
        
        frn_lb[t] <- quantile(frn_smpl[t,], 0.05)
        frn_ub[t] <- quantile(frn_smpl[t,], 0.95)
        ftp_lb[t] <- quantile(ftp_smpl[t,], 0.05)
        ftp_ub[t] <- quantile(ftp_smpl[t,], 0.95)
        
      }
    }
    
    rmse[iEst] <- format(sqrt(mean((100*mfit-100*m)^2)), digits=3)
    maxeig[iEst] <- format(max(abs(eigen(phi)$values)), digits=4)
    maxeigQ[iEst] <- format(max(abs(eigen(phiQ[1:3,1:3])$values)), digits=4)
    
    hmax <- 200
    h <- 20
    irf <- irf_var1(phi, hmax)
    irf_5y[iEst] <- round(irf[h], digits=3)
    
    irfQ <- irf_var1(phiQ, hmax)
    irfQ_5y[iEst] <- round(irfQ[h], digits=3)
    
    vol_f[iEst] <- round(sd(fhat[,iEst]), digits=3)
    vol_frn[iEst] <- round(sd(frn[,iEst]), digits=3) 
    vol_ftp[iEst] <- round(sd(ftp[,iEst]), digits=3)
     
  }
  
  if(k==4|k==5){  #Norway and Sweden
    
    fhat <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(fhat)), fhat) 
    fhat.LB <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(fhat.LB)), fhat.LB) 
    frn <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(frn)), frn)
    ftp <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(ftp)), ftp)
    forw <- rbind(matrix(NA,nrow(dtplot)-nrow(m),1), matrix(forw,length(forw),1))
    frn_lb <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(frn_lb)), frn_lb)
    frn_ub <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(frn_ub)), frn_ub)
    ftp_lb <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(ftp_lb)), ftp_lb)
    ftp_ub <- rbind(matrix(NA,nrow(dtplot)-nrow(m),ncol(ftp_ub)), ftp_ub)
    
  } 
  
  tp <- fhat - frn        #a little bit different with matlab because frn
  
  # averages over early and late part of the sample
  d_forw[i] <- round(100*(mean(forw[ind_late]) - mean(forw[ind_early])),digits=0)
  d_fhat[i,] <- round(100*(colMeans(fhat[ind_late,], na.rm = TRUE) - colMeans(fhat[ind_early,], na.rm = TRUE)),digits=0)
  d_exp[i,] <- round(100*(colMeans(frn[ind_late,], na.rm = TRUE) - colMeans(frn[ind_early,], na.rm = TRUE)),digits=0)
  d_ftp[i,] <- round(100*(colMeans(tp[ind_late,], na.rm = TRUE) - colMeans(tp[ind_early,], na.rm = TRUE)),digits=0)
  
  # obtain start and end dates for recessions
  j <- 0                # recession
  inRec <- FALSE        # currently in recession
  startRec <- endRec <- c()
  
  for(t in 1:Time){
    
    if(ecri[t,i]==1 & inRec==FALSE){
      # new recession
      j <- j + 1
      inRec <- TRUE
      startRec <- cbind(startRec,dtplot[t])
    } 
    else if(ecri[t,i]==0 & inRec==TRUE){
      # end of recession
      inRec <- FALSE
      endRec <- cbind(endRec,dtplot[t-1])
    }
  }
  
  if(inRec==TRUE){
    endRec[j] <- max(dtplot)
  }
  
  nRec <- j      # number of recessions
  
  # save term premium ====================================================
  if(iEst==1){
    
    tp.ols <- cbind(dtplot,tp)
    tp.ols <- rbind(c(cvec[i], "",""), tp.ols)
    term_premium.ols <- cbind(term_premium.ols, tp.ols)
    
  } else {
    
    tp.bc2 <- cbind(dtplot,tp)
    tp.bc2 <- rbind(c(cvec[i], "",""), tp.bc2)
    term_premium.bc2 <- cbind(term_premium.bc2, tp.bc2)
    
  }

  # plot Figure 1 -- risk-neutral rates ==================================
  # fitted
  data.temp1 <- cbind(dtplot,fhat[,1],fhat.LB[,2],frn[,2],frn_lb,frn_ub)
  h1 <- fig1(data.temp1,k)
  print(h1)
  figure1[[i]] <- h1

  # plot Figure 2 -- term premia =========================================
  data.temp2 <- cbind(dtplot,ftp[,1],ftp[,2])
  h2 <- fig2(data.temp2,startRec,endRec,k)
  
  print(h2)
  figure2[[i]] <- h2
  remove("h1", "h2", "data.temp1", "data.temp2")

  # Table 1 -- summary statistics ========================================
  for(iEst in 1:2){
  table1 <- rbind(table1, c(est_vec[iEst],rmse[iEst],maxeig[iEst],maxeigQ[iEst],irf_5y[iEst],
                                                irfQ_5y[iEst],vol_f[iEst],vol_frn[iEst],vol_ftp[iEst])) 
  }
  
  cat("Completed", i, "/10 \n")
}

# # Reshape the 10 countries panel
# colnames(m10.panel) <- c("countries","T","","","","","","","","","","","","","","")
# rownames(m10.panel) <- NULL
# temp1 <- as.data.frame(m10.panel)
# temp2 <- temp1[order(temp1$T),]
# temp3 <- as.matrix(temp2[,-c(1,2)])
# 
# colnames(m_m10.panel) <- c("countries","T","","")
# rownames(m_m10.panel) <- NULL
# temp4 <- as.data.frame(m_m10.panel)
# temp5 <- temp4[order(temp4$T),]
# temp6 <- as.matrix(temp5[,-c(1,2)])
# temp7 <- agg.countries.lamda(temp3, temp6, ns)
# 
# agg.AcPL <- temp7$AcPL
# agg.BcPL <- temp7$BcPL
# remove(temp1,temp2,temp3,temp4,temp5,temp6,temp7)


# Figures output =======================================================
pdf("Figure 0.pdf", width=8, height=12)
grid.arrange(figure0[[1]], figure0[[2]]+theme(legend.position = c(0.7,0.75),legend.title = element_blank(),legend.direction = "horizontal",
             legend.background = element_rect(colour = "black",size=0.1),plot.title = element_text(hjust=0.5)), figure0[[3]], figure0[[4]],
             figure0[[5]],figure0[[6]], figure0[[7]], figure0[[8]], figure0[[9]], figure0[[10]], ncol=2)

dev.off()


pdf("Figure 1.pdf", width=8, height=12)
grid.arrange(figure1[[1]], figure1[[2]]+theme(legend.position = c(0.75,0.65),legend.title = element_blank(),legend.direction = "vertical",
       legend.background = element_rect(colour = "black",size=0.2),plot.title = element_text(hjust=0.5)), 
       figure1[[3]], figure1[[4]], figure1[[5]],figure1[[6]], figure1[[7]], 
       figure1[[8]], figure1[[9]], figure1[[10]],ncol=2)


dev.off()
pdf("Figure 2.pdf", width=8, height=12)
grid.arrange(figure2[[1]], figure2[[2]]+theme(legend.position = c(0.8,0.8),legend.title = element_blank(),legend.direction = "horizontal",
            legend.background = element_rect(colour = "black",size=0.3),plot.title = element_text(hjust=0.5)), figure2[[3]], figure2[[4]],
            figure2[[5]], figure2[[6]], figure2[[7]], figure2[[8]], figure2[[9]], figure2[[10]], ncol=2)

dev.off()

# Table 2 -- downward trend in long-term rates =======================================
# -> inflation expectations vs. risk-neutral rates

# Consensus survey: first obs is October 1990, last obs is April 2009
consensus0 <- read_excel("10 countries/Longtermsurveydata.xls", sheet="consensus", range = "C2:AF40", col_names=TRUE)

#US, Japan, Germany, UK, Canada, Norway, Sweden, Switzerland, Australia, New Zealand
pi_exp <- as.data.frame(consensus0[,((1:10)-1)*3+2])
g_exp <- as.data.frame(consensus0[,((1:10)-1)*3+1])

for(i in c(1:5,9)){
  
  k <- ki[i]
  g_diff <- round(100*(mean(g_exp[(nrow(g_exp)-2):nrow(g_exp),i], na.rm=TRUE) - mean(g_exp[1:3,i], na.rm=TRUE)), digits=0)
  pi_diff <- round(100*(mean(pi_exp[(nrow(pi_exp)-2):nrow(pi_exp),i], na.rm=TRUE) - mean(pi_exp[1:3,i], na.rm=TRUE)), digits=0)
  
  # for Australia first obs is NA
  if (i==9){
    
    g_diff <- round(100*(mean(g_exp[(nrow(g_exp)-2):nrow(g_exp),i], na.rm=TRUE) - mean(g_exp[2:3,i], na.rm=TRUE)), digits=0)
    pi_diff <- round(100*(mean(pi_exp[(nrow(pi_exp)-2):nrow(pi_exp),i], na.rm=TRUE) - mean(pi_exp[2:3,i], na.rm=TRUE)), digits=0)
    
  }
  
  for(iEst in 1:2){
    if(iEst==1){
      table2 <- rbind(table2, c(cvec[k],pi_diff,g_diff,d_forw[i],est_vec[iEst],
                                             d_fhat[i,iEst],d_exp[i,iEst],d_ftp[i,iEst]))
    } else {
      table2 <- rbind(table2, c("","","","",est_vec[iEst],d_fhat[i,iEst],d_exp[i,iEst],d_ftp[i,iEst]))
    }
  }
  
}

colnames(table1) <- c(" ", "RMSE", "Max.eigen(P)", "Max.eigen(Q)","IRF(5y)(P)","IRF(5y)(Q)","Vfr(Fitted)","Vfr(Exp)","Vfr(FTP)")
colnames(table2) <- c(" ", "Inflation Exp", "Growth Exp", "FRA"," ","Fitted","Exp","FTP")
xtable(table1)
xtable(table2)

write.xlsx(table1, file = "Output.xlsx", sheetName = "table 1", row.names = FALSE, col.names = TRUE, append = TRUE)
write.xlsx(table2, file = "Output.xlsx", sheetName = "table 2", row.names = FALSE, col.names = TRUE, append = TRUE)
write.xlsx(term_premium.ols, file = "Output.xlsx", sheetName = "term_premium_ols", row.names = FALSE, col.names = FALSE, append = TRUE)
write.xlsx(term_premium.bc2, file = "Output.xlsx", sheetName = "term_premium_bc2", row.names = FALSE, col.names = FALSE, append = TRUE)




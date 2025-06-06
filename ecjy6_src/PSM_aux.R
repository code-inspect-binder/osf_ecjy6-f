#Appendix. Auxiliary R subroutines (PSM_aux.R)
library(compiler)
library(boot)
library(foreach)
library(doMC)

#initialize model-generating functions
SDE_lm <- function(nTransit=0){
  PK <- list()
  if(nTransit == 0){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$ka,0, phi$ks,-phi$ke), nrow=2, ncol=2, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,phi$ka, 0,0), nrow=2, ncol=2, byrow=T)
      matC <- matrix(c(0,phi$scal), nrow=1, ncol=2)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      matrix(c(phi$sigma,0, 0,0), nrow=2, byrow=T)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(1,phi$init), nrow=2)
    }
  }else if(nTransit == 1){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0, phi$kt,-phi$ka,0, 0,phi$ks,-phi$ke), nrow=3, ncol=3, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,phi$ka, 0,0), nrow=3, ncol=2, byrow=T)
      matC <- matrix(c(0,0,phi$scal), nrow=1, ncol=3)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      matrix(c(0,0,0, phi$sigma,0,0, 0,0,0), nrow=3, byrow=T)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,1,phi$init), nrow=3)
    }
  }else if(nTransit == 2){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0,0, phi$kt,-phi$kt,0,0, 0,phi$kt,-phi$ka,0, 0,0,phi$ks,-phi$ke), nrow=4, ncol=4, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,0, 0,phi$ka, 0,0), nrow=4, ncol=2, byrow=T)
      matC <- matrix(c(0,0,0,phi$scal), nrow=1, ncol=4)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      matrix(c(0,0,0,0, 0,0,0,0, phi$sigma,0,0,0, 0,0,0,0), nrow=4, byrow=T)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,0,1,phi$init), nrow=4)
    }
  }else if(nTransit == 3){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0,0,0, phi$kt,-phi$kt,0,0,0, 0,phi$kt,-phi$kt,0,0, 0,0,phi$kt,-phi$ka,0, 0,0,0,phi$ks,-phi$ke), nrow=5, ncol=5, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,0, 0,0, 0,phi$ka, 0,0), nrow=5, ncol=2, byrow=T)
      matC <- matrix(c(0,0,0,0,phi$scal), nrow=1, ncol=5)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      matrix(c(0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, phi$sigma,0,0,0,0, 0,0,0,0,0), nrow=5, byrow=T)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,0,0,1,phi$init), nrow=5)
    }
  }else if(nTransit == 4){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0,0,0,0, phi$kt,-phi$kt,0,0,0,0, 0,phi$kt,-phi$kt,0,0,0, 0,0,phi$kt,-phi$kt,0,0, 0,0,0,phi$kt,-phi$ka,0, 0,0,0,0,phi$ks,-phi$ke), nrow=6, ncol=6, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,0, 0,0, 0,0, 0,phi$ka, 0,0), nrow=6, ncol=2, byrow=T)
      matC <- matrix(c(0,0,0,0,0,phi$scal), nrow=1, ncol=6)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      SIG <- matrix(0, nrow=6, ncol=6)
      SIG[5,1] <- phi$sigma
      return(SIG)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,0,0,0,1,phi$init), nrow=6)
    }
  }else if(nTransit == 5){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0,0,0,0,0, phi$kt,-phi$kt,0,0,0,0,0, 0,phi$kt,-phi$kt,0,0,0,0, 0,0,phi$kt,-phi$kt,0,0,0, 0,0,0,phi$kt,-phi$kt,0,0, 0,0,0,0,phi$kt,-phi$ka,0, 0,0,0,0,0,phi$ks,-phi$ke), nrow=7, ncol=7, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,0, 0,0, 0,0, 0,0, 0,phi$ka, 0,0), nrow=7, ncol=2, byrow=T)
      matC <- matrix(c(0,0,0,0,0,0,phi$scal), nrow=1, ncol=7)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      SIG <- matrix(0, nrow=7, ncol=7)
      SIG[6,1] <- phi$sigma
      return(SIG)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,0,0,0,0,1,phi$init), nrow=7)
    }
  }else if(nTransit == 6){
    PK$Matrices <- function(phi) {
      matA <- matrix(c(-phi$kt,0,0,0,0,0,0,0, phi$kt,-phi$kt,0,0,0,0,0,0, 0,phi$kt,-phi$kt,0,0,0,0,0, 0,0,phi$kt,-phi$kt,0,0,0,0, 0,0,0,phi$kt,-phi$kt,0,0,0, 0,0,0,0,phi$kt,-phi$kt,0,0, 0,0,0,0,0,phi$kt,-phi$ka,0, 0,0,0,0,0,0,phi$ks,-phi$ke), nrow=8, ncol=8, byrow=T)
      diag(matA) <- diag(matA) + 1e-10
      matB <- matrix(c(phi$ka*phi$stress,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,phi$ka, 0,0), nrow=8, ncol=2, byrow=T)
      matC <- matrix(c(0,0,0,0,0,0,0,phi$scal), nrow=1, ncol=8)
      matD <- matrix(c(0,0), nrow=1, ncol=2)
      list(matA=matA, matB=matB, matC=matC, matD=matD)
    }
    PK$SIG = function(phi) { #State variance
      SIG <- matrix(0, nrow=8, ncol=8)
      SIG[7,1] <- phi$sigma
      return(SIG)
    }
    PK$X0 = function(Time=NA, phi, U=NA) {
      matrix(c(0,0,0,0,0,0,1,phi$init), nrow=8)
    }
  }else{stop("nTransit < 0 or nTransit > 6")}
  
  PK$h = function(eta, theta, covar) {
    phi <- theta
    phi$scal <- exp(theta$scal*covar[1])
    phi$stress <- theta$stress*exp(-theta$sex*covar[2])
    phi$init <- (theta$ks/theta$ke)#*exp(theta$init)
    phi
  }
  PK$S = function(phi) { #Error variance
    matrix(phi[["S"]])
  }
  PK$ModelPar = function(THETA){
    list(
      theta=list(stress = THETA['stress'], ka = THETA['ka'], kt = THETA['kt'], ks= THETA['ks'], ke = THETA['ke'], init = THETA['init'], sigma=THETA['sigma'], S=THETA['S'], scal = THETA['scal'], sex = THETA['sex']),
      OMEGA=NULL
    )
  }
  return(PK)
}
ODE_lm <- function(nTransit=0){
  PK <- SDE_lm(nTransit = nTransit)
  PK$SIG <- function(phi) matrix(0, nrow=2+nTransit, ncol=2+nTransit)
  return(PK)
}

SDE_lme <- function(nTransit=0, nOMEGAs=0){
  PK <- SDE_lm(nTransit = nTransit)
  if(nOMEGAs==1){
    PK$h <- function(eta, theta, covar){
      phi <- theta
      phi$scal <- exp(theta$scal * covar[1])
      phi$stress <- theta$stress * exp(eta[1] -theta$sex * covar[2])
      phi$init <- (theta$ks/theta$ke)
      phi
    }
    PK$ModelPar <- function(THETA){
      list(
        theta=list(stress = THETA['stress'], ka = THETA['ka'], kt = THETA['kt'], ks= THETA['ks'], ke = THETA['ke'], init = THETA["init"], sigma=THETA['sigma'], S=THETA['S'], scal = THETA['scal'], sex = THETA["sex"]),
        OMEGA=matrix(THETA['OMEGA_stress'], ncol=1, nrow=1)
      )
    }
  }else if(nOMEGAs==2){
    PK$h <- function(eta, theta, covar){
      phi <- theta
      phi$scal <- exp(theta$scal * covar[1])
      phi$stress <- theta$stress * exp(eta[1] -theta$sex * covar[2])
      phi$ke <- theta$ke * exp(eta[2])
      phi$init <- (theta$ks/theta$ke)
      phi
    }
    PK$ModelPar <- function(THETA){
      list(
        theta=list(stress = THETA['stress'], ka = THETA['ka'], kt = THETA['kt'], ks= THETA['ks'], ke = THETA['ke'], init = THETA["init"], sigma=THETA['sigma'], S=THETA['S'], scal = THETA['scal'], sex = THETA["sex"]),
        OMEGA=diag(c(THETA['OMEGA_stress'], THETA['OMEGA_ke']))
      )
    }
  }else if(nOMEGAs==3){
    PK$h <- function(eta, theta, covar){
      phi <- theta
      phi$scal <- exp(theta$scal * covar[1])
      phi$stress <- theta$stress * exp(eta[1] -theta$sex * covar[2])
      phi$ke <- theta$ke * exp(eta[2])
      phi$kt <- theta$kt * exp(eta[3])
      phi$init <- (theta$ks/theta$ke)
      phi
    }
    PK$ModelPar <- function(THETA){
      list(
        theta=list(stress = THETA['stress'], ka = THETA['ka'], kt = THETA['kt'], ks= THETA['ks'], ke = THETA['ke'], init = THETA["init"], sigma=THETA['sigma'], S=THETA['S'], scal = THETA['scal'], sex = THETA["sex"]),
        OMEGA=diag(c(THETA['OMEGA_stress'], THETA['OMEGA_ke'], THETA['OMEGA_kt']))
      )
    }
  }else if(nOMEGAs==4){
    PK$h <- function(eta, theta, covar){
      phi <- theta
      phi$scal <- exp(theta$scal * covar[1])
      phi$stress <- theta$stress * exp(eta[1] -theta$sex * covar[2])
      phi$ke <- theta$ke * exp(eta[2])
      phi$kt <- theta$kt * exp(eta[3])
      phi$init <- (theta$ks/theta$ke)*exp(eta[4])
      phi
    }
    PK$ModelPar <- function(THETA){
      list(
        theta=list(stress = THETA['stress'], ka = THETA['ka'], kt = THETA['kt'], ks= THETA['ks'], ke = THETA['ke'], init = THETA["init"], sigma=THETA['sigma'], S=THETA['S'], scal = THETA['scal'], sex = THETA["sex"]),
        OMEGA=diag(c(THETA['OMEGA_stress'], THETA['OMEGA_ke'], THETA['OMEGA_kt'], THETA['OMEGA_init']))
      )
    }
  }
  return(PK)
}
ODE_lme <- function(nTransit=0, nOMEGAs=0){
  PK <- SDE_lme(nTransit = nTransit, nOMEGAs = nOMEGAs)
  PK$SIG <- function(phi) matrix(0, nrow=2+nTransit, ncol=2+nTransit)
  return(PK)
}

#revised function to estimate parameters and tune initial values
PSM.estim <- function(Model, Data, Par, trace = 0, control = NULL, fast = T){
  at_boundary <- TRUE
  while(at_boundary){
    fit <- PSM.estimate(Model = Model, Data = Data, Par = Par, CI = FALSE, trace = trace, control = control, fast = fast)
    
    if (!is.null(Par$LB)) {
      UBfac <- Par$UB[1]/Par$Init[1]
      LBfac <- Par$LB[1]/Par$Init[1]
      
      closeUB <- (fit$THETA/Par$UB) > 0.95
      Par$Init[which(closeUB)] <- Par$UB[which(closeUB)]
      Par$UB[which(closeUB)] <- Par$Init[which(closeUB)]*UBfac
      Par$LB[which(closeUB)] <- Par$Init[which(closeUB)]*LBfac
      
      closeLB <- (fit$THETA/Par$LB) < 1.05
      Par$Init[which(closeLB)] <- Par$LB[which(closeLB)]
      Par$LB[which(closeLB)] <- Par$Init[which(closeLB)]*LBfac
      Par$UB[which(closeLB)] <- Par$Init[which(closeLB)]*UBfac
      
      if(sum(closeUB) > 0 & trace > 0) cat(paste("re-fit: parameter(s)", paste(names(Par$UB[which(closeUB)]), collapse = ", "), "close to upper boundary.\n"))
      if(sum(closeLB) > 0 & trace > 0) cat(paste("re-fit: parameter(s)", paste(names(Par$LB[which(closeLB)]), collapse = ", "), "close to lower boundary.\n"))
      if(sum(closeLB,closeUB) == 0){
        if(trace > 0) cat("success: all parameters distant from boundaries.\n")
        at_boundary <- FALSE
      }
    }else{at_boundary <- FALSE}
  }
  return(fit)
}

#adjunct function to estimate or bootstrap confidence intervals for all parameters
PSM.ci <- function(PK.fit, PK.Model, PK.Data, bootstrap=NULL, optimizer="optim", cores=NULL, numderiv=T){
  if(!is.null(bootstrap)){
    #pb <- txtProgressBar(style = 3, min=1, max=bootstrap)
    
    boot <- function(PK.Model, PK.Data, PK.fit, bootstrap){
      res <- matrix(NA, ncol=length(PK.fit$THETA), nrow=bootstrap)
      colnames(res) <- names(PK.fit$THETA)
      for(i in 1:bootstrap){
        Pars <- list(Init=round(PK.fit$THETA,3), LB=round(PK.fit$THETA,3)*.5, UB=round(PK.fit$THETA,3)*1.5)
        samp <- sample(1:length(PK.Data), replace = T)
        tempdat <- Vectorize(function(x) list(PK.Data[[x]]))(samp)
        res[i,] <- try(PSM.estimate(PK.Model, tempdat, Pars, trace = 0, control=list(optimizer=optimizer))$THETA)
        #setTxtProgressBar(pb, i)
      }; return(res)
    }
    
    if(is.null(cores)){
      res <- boot(PK.Model, PK.Data, PK.fit, bootstrap)
    }else{
      registerDoMC(cores=cores)
      res <- foreach(n = 1:cores, .combine = rbind, .verbose = F) %dopar% {
        boot(PK.Model, PK.Data, PK.fit, bootstrap)
      }
    }
    
    ci <- apply(res,2, quantile, probs=c(.025,.5,.975))
    SD <- apply(res,2,sd)
    return(list(CI=ci, SD=SD, replicates=res))
  }else{
    tmpfun <- function(THETA, Model, Data, THETAnames) {
      names(THETA) <- THETAnames
      PSM:::APL.KF(THETA, Model, Data, GUIFlag = 0)
    }
    #tmpfun(as.numeric(PK.fit$THETA),  PK.Model, PK.Data, names(PK.fit$THETA))
    
    if(numderiv){
      Hess <- numDeriv::hessian(tmpfun, as.numeric(PK.fit$THETA), method = "Richardson", Model = PK.Model, Data = PK.Data, THETAnames = names(PK.fit$THETA))
    }else{
      Hess <- pracma::hessian(tmpfun, as.numeric(PK.fit$THETA), Model = PK.Model, Data = PK.Data, THETAnames = names(PK.fit$THETA))  
    }
    COV <- MASS::ginv(Hess)
    SD <- matrix(sqrt(diag(COV)), nrow = 1)
    SDmat <- t(SD) %*% SD
    COR <- 1/SDmat * COV
    colnames(SD) <- colnames(COR) <- rownames(COR) <- names(PK.fit$THETA)
    ci <- matrix(c(PK.fit$THETA - 1.96 * sqrt(diag(solve(Hess))), 
                   PK.fit$THETA, PK.fit$THETA + 1.96 * sqrt(diag(solve(Hess)))), 
                 nrow = 3, byrow = TRUE)
    rownames(ci) <- c("Lower CI95", "MLE", "Upper CI95")
    colnames(ci) <- names(PK.fit$THETA)
    return(list(CI = ci, SD = SD, COR = COR))
  }
}
PSM.ci <- cmpfun(PSM.ci)

#functions for visual predictive checking
mA <- function(x){ #calculate moving average
  y <- x
  for(i in 2:(length(x)-1) ) y[i] <- mean(x[(i-1):(i+1)])
  return(y)
}

library(vioplot) #truncated violinplot
vioplotm <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, horizontal = FALSE, col = "white", border = "black", lty = 1, lwd = 1, rectCol = "black", colMed = "black", pchMed = 19, at, add = FALSE, wex = 1, drawRect = T){
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- quantile(data, .05)
    data.max <- quantile(data, .95)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], rev(at[i] + height[[i]])), col = col, border = border, lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, q1 = q1, q3 = q3))
} 

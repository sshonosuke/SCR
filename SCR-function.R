###-----------------------------------------------------###
###  Functions for spatially clustered regression (SCR)  ###
###-----------------------------------------------------###
## This code implements the following two functions for SCR/SFCR
# 'SCR': SCR/SFCR with fixed G and fixed lambda
# 'SCR.select': tuning parameter (G and lambda) selection via information criteria

## packages
library(glmnet)
library(SparseM)


###  Spatially clustered regression (with LASSO)   ###
## IMPUT
# Y: n-dimensional response vector 
# X: (n,p)-matrix of covariates (p: number of covariates)
# W: (n,n)-matrix of spatial weight
# Sp: (n,2)-matrix of location information 
# lambda: tuning parameter for L1-reguralization (lambda=0 leads to standard regression) 
# G: number of groups 
# Phi: tuning parameter for spatial similarity 
# offset: n-dimensional vector of offset term (applicable only to "poisson" case)
# fuzzy: if True, SFCR is applied
# maxitr: maximum number of iterations
# family: distribution family ("gaussian" or "poisson")

## Output
# Beta: (G,p)-matrix of group-wise regression coeffieicnts
# Sig: G-dimensional vector of group-wise standard deviations (only for "gaussian" case)
# group: n-dimensional vector of group assignment 
# sBeta: (n,p)-matrix of location-wise regression coeffieicnts
# sSig: n-dimensional vector of location-wise standard deviations (only for "gaussian" case)
# ML: maximum log-likelihood
# itr: number of iterations 

## Remark
# matrix X should not include an intercept term
# initial grouping is determined by K-means of spatial locations


## Main function 
SCR <- function(Y, X, W, Sp, lambda=0, G=5, Phi=1, offset=NULL, fuzzy=F, maxitr=100, family="gaussian"){
  ## Preparations
  ep <- 10^(-5)      # convergence criterion 
  X <- as.matrix(X)
  n <- dim(X)[1]     # number of samples
  p <- dim(X)[2]+1    # number of regression coefficients
  XX <- as.matrix( cbind(1,X) )
  W <- as(W, "sparseMatrix")
  if(is.null(offset)){ offset <- rep(0, n) }
  
  ## Initial values
  Ind <- kmeans(Sp, G)$cluster
  Pen <- rep(0, G)
  Beta <- matrix(0, p, G)
  dimnames(Beta)[[2]] <- paste0("G=",1:G)
  Sig <- rep(1, G)    # not needed under non-Gaussian case
  
  ## iterative algorithm 
  val <- 0
  mval <- 0
  for(k in 1:maxitr){
    cval <- val
    
    ## penalty term
    Ind.mat <- matrix(0, n, G)
    for(g in 1:G){
      Ind.mat[Ind==g, g] <- 1
    }
    Ind.mat <- as(Ind.mat, "sparseMatrix")
    Pen <- W%*%Ind.mat     # penalty term
    
    ## model parameters (clustered case)
    if(fuzzy==F){
      for(g in 1:G){
        if(length(Ind[Ind==g])>p+1){
          # gaussian
          if(family=="gaussian"){
            if(lambda==0){
              fit <- lm(Y[Ind==g]~X[Ind==g,])
            }else{
              fit <- glmnet(x=X[Ind==g,], y=Y[Ind==g], family="gaussian", lambda=lambda)
            }
            Beta[,g] <- as.vector( coef(fit) )
            resid <- Y-as.vector(XX%*%Beta[,g])
            Sig[g] <- sqrt(mean(resid[Ind==g]^2))
            Sig[g] <- max(Sig[g], 0.1)
          }
          # poisson
          if(family=="poisson"){
            fit <- glmnet(x=X[Ind==g,], y=Y[Ind==g], offset=offset[Ind==g], family="poisson", lambda=lambda)
            Beta[,g] <- as.vector( coef(fit) )
          }
        }
      }
    }
    
    ## model parameters (fuzzy case)
    if(fuzzy==T){
      # Gaussian
      if(family=="gaussian"){
        Mu <- XX%*%Beta      # (n,G)-matrix
        ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
        log.dens <- log(dnorm(Y,Mu,ESig)) + Phi*Pen
        mval <- apply(log.dens, 1, max)
        log.denom <- mval + log(apply(exp(log.dens-mval), 1, sum))
        PP <- exp(log.dens-log.denom)     # weight
        for(g in 1:G){ 
          if(sum(PP[,g])>0.1){
            fit <- glmnet(x=X, y=Y, lambda=lambda, weights=PP[,g], family="gaussian")
            Beta[,g] <- as.vector( coef(fit) )
            resid <- Y-as.vector(XX%*%Beta[,g])
            Sig[g] <- sqrt( sum(PP[,g]*resid^2)/sum(PP[,g]) )
            Sig[g] <- max(Sig[g], 0.1)
          }
        }
      }
      
      # Poisson
      if(family=="poisson"){
        Mu <- exp(offset + XX%*%Beta)    # (n,G)-matrix
        log.dens <- log(dpois(Y, Mu)) + Phi*Pen 
        mval <- apply(log.dens, 1, max)
        log.denom <- mval + log(apply(exp(log.dens-mval), 1, sum))
        PP <- exp(log.dens-log.denom)     # weight
        for(g in 1:G){
          if(sum(PP[,g])>0.1){
            fit <- glmnet(x=X, y=Y, offset=offset, lambda=lambda, weights=PP[,g], family="poisson")
            Beta[,g] <- as.vector( coef(fit) )
          }
        }
      }
    }
    
    ## Grouping (clustered case)
    if(fuzzy==F){
      # Gaussian
      if(family=="gaussian"){
        Mu <- XX%*%Beta      # (n,G)-matrix
        ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
        Q <- dnorm(Y, Mu, ESig, log=T) + Phi*Pen     # penalized likelihood
      }
      if(family=="poisson"){ 
        Mu <- exp(offset + XX%*%Beta) 
        Q <- dpois(Y, Mu, log=T)  + Phi*Pen    # penalized likelihood
      }
      Ind <- apply(Q, 1, which.max)
    }
    
    ## Grouping (fuzzy case)
    if(fuzzy==T){
      # Gaussian
      if(family=="gaussian"){
        Mu <- XX%*%Beta      # (n,G)-matrix
        ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
        Q <- dnorm(Y, Mu, ESig, log=T) + Phi*Pen   # penalized likelihood
        mval <- apply(Q, 1, max)
        log.denom <- mval + log(apply(exp(Q-mval), 1, sum))
        PP <- exp(Q-log.denom)
      }
      # Poisson
      if(family=="poisson"){ 
        Mu <- exp(offset + XX%*%Beta)    # (n,G)-matrix
        Q <- dpois(Y, Mu, log=T)  + Phi*Pen   # penalized likelihood
        mval <- apply(Q, 1, max)
        log.denom <- mval + log(apply(exp(Q-mval), 1, sum))
        PP <- exp(Q-log.denom)
      }
      Ind <- apply(PP, 1, which.max)
    }
    
    ## Value of objective function
    val <- sum( apply(Q, 1, max) ) + lambda*sum(abs(Beta))
    dd <- abs(cval-val)/abs(val)
    mval <- max(mval, cval)
    if( dd<ep | abs(mval-val)<ep ){ break }
  }
  
  ## varying parameters
  if(fuzzy==F){  sBeta <- t(Beta[,Ind]) }
  if(fuzzy==T){  sBeta <- PP%*%t(Beta) }
  sSig <- Sig[Ind]     # location-wise error variance

  ## maximum likelihood 
  if(family=="gaussian"){
    hmu <- apply(XX*sBeta, 1, sum)
    ML <- sum( dnorm(Y, hmu, sSig, log=T) ) 
  }
  if(family=="poisson"){
    hmu <- exp(offset + apply(XX*sBeta, 1, sum))
    ML <- sum( dpois(Y, hmu, log=T) ) 
  }
  
  ## Results
  result <- list(Beta=Beta, Sig=Sig, group=Ind, sBeta=sBeta, sSig=sSig, ML=ML, itr=k)
  return(result)
}










###  Selection of tuning parameters  ###
## Imput
# most of imputs are the same as 'SCR' 
# G.set: vector of candidates for G
# Lam.set: vector of candidates for lambda
# print: if True, interim progress is reported 

## Output
# AIC: AIC-type criteria
# BIC: BIC-type criteria
# select: selection results 

## Main function 
SCR.select <- function(Y, X, W, Sp, G.set=NULL, Lam.set=NULL, Phi=1, offset=NULL, maxitr=50, print=T, family="gaussian"){
  ## Preparations
  if(is.null(G.set)){ G.set <- seq(10, 40, by=5) }
  if(is.null(Lam.set)){ Lam.set <- 0 }
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]+1
  L <- length(G.set)
  M <- length(Lam.set)
  
  ## computing information criteria
  AIC <- matrix(NA, M, L)
  BIC <- matrix(NA, M, L)
  for(m in 1:M){
    for(l in 1:L){
      fit <- SCR(Y, X, W, Sp, offset=offset, G=G.set[l], lambda=Lam.set[m], Phi=Phi, maxitr=maxitr, family=family)
      pp <- sum(fit$Beta!=0) 
      if(family=="gaussian"){ pp <- pp + G.set[l] }
      AIC[m, l] <- -2*fit$ML + 2*pp
      BIC[m, l] <- -2*fit$ML + log(n)*pp
      if(print){ print( paste0("G=",G.set[l], ", Lam=", Lam.set[m], ", iteration=", fit$itr) ) }
    }
  }
  dimnames(AIC)[[1]] <- dimnames(BIC)[[1]] <- paste0("lambda=", Lam.set)
  dimnames(AIC)[[2]] <- dimnames(BIC)[[2]] <- paste0("G=", G.set)
  
  ## selection
  if(M>1){
    num <- which.min(AIC)
    G.AIC <- G.set[ceiling(num/M)]
    Lam.AIC <- Lam.set[num-floor(num/M)*M]
    num <- which.min(BIC)
    G.BIC <- G.set[ceiling(num/M)]
    Lam.BIC <- Lam.set[num-floor(num/M)*M]
    select <- c(G.AIC, Lam.AIC, G.BIC, Lam.BIC)
    names(select) <- c("AIC.G", "AIC.Lam", "BIC.G", "BIC.Lam")
  }
  if(M==1){
    G.AIC <- G.set[which.min(AIC)]
    G.BIC <- G.set[which.min(BIC)]
    select <- c(G.AIC, G.BIC)
    names(select) <- c("AIC.G", "BIC.G")
  }
  
  ## result
  IC <- list(AIC=AIC, BIC=BIC, select=select)
  return(IC)
}



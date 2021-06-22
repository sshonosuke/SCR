###-----------------------------------------------------###
###  Functions for spatially clustered regression (SCR)  ###
###-----------------------------------------------------###
## This code implements the following two functions for SCR/SFCR
# 'SCR': SCR/SFCR with fixed G 
# 'SCR.select': tuning parameter (G) selection via BIC-type criteria

## packages
library(SparseM)
library(MASS)


###  Spatially clustered regression (with LASSO)   ###
## Input
# Y: n-dimensional response vector 
# X: (n,p)-matrix of covariates (p: number of covariates)
# W: (n,n)-matrix of spatial weight
# Sp: (n,2)-matrix of location information 
# G: number of groups 
# Phi: tuning parameter for spatial similarity 
# offset: n-dimensional vector of offset term (applicable only to "poisson" and "NB")
# fuzzy: if True, SFCR is applied
# maxitr: maximum number of iterations
# family: distribution family ("gaussian", "poisson" or "NB)

## Output
# Beta: (G,p)-matrix of group-wise regression coefficients
# Sig: G-dimensional vector of group-wise standard deviations (only for "gaussian" and "NB")
# group: n-dimensional vector of group assignment 
# sBeta: (n,p)-matrix of location-wise regression coefficients
# sSig: n-dimensional vector of location-wise standard deviations (only for "gaussian")
# s: n-dimensional vector of location-wise standard deviations (only for "gaussian")
# ML: maximum log-likelihood
# itr: number of iterations 

## Remark
# matrix X should not include an intercept term
# initial grouping is determined by K-means of spatial locations


## Main function 
SCR <- function(Y, X, W, Sp, G=5, Phi=1, offset=NULL, fuzzy=F, maxitr=100, delta=1, family="gaussian"){
  ## Preparations
  ep <- 10^(-5)      # convergence criterion 
  X <- as.matrix(X)
  n <- dim(X)[1]     # number of samples
  p <- dim(X)[2]+1    # number of regression coefficients
  XX <- as.matrix( cbind(1,X) )
  W <- as(W, "sparseMatrix")
  if(is.null(offset)){ offset <- rep(0, n) }
  nmax <- function(x){ max(na.omit(x)) }   # new max function 
  
  ## Initial values
  M <- 20  # the number of initial values of k-means
  WSS <- c()
  CL <- list()
  for(k in 1:M){
    CL[[k]] <- kmeans(Sp, G)
    WSS[k] <- CL[[k]]$tot.withinss
  }
  Ind <- CL[[which.min(WSS)]]$cluster
  Pen <- rep(0, G)
  Beta <- matrix(0, p, G)
  dimnames(Beta)[[2]] <- paste0("G=",1:G)
  Sig <- rep(1, G)    # not needed under non-Gaussian case
  Nu <- rep(1, G) 
  
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
            fit <- lm(Y[Ind==g]~X[Ind==g,])
            Beta[,g] <- as.vector( coef(fit) )
            resid <- Y-as.vector(XX%*%Beta[,g])
            Sig[g] <- sqrt(mean(resid[Ind==g]^2))
            Sig[g] <- max(Sig[g], 0.1)
          }
          # poisson
          if(family=="poisson"){
            x <- X[Ind==g,]
            y <- Y[Ind==g]
            off <- offset[Ind==g]
            fit <- glm(y~x, offset=off, family="poisson")
            Beta[,g] <- as.vector( coef(fit) )
          }
          # NB
          if(family=="NB"){
            x <- X[Ind==g,]
            y <- Y[Ind==g]
            off <- offset[Ind==g]
            fit <- glm.nb(y~x+offset(off))
            Beta[,g] <- as.vector( coef( fit ) )
            Nu[g] <- fit$theta
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
            fit <- lm(Y~X, weights=PP[,g])
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
            fit <- glm(Y~X, offset=offset, weights=PP[,g], family="poisson")
            Beta[,g] <- as.vector( coef(fit) )
          }
        }
      }
      # NB
      if(family=="NB"){
        Mu <- exp(offset + XX%*%Beta)    # (n,G)-matrix
        log.dens <- dnbinom(Y, size=Nu, prob=Nu/(Nu+Mu), log=T) + Phi*Pen 
        mval <- apply(log.dens, 1, max)
        log.denom <- mval + log(apply(exp(log.dens-mval), 1, sum))
        PP <- exp(log.dens-log.denom)     # weight
        for(g in 1:G){
          if(sum(PP[,g])>0.1){
            fit <- glm.nb(Y~X+offset(offset), weights=PP[,g])
            Beta[,g] <- as.vector( coef(fit) )
            Nu[g] <- fit$theta
          }
        }
      }
    }
    
    ## Grouping (clustered case)
    if(fuzzy==F){
      if(family=="gaussian"){
        Mu <- XX%*%Beta      # (n,G)-matrix
        ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
        Q <- dnorm(Y, Mu, ESig, log=T) + Phi*Pen     # penalized likelihood
      }
      if(family=="poisson"){ 
        Mu <- exp(offset + XX%*%Beta) 
        Q <- dpois(Y, Mu, log=T)  + Phi*Pen    # penalized likelihood
      }
      if(family=="NB"){ 
        Mu <- exp(offset + XX%*%Beta) 
        Q <- dnbinom(Y, size=Nu, prob=Nu/(Nu+Mu), log=T)  + Phi*Pen    # penalized likelihood
      }
      Ind <- apply(Q, 1, which.max)
    }
    
    ## Grouping (fuzzy case)
    if(fuzzy==T){
      if(family=="gaussian"){
        Mu <- XX%*%Beta      # (n,G)-matrix
        ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
        Q <- delta*(dnorm(Y, Mu, ESig, log=T) + Phi*Pen)   # penalized likelihood
        mval <- apply(Q, 1, max)
        log.denom <- mval + log(apply(exp(Q-mval), 1, sum))
        PP <- exp(Q-log.denom)
      }
      if(family=="poisson"){ 
        Mu <- exp(offset + XX%*%Beta)    # (n,G)-matrix
        Q <- delta*(dpois(Y, Mu, log=T)  + Phi*Pen)   # penalized likelihood
        mval <- apply(Q, 1, max)
        log.denom <- mval + log(apply(exp(Q-mval), 1, sum))
        PP <- exp(Q-log.denom)
      }
      if(family=="NB"){ 
        Mu <- exp(offset + XX%*%Beta)    # (n,G)-matrix
        Q <- delta*(dnbinom(Y, size=Nu, prob=Nu/(Nu+Mu), log=T) + Phi*Pen)   # penalized likelihood
        mval <- apply(Q, 1, max)
        log.denom <- mval + log(apply(exp(Q-mval), 1, sum))
        PP <- exp(Q-log.denom)
      }
      Ind <- apply(PP, 1, which.max)
    }
    
    ## Value of objective function
    val <- sum( apply(Q, 1, nmax) ) 
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
  if(family=="NB"){
    sNu <- Nu[Ind]
    hmu <- exp(offset + apply(XX*sBeta, 1, sum))
    ML <- sum( dnbinom(Y, size=sNu, prob=sNu/(sNu+hmu), log=T) ) 
  }
  
  ## Results
  result <- list(Beta=Beta, Sig=Sig, Nu=Nu, group=Ind, sBeta=sBeta, sSig=sSig, ML=ML, itr=k)
  return(result)
}










###  Selection of tuning parameters  ###
## Imput
# most of inputs are the same as 'SCR' 
# G.set: vector of candidates for G
# print: if True, interim progress is reported 

## Output
# BIC: BIC-type criteria
# select: selection results 

## Main function 
SCR.select <- function(Y, X, W, Sp, G.set=NULL, Phi=1, offset=NULL, maxitr=50, print=T, family="gaussian"){
  ## Preparations
  if(is.null(G.set)){ G.set <- seq(10, 40, by=5) }
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]+1
  L <- length(G.set)
  
  ## computing information criteria
  BIC <- c()
  for(l in 1:L){
    fit <- SCR(Y, X, W, Sp, offset=offset, G=G.set[l], Phi=Phi, maxitr=maxitr, family=family)
    pp <- length(fit$Beta) 
    if(family=="gaussian"| family=="NB"){ pp <- pp + G.set[l] }
    BIC[l] <- -2*fit$ML + log(n)*pp
    if(print){ print( paste0("G=",G.set[l], ", iteration=", fit$itr) ) }
  }
  names(BIC) <- paste0("G=", G.set)
  
  ## selection
  hG <- G.set[which.min(BIC)]

  ## result
  Result <- list(BIC=BIC, G=hG)
  return(Result)
}



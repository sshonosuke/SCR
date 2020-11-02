###-----------------------------------------------------###
###     one-shot simulation experiment (fitting)        ###
###-----------------------------------------------------###
set.seed(1)

library(spgwr)   
library(glmnet) 
library(ape)    
library(MASS)
source("SCR-function.R")


## load simulated dataset
scenario <- 2    # 1 (clustered case) or 2 (smoothed case)
load(paste0("simdata", scenario, ".RData"))


## GWR: geographically weighted regression
b.opt <- gwr.sel(Y~X, coords=Sp, verbose=F)
fit <- gwr(Y~X, coords=Sp, bandwidth=b.opt)
est0 <- cbind(fit$SDF$`(Intercept)`, fit$SDF$`Xx1`, fit$SDF$`Xx2`)

## SCR: spatially clustered regression
G.set <- seq(5, 30, by=5)
IC <- SCR.select(Y, X, W, Sp, G.set=G.set, Phi=1, print=F, family="gaussian")
hG <- IC$select[2]     # optimal number of group via BIC

# SCR
fit1 <- SCR(Y, X, W, Sp, G=hG, Phi=1)    
est1 <- fit1$sBeta

# SFCR
fit2 <- SCR(Y, X, W, Sp, G=hG, Phi=1, fuzzy=T)
est2 <- fit2$sBeta


## SHP: spatially homogeneity pursuit method
n <- length(Y)
MST <- mst(dist(Sp))    # minimum spanning tree
H <- c()
for(i in 1:n){
  for(j in 1:n){
    if(i<j){
      if(MST[i,j]==1){
        h <- rep(0, n)
        h[i] <- 1
        h[j] <- -1
        H <- rbind(H, h)
      }
    }
  }
}
HH <- rbind(H, rep(1/n, n))
invH <- solve(HH)
XX <- cbind(diag(n), diag(X[,1]), diag(X[,2]))
invHH <- matrix(0, 3*n, 3*n)
for(k in 1:3){
  sub <- (1+n*(k-1)):(n*k)
  invHH[sub, sub] <- invH
}
XH <- XX%*%invHH
pen <- rep(c(rep(1, n-1), 0), 3)
fit.SHP <- glmnet(x=XH, y=Y, family="gaussian", intercept=F, penalty.factor=pen)
BIC <- c()
for(j in 1:100){
  mu <- as.vector(XH%*%as.vector(fit.SHP$beta[,j]))
  sig <- sqrt(mean((Y-mu)^2))
  BIC[j] <- -2*sum(dnorm(Y, mu, sig, log=T)) + log(n)*fit.SHP$df[j]
}
opt <- which.min(BIC)
est <- as.vector(invHH%*%as.vector(fit.SHP$beta[,opt]))
est.SHP <- matrix(est, n, 3)



## mean squared errors
# GWR
apply((est0-Beta)^2, 2, mean)
# SHP
apply((est.SHP-Beta)^2, 2, mean)
# SCR
apply((est1-Beta)^2, 2, mean)
# SFCR
apply((est2-Beta)^2, 2, mean)


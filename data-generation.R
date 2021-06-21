###-----------------------------------------------------###
###   one-shot simulation experiment (data generation)  ###
###-----------------------------------------------------###
set.seed(1)

library(MASS)
scenario <- 2   # 1 (clustered case) or 2 (smoothed case)
phi <- 0.6    # range parameter for covariates
n <- 1000   # sample size


##  generation of sampling locations
Sp <- matrix(NA,n,2)
for(i in 1:n){
  dd <- 0
  while(dd<(0.5^2)){
    rn <- c(runif(1,-1,1),runif(1,0,2))
    dd <- sum(rn[1]^2+0.5*rn[2]^2)
  }
  Sp[i,] <- rn
}



## covariates
dd <- as.matrix(dist(Sp))
mat <- exp(-dd/phi)

z1 <- mvrnorm(1, rep(0, n), mat)
z2 <- mvrnorm(1, rep(0, n), mat)
x1 <- z1
rr <- 0.75
x2 <- rr*z1 + sqrt(1-rr^2)*z2
X <- cbind(x1,x2)

## correlation matrix for regression coefficients
V1 <- exp(-dd/1)
V2 <- exp(-dd/2)
V3 <- exp(-dd/3)


## Weight matrix for SCR
DD <- as.matrix(dist(Sp))
th <- 0.1
W <- matrix(0, n, n)
for(i in 1:n){
  W[i,] <- ifelse(DD[i,]<th, 1, 0)
}
diag(W) <- 0



# Parameters (clustered case)
if(scenario==1){
  Beta0 <- c()
  Beta1 <- c()
  Beta2 <- c()
  Sig <- c()
  M1 <- 2
  M2 <- 3
  Grid1 <- seq(-1, 1, length=M1+1)
  Grid2 <- seq(0, 2, length=M2+1)
  for(k in 1:M1){
    for(j in 1:M2){
      Ind <- (1:n)[Grid1[k]<Sp[,1] & Sp[,1]<Grid1[k+1] & Grid2[j]<Sp[,2] & Sp[,2]<Grid2[j+1]]
      Beta0[Ind] <- 2*(Grid1[k] + Grid2[j])
      Beta1[Ind] <- (Grid1[k]^2 + Grid2[j]^2)
      Beta2[Ind] <- -(Grid1[k] + Grid2[j])
      Sig[Ind] <- 0.5 + 0.2*abs(Grid1[k] - Grid2[j])
    }
  }
  Beta <- cbind(Beta0, Beta1, Beta2)
}


# Parameters (smoothed case)
if(scenario==2){
  Beta0 <- mvrnorm(1, rep(0,n), 2*V1)
  Beta1 <- mvrnorm(1, rep(0,n), 2*V2)
  Beta2 <- mvrnorm(1, rep(0,n), 2*V3)
  Sig <- 0.2*exp(mvrnorm(1, rep(0,n), 0.5*V3))
  Beta <- cbind(Beta0, Beta1, Beta2)
}




## Data
Mu <- apply(cbind(1,X)*Beta, 1, sum)
Y <- rnorm(n, Mu, Sig)



## save
save(Y, X, W, DD, Sp, Beta, file=paste0("simdata", scenario, ".RData"))
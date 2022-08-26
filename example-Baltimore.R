###   Example: Boltimore data   ###
rm(list=ls())

## dataset
library(spdep)
library(spatialreg)
K <- 5    # number of nearest neighbour
Y <- as.vector(baltimore$PRICE)
X <- as.matrix(cbind(baltimore$NROOM, baltimore$DWELL, baltimore$NBATH, baltimore$PATIO, 
                     baltimore$FIREPL, baltimore$AC, baltimore$BMENT, baltimore$NSTOR,
                     baltimore$GAR, baltimore$AGE, baltimore$CITCOU, baltimore$LOTSZ, baltimore$SQFT))
Sp <- as.matrix(cbind(baltimore$X, baltimore$Y))
knn <- knn2nb(knearneigh(Sp, k=K))
listw_10nn_dates <- nb2listw(knn)
W <- as(as_dgRMatrix_listw(listw_10nn_dates), "CsparseMatrix")*K


## fitting
source("SCR-function.R")
select <- SCR.select(Y, X, W, Sp, G.set=2:10)   # selection of the number of clusters
select$G

fit <- SCR(Y, X, W, Sp, G=select$G)  # estimation with the selected G

plot(Sp, col=fit$group)  




  

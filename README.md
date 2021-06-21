# SCR: Spatially Clustered Regression

This repository provides R code implementing spatially clustered regression for spatial data analysis, as proposed by the following paper.

[Sugasawa, S. and Murakami, D. (2021). Spatially Clustered Regression. *Spatial Statistics* in press.](https://doi.org/10.1016/j.spasta.2021.100525)

(arXiv version: https://arxiv.org/abs/2011.01493)

The repository includes the following 5 files.

* SCR-function.R : Script implementing the proposed method
* simulation.R : Script applying the proposed method to two simulated datasets 
* simdata1.RData: Simulated data 1
* simdata2.RData: Simulated data 2
* data-generation.R: Script for generating the two simulated datasets
* example-Boltimore.R: Script for applying SCR to boltimore dataset (available from `spdep` package)

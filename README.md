# CSurf: Change Surface (CSurf) regression method for nonlinear subgroup identification

Pharmacogenomics stands as a pivotal driver toward personalized medicine, aiming to optimize drug efficacy while minimizing adverse effects by uncovering the impact of genetic variations on inter-individual outcome variability. Despite its promise, the intricate landscape of drug metabolism introduces complexity, where the correlation between drug response and genes can be shaped by numerous non-genetic factors, often exhibiting heterogeneity across diverse subpopulations. This challenge is particularly pronounced in datasets like the International Warfarin Pharmacogenetic Consortium (IWPC), which encompasses diverse patient information from multiple nations. To capture the between-patient heterogeneity in dosing requirement, we formulate a novel change \tcr{surface} model as a model-based approach for multiple subgroup identification in complex datasets. A key feature of our approach is its ability to accommodate nonlinear subgroup divisions, providing a clearer understanding of dynamic drug-gene associations. Furthermore, our model effectively handles high-dimensional data through a doubly penalized approach, ensuring both interpretability and adaptability. We propose an iterative two-stage method that combines a change point detection technique in the first stage with a smoothed local adaptive majorize-minimization algorithm for \tcr{surface} regression in the second stage. Performance of the proposed methods is evaluated through extensive numerical studies. Application of our method to the IWPC dataset leads to significant new findings, where three subgroups subject to different pharmacogenomic relationships are identified, contributing valuable insights into the complex dynamics of drug-gene associations in patients.

## Installation
You can install CSurf from GitHub using devtools. It should install in 1 minute, but may take longer if you need to update dependencies.  

``` r
library(devtools)
devtools::install_github("chrystal23/CSurf")
```

## Tutorial

The main function CSurf requires the following inputs:

### Inputs

- `y`: A vector (length n) of the response for subgroup identification.
- `X`: A matrix (n x p) of the covariates for subgroup identification. The covariates are also used as potential thresholding variables for subgroup identification.
- `f_hat0`: A vector (length n) of the initialization for the threshold function (change surface).
- `Tighten`: A boolean indicating whether to include the tightening stage in optimization.
- `tol_out`: The tolerance value for optimization (stopping criterion).
- `maxiter_out`: The maximal number of iterations allowed for optimization.
- `m`: An integer indicating the power of the threshold function.
- `lam2rate`: The ratio between the two tunning parameters, lambda1 and lambda2.

The CSurf function returns a list of estimation results: 

 ### Outputs

- `f_all.hat`: The estimated threshold functions. Each column represent one component function.
- `tau.hat`: The estimated change points over the threshold functions.
- `alp.hat`: The estimated regression coefficients.
- `X_ord`: the re-ordered covariate matrix, which corresponds to above estimation results
- `y_ord`: the re-ordered response vector, which corresponds to above estimation results
- `coverged_out`: An boolean indicating whether the algorithm successfully converged.
- `iter_out`: The number of iterations it took to converge.
- `...`: Other estimation results.

An example input data is can be loaded using
```r
data(syn_data)
```
where an example response vector (length 500) `y`, a covariate matrix (500 x 30) `X`, and a vector of threshold function initialization (length 500) `f_hat0` are included. To implement our CSurf method on this example data to identify subgroups, we run
```r
results <- CSurf(y = y, X = X, f_hat0 = f_hat0, h = .5, maxiter_out = 10, tol_out = 0.0001)
```
where `h = .5, maxiter_out = 10, tol_out = 0.0001` is exactly the default setting but can be specified by users based on their requirement for the optimization.
The estimation results are
```r
f_all.hat = results$f_all.hat
alp.hat<- results$alp.hat
tau.hat<- results$tau.hat
Xord = results$Xord
```
The component of threshold functions that characterize the subgroups can be visualized:
```r
par(mfrow = c(1,3))
for (j in 1:3) {
   plot(y = f_all.hat[order(Xord[,j]),j], 
        x = sort(Xord[,j]), 
        ylim = range(f_all.hat[,j])*2,
        type = "l", 
        xlab = paste("f", j, sep =  "_"), 
        ylab = paste("f", j, sep =  "_")) 
 }
```
<img src="man/figures/results_fhat.png" width="80%" />

We can obtain the subgroup membership, denoted by `subgr.hat`, of each observation (subject) from the results
```r
subgr.hat <- rep(1,dim(Xord)[1])
for (k in 1:length(tau.hat)) {
  id <- which( rowSums(f_all.hat) > tau.hat[k] )
  subgr.hat[id] <- subgr.hat[id] + 1
}
table(subgr.hat)
```
| subgr.hat | 1 | 2 | 3 |
|---|-------------|----------|-------------------------|
|  | 116 | 196 | 188 |



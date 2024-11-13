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
- 
 f_hat0 A vector (length n) of the initialization for the threshold function (change surface).
 Tighten A boolean indicating whether to include the tightening stage in optimization.
tol_out The tolerance value for optimization (stopping criterion).
maxiter_out The maximal number of iterations allowed for optimization.
m An integer indicating the power of the threshold function.
lam2rate The ratio between the two tunning parameters, lambda1 and lambda2.

The CSurf function returns a list of estimation results: 

 ### Outputs

#' \item{f_all.hat}{The estimated threshold functions. Each column represent one component function.}
#' \item{tau.hat}{The estimated change points over the threshold functions.}
#' \item{alp.hat}{The estimated regression coefficients.}
#' \item{coverged_out}{An boolean indicating whether the algorithm successfully converged.}
#' \item{iter_out}{The number of iterations it took to converge.}
#' \item{...}{Other estimation results.}

An example input data is can be loaded using
```r
data(syn_data)
results <- CSurf(y = y, X = X, f_hat0 = f_hat0, h = .5, maxiter_out = 10, tol_out = 0.01)
```
 

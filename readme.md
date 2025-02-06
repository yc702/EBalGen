
# EBalGen: Entropy balancing for covariate shift causal generalization

## Overview

### Overview and introduction

The R package *EBalGen* provides tools for implementing exact and
approximate balancing methods for causal generalization based on the
work of Chen, Chen and Yu (2023) and Chen, Chen and Yu (2025+). The key
idea is that for causal generalization, differences in the distributions
of treatment effect modifiers between these populations, known as
covariate shift, can lead to varying ATEs. Our methods use only
summary-level information from a target sample while accounting for the
possible covariate shifts.

### Key Functionalities

*EBalGen* provides the following functionalities:

- Weight: Computes both exact and approximate balancing weights to
  account for covariate differences between the source and target
  populations and within source population.

- ATE: Estimates the average treatment effect (ATE) using the computed
  weights.

- Confidence Interval (CI): Constructs confidence intervals for the
  estimated ATE using resampling-based perturbation method (RPM).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yc702/EBalGen")
```

Load the package with

``` r
library(EBalGen)
```

## Example 1: Source and target population have a good overlap.

### Data Set up

We set the total sample size $n = n_s +n_t= 800$ and is split into
source $n_s = 401$ and target $n_t = 399$ samples. We generate 5
covariates $X = (X_1, \ldots, X_5)$ from a uniform distribution
$U(-2, 2)$. The source/target participation probability $\rho (x)$
follows $\text{logit} \{ \rho (x) \} = 0.4x_1+ 0.3x_2-0.2x_4$. That is,
there is shift in the distribution of $(X_1, X_2, X_4)$.

Propensity score $(\pi(x))$ model, we assume when the treatment
assignment is related to $H$ linearly with
$\text{logit}\{\pi(x)\} = 0.7x_2 + 0.5 x_3$. In this case, all the
confounders are included in $H$, and it is enough that we only balance
on $H$ to account for confounding.

Outcome model, we assume
$Y_i = m(X_i) + (A_i-0.5) \tau(X_i) + \epsilon_i$ with
$\epsilon_i \overset{\text{i.i.d}}{\sim} N(0,1)$.

CATE function we assume \$(x) = x_1 - 0.6x_2 - 0.4 x_3 \$. For the main
effect $m(x)$, it has the form of
$m(x) = 0.5x_1 + 0.3x_2 +0.3x_3 - 0.4 x_4 - 0.7 x_5$.

In this setting, the Average treatment effect for target (ATET) is -0.14

### Visual check of the propensity scores between source and target samples

Here is the density plot of the propensity scores of source and target
samples fitted using simple logistic regression using all 5 covariates.
The distribution of propensity scores in both populations shows a
substantial degree of overlap, indicating that the covariate
distributions between the two samples are sufficiently similar. This
overlap suggests that the generalization of treatment effects from the
source population to the target population is more reliable and exact
balancing could be achieved.

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%"  />

### Compute exact balancing weights and generalized ATE

For causal generalization, assume in the target sample we only have
summary-level information of $X_1, X_2,$ and $X_3$. We set
$H(x) = (1, x_1, x_2, x_3)$ and $G(x) = (x_4, x_5)$. We consider
balancing on the first moments of $X_1, X_2,$ and $X_3$.

**Data input:**

- `xs` A data matrix for the source sample. Each column represents
  source sample covariate and each row represents an observation.

- `ys` A vector of the source sample response values.

- `trts` A vector of 0, 1 or TRUE/FALSE of treatment assignment for the
  source sample.

- `H_vars` A vector of numbers specifying which covariate in the source
  sample need to be balanced between source and target samples. Here we
  balance on 1,2,3 covariates.

- `target_moments`A vector of first moments of the target sample
  covariates that needs to be balanced between source and target. Here
  the values are -0.2390874, -0.1870105, -0.0581574.

- `delta` A vector specifying the approximate balancing tolerance
  margin. Here the values are (0, 0, 0, 0, 0, 0, 0, 0). For the delta’s
  length of 8, it is balancing between source and target and within
  source: 3+3+2.

Here is the summary statistics of the weights for the source sample.

``` r
library(EBalGen)

## Source sample
wts_gen <- ebal_wts(xs, trts,H_vars, target_moments,H_add_intercept = TRUE,delta )$w
summary(wts_gen)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2741  0.9978  1.5188  2.0000  2.6778 11.5155
```

Here is the generalized ATE of source sample on the target sample

``` r

ebal_ATE(xs,ys,trts,H_vars, target_moments, H_add_intercept=TRUE,delta)$ATE
#>      value 
#> -0.1051732
```

### Compute exact balancing CI

Use resampling-based perturbation `RPM_CI()` with additional input of
`target_sd` and the number of bootstrap iteration 300.

- `target_sd` A vector of standard deviation of the target sample
  covariates that needs to be balanced between source and target.

``` r
## CI construction
target_sd = colStdevs(xt)[H_vars]
ATE_CI = RPM_CI(xs, ys,trts, H_vars=H_vars,target_mean=target_moments,target_sd=target_sd,num_sim=300, H_add_intercept=TRUE,cluster=5, set_seed=100)

## Lower bound of 95% CI
ATE_CI$lb_ATE
#>      2.5% 
#> -1.906325

## Upper bound of 95% CI
ATE_CI$ub_ATE
#>    97.5% 
#> 1.918743
```

------------------------------------------------------------------------

## Example 2: Source and target population have bad overlap

### Data Set up

We set the total sample size $n = n_s +n_t= 400$ and is split into
source $n_s = 281$ and target $n_t = 119$ samples. We generate 5
covariates $X = (X_1, \ldots, X_5)$ from a uniform distribution
$U(-2, 6)$. The remaining settings are identical to those in the
previous example.

In this setting, the Average treatment effect for target (ATET) is -0.64

### Visual check of the propensity scores between source and target samples

Here is the density plot of the propensity scores of source and target
samples fitted using simple logistic regression using all 5 covariates.
The distribution of propensity scores in both populations shows a
limited degree of overlap, indicating that the covariate distributions
between the two samples are quite different. This overlap suggests that
approximate balancing should be used.

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%"  />

### Compute exact balancing weights and generalized ATE

For data input, everything are the same except we now have the
approximate balancing tolerance margin `delta` as all 0.1:(0.1, 0.1,
0.1, 0.1, 0.1, 0.1, 0.1, 0.1).

Here is the summary statistics of the weights for the source sample.

``` r
library(EBalGen)

## Source sample

wts_gen <- ebal_wts(xs, trts,H_vars, target_moments,H_add_intercept = TRUE,delta=numeric(8)+0.1 )$w
summary(wts_gen)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>  0.00104  0.17267  0.47180  2.00001  1.42941 76.36493
```

Here is the generalized ATE of source sample on the target sample

``` r

ebal_ATE(xs,ys,trts,H_vars, target_moments, H_add_intercept=TRUE,delta=numeric(8)+0.1)$ATE
#>      value 
#> -0.6982687
```

### Compute approximate balancing confidence interval

Use resampling-based perturbation `RPM_AB()` with additional input of
`target_sd` and the number of bootstrap iteration 300.

``` r
## CI construction
target_sd = colStdevs(xt)[H_vars]
ATE_CI = RPM_AB(xs, ys,trts, H_vars=H_vars,target_mean=target_moments,target_sd=target_sd,num_sim=300, H_add_intercept=TRUE,cluster=5, set_seed=100)

## Lower bound of 95% CI
ATE_CI$lb_ATE
#>      2.5% 
#> -2.455437

## Upper bound of 95% CI
ATE_CI$ub_ATE
#>    97.5% 
#> 3.743462

## Number of simulations that uses exact balancing over 300
ATE_CI$use_exact
#> [1] 64
```

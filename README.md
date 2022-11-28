# Fit negative-binomial regression treating known/estimated exposure quantile functions as the covariate

In this document, we aim to demonstrate: (1) the use of the function to
fit negative-binomial (NB) regression models using known exposure
quantile functions as the covariate to estimate short-term effects of
environmental exposures in a time-series design, (2) the use of function
to estimate exposure quantile functions using previously proposed
semiparametric Bayesian approaches based on individual-level exposures,
(3) the use of to fit NB regression models using estimated exposure
quantile functions while accounting for uncertainties associated with
estimating quantile functions.

## The health model

The proposed NB regression models using exposure quantile functions as
the functional covariate is given below

$$
\\begin{align\*}
Y_i &\\sim NB(\\eta_i, \\xi) \\\\
E\[Y_i\] &= \\xi \\exp(\\eta_i) \\\\
\\eta_i &= \\int_0^1 \\beta(\\tau) Q_i(\\tau) + \\boldsymbol{\\gamma}^T\\boldsymbol{Z}\_i + \\boldsymbol{\\epsilon}\_i
\\end{align\*}
$$

## Example 1: Exposure quantile functions are independent across time points

Let’s simulate individual-level exposures and aggregate health outcomes
assuming data are collected over 1000 time points.

We begin with simulating individual-level exposures assuming
time-specific quantile functions are independent.

``` r
## the number of time points
num.time = 1000
## the number of individuals within each time point
num.size.time = 100
## simulate time-specific means and SDs of exposure distributions
mean.vec = rnorm(num.time, mean = 7.2, sd = 1)
sd.vec = rnorm(num.time, mean = 1, sd = 0.2)
## simulate individual-level exposures
x.sim.mat <- do.call(rbind, lapply(1:num.time, function(x) 
                                     rnorm(num.size.time, 
                                           mean.vec[x], sd.vec[x])))
# kbl(nobs.vec.print, booktabs = T) %>% 
#    kable_styling(latex_options = c("hold_position")) 
```

We then generate aggregate health outcomes by assuming *β*(*τ*) = *τ*.

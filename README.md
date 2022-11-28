# Fit negative-binomial regression treating known/estimated exposure quantile functions as the covariate

In this document, we aim to demonstrate: (1) the use of the function
`fit.health.knownquan` to fit negative-binomial (NB) regression models
using known exposure quantile functions as the covariate to estimate
short-term effects of environmental exposures in a time-series design,
(2) the use of function `fit.exposure` to estimate exposure quantile
functions using previously proposed semiparametric Bayesian approaches
based on individual-level exposures, (3) the use of
`fit.health.quan.errors` to fit NB regression models using estimated
exposure quantile functions while accounting for uncertainties
associated with estimating quantile functions.

## The health model

The proposed NB regression models using exposure quantile functions as
the functional covariate is given below

$$
\\begin{align\*}
&Y_i \\sim NB(\\eta_i, \\xi) , \\; \\; \\; E\[Y_i\] = \\xi \\exp(\\eta_i), \\\\
&\\eta_i = \\int_0^1 \\beta(\\tau) Q_i(\\tau) + \\boldsymbol{\\gamma}^T\\boldsymbol{Z}\_i + \\epsilon_i,
\\end{align\*}
$$
where *y*<sub>*i*</sub> denotes aggregated counts (i.e., the number
deaths) observed at the time point *i*, *Q*<sub>*i*</sub>(*τ*) is the
exposure quantile function of a continuous exposure at the time point
*i*, **Z**<sub>*i*</sub> is a vector of other covariates,
*ϵ*<sub>*i*</sub> represents a mean-zero residual process, and *ξ* is
the over-dispersion parameter.

## Example 1: Exposure quantile functions are independent across time points

### Simulate data

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
```

We then generate aggregate health outcomes by assuming
*η*<sub>*i*</sub> = *β*<sub>0</sub> + ∫<sub>0</sub><sup>1</sup>*β*(*τ*)*Q*<sub>*i*</sub>(*τ*)
and *β*(*τ*) = *τ*.

``` r
## specify the beta(tau) 
betatau.true <- function(x) {x}
## simulate aggregate health outcomes 
integral <- function(x, betatau.fun, mean, sd){
  re = betatau.fun(x)*qnorm(x, mean = mean, sd = sd)
  return(re)
}
integration.vec <- do.call(rbind,
                           lapply(1:num.time,
                                  function(y)
                                    integrate(integral, lower = 0, upper = 1,
                                              betatau.fun = betatau.true,
                                              mean = mean.vec[y],
                                              sd = sd.vec[y])$value))

## intercept in the health model
beta0 = -3.5
## the over-dispersion parameter 
xi = 20 
eta.vec.true = beta0 + integration.vec
q.nb = 1/(1+exp(eta.vec.true))
y.vec <- apply(cbind(rep(xi, num.time), q.nb), 1,
               function(x) rnbinom(1, size = x[1], prob = x[2]))
```

### Fit the health model assuming exposure quantile functions are known

The *β*(*τ*) is modeled via the basis expansion. Specifically, *β*(*τ*)
is expanded using orthonormal Bernstein polynomials of degree *n*, in
this example we set *n* = 1. To approximate *β*(*τ*) = *τ* using
orthonormal Bernstein polynomials of degree 1, true values of basis
coefficients are (0.288,0.500).

``` r
re.fit.knownquan <- fit.health.knownquan(y = y.vec,
                                         Z.design = NULL,
                                         n = 1, dist.type = "normal",
                                         mean = mean.vec, sd = sd.vec,
                                         rand.int = FALSE,
                                         niter = 5000, burn_in = 2500)
```

Check trace plots of the intercept, basis coefficients, and the
over-dispersion parameter. ![This is an
image](./results/trace_knownquan.png){height = 50%}

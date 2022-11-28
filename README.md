README
================

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

![
\\begin{align\*}
&Y_i \\sim NB(\\eta_i, \\xi) , \\; \\; \\; E\[Y_i\] = \\xi \\exp(\\eta_i), \\\\
&\\eta_i = \\int_0^1 \\beta(\\tau) Q_i(\\tau) + \\boldsymbol{\\gamma}^T\\boldsymbol{Z}\_i + \\epsilon_i,
\\end{align\*}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Balign%2A%7D%0A%26Y_i%20%5Csim%20NB%28%5Ceta_i%2C%20%5Cxi%29%20%2C%20%5C%3B%20%5C%3B%20%5C%3B%20E%5BY_i%5D%20%3D%20%5Cxi%20%5Cexp%28%5Ceta_i%29%2C%20%5C%5C%0A%26%5Ceta_i%20%3D%20%5Cint_0%5E1%20%5Cbeta%28%5Ctau%29%20Q_i%28%5Ctau%29%20%2B%20%5Cboldsymbol%7B%5Cgamma%7D%5ET%5Cboldsymbol%7BZ%7D_i%20%2B%20%5Cepsilon_i%2C%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
&Y_i \sim NB(\eta_i, \xi) , \; \; \; E[Y_i] = \xi \exp(\eta_i), \\
&\eta_i = \int_0^1 \beta(\tau) Q_i(\tau) + \boldsymbol{\gamma}^T\boldsymbol{Z}_i + \epsilon_i,
\end{align*}
")

where ![y_i](https://latex.codecogs.com/png.latex?y_i "y_i") denotes
aggregated counts (i.e., the number deaths) observed at the time point
![i](https://latex.codecogs.com/png.latex?i "i"),
![Q_i(\\tau)](https://latex.codecogs.com/png.latex?Q_i%28%5Ctau%29 "Q_i(\tau)")
is the exposure quantile function of a continuous exposure at the time
point ![i](https://latex.codecogs.com/png.latex?i "i"),
![\\boldsymbol{Z}\_i](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7BZ%7D_i "\boldsymbol{Z}_i")
is a vector of other covariates,
![\\epsilon_i](https://latex.codecogs.com/png.latex?%5Cepsilon_i "\epsilon_i")
represents a mean-zero residual process, and
![\\xi](https://latex.codecogs.com/png.latex?%5Cxi "\xi") is the
over-dispersion parameter.

## Example 1: Exposure quantile functions are independent across time points

### Simulate data

Let us simulate individual-level exposures and aggregate health outcomes
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
![\\eta_i = \\beta_0 + \\int_0^1 \\beta(\\tau) Q_i(\\tau)](https://latex.codecogs.com/png.latex?%5Ceta_i%20%3D%20%5Cbeta_0%20%2B%20%5Cint_0%5E1%20%5Cbeta%28%5Ctau%29%20Q_i%28%5Ctau%29 "\eta_i = \beta_0 + \int_0^1 \beta(\tau) Q_i(\tau)")
and
![\\beta(\\tau) = \\tau](https://latex.codecogs.com/png.latex?%5Cbeta%28%5Ctau%29%20%3D%20%5Ctau "\beta(\tau) = \tau").

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

The
![\\beta(\\tau)](https://latex.codecogs.com/png.latex?%5Cbeta%28%5Ctau%29 "\beta(\tau)")
is modeled via the basis expansion. Specifically,
![\\beta(\\tau)](https://latex.codecogs.com/png.latex?%5Cbeta%28%5Ctau%29 "\beta(\tau)")
is expanded using orthonormal Bernstein polynomials of degree
![n](https://latex.codecogs.com/png.latex?n "n"), in this example we set
![n = 1](https://latex.codecogs.com/png.latex?n%20%3D%201 "n = 1"). To
approximate
![\\beta(\\tau) = \\tau](https://latex.codecogs.com/png.latex?%5Cbeta%28%5Ctau%29%20%3D%20%5Ctau "\beta(\tau) = \tau")
using orthonormal Bernstein polynomials of degree
![1](https://latex.codecogs.com/png.latex?1 "1"), true values of basis
coefficients are
![(0.288, 0.500)](https://latex.codecogs.com/png.latex?%280.288%2C%200.500%29 "(0.288, 0.500)").

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
image](./results/trace_knownquan.png)

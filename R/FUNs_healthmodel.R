############################################################################
### Functions used in fitting health model with known/unknown quantile funs
### > aggregate health outcomes are assumed for follow NB dist
### > beta(tau) is expanded using orthonormal Bernstein polynomials
### > Gibbs sampling is implemented for updating regression coefs by 
### introducing Polya-Gamma latent random variables 
### yuzi.zhang@emory.edu
############################################################################

### Reference about Bernstein polynomials ###
# Bellucci, M. A. (2014). 
# On the explicit representation of orthonormal Bernstein polynomials. 
# arXiv preprint arXiv:1404.2293.

### orthonormal Bernstein basis polynomials of degree n ###
### Eqn. (7)
bernstein_orth <- function(tau, j, n){
  comp1 = sqrt(2*(n-j)+1)
  comp2 = (1-tau)^(n-j)
  k.vec = 0:j
  comp3 <- function(x){
    sum((-1)^(k.vec)*choose(2*n+1-k.vec, j-k.vec)*choose(j, k.vec)*x^(j-k.vec))
  } 
  comp4  = sapply(tau, comp3, simplify = T)
  return(comp1*comp2*comp4)
}


## A function to compute int_beta(tau)dtau 
## based on orthonormal Bernstein basis polynomials ##
get_int_beta_bernstein <- function(n, coef){
  ## INPUTS:
  ## n: degrees of orthonormal Bernstein basis functions
  ## coef: basis coefficients
  ## OUTPUT: a scalar int_beta(tau)dtau 
  intergral <- function(x, n, coef){
    X <- NULL
    for(j in 0:n){
      X <- cbind(X, bernstein_orth(tau = x, j = j, n = n))
    }
    return( X%*%matrix(coef, ncol = 1) )
  }
  re <- integrate(intergral, lower = 0, upper = 1,
                  n = n, coef = coef)$value
  return(re)
}


## A function to compute beta(tau) 
## based on orthonormal Bernstein polynomials of degree n ##
beta.tau.bernstein <- function(tau, n, beta.vec){
  ## INPUTS:
  ## tau: the percentile level that beta(tau) is evaluated at
  ## n: degrees of orthonormal Bernstein basis functions
  ## beta.vec: basis coefficients
  ## OUTPUT: a scalar beta(tau)
  X <- NULL
  for(j in 0:n){
    X <- cbind(X, bernstein_orth(tau = tau, j = j, n = n))
  }
  return(X%*%matrix(beta.vec, ncol = 1))
}


## A function to compute "exposure covariate" 
## when exposure quantile functions are known
## i.e., compute integration of int_K(.,tau)Q(tau)dtau
get_Xdesign_knownquan <- function(n, dist.type, 
                                  mean = NULL, sd = NULL,
                                  meanlog = NULL, sdlog = NULL,
                                  shape = NULL, rate = NULL, ...){
  ## INPUTS:
  ## n: degrees of orthonormal Bernstein basis functions
  ## dist.type: a character specifies the distribution assumptions 
  ## that you are willing to make for the exposure of interest
  ## three options are available:
  ## "normal", "lognormal", "gamma"
  ## mean/sd: a vector of length that equals to the number of time points, 
  ## this vector contains the group-specific mean/standard deviation of 
  ## the assumed normal dist
  ## logmean/logsd: 
  ## this vector contains the group-specific mean/sd on the log scale of 
  ## the assumed log-normal dist
  ## shape/rate: 
  ## this vector contains the group-specific shape/rate of the assumed gamma dist
  ## OUTPUT:
  ## a design matrix of dimension (N x (n+1)), 
  ## N is the number of observations
  ## n+1 is the number of basis functions for expanding beta(tau) 
  
  if(! dist.type %in% c("normal", "lognormal", "gamma")){
    stop("The exposure distribution can only be assumed as normal/lognormal/gamma")
  }else if(dist.type == "normal" & (is.null(mean) | is.null(sd))){
    stop("Please provide means and sds for assumed normal distributions using arguments mean and sd")
  }else if(dist.type == "lognormal" & (is.null(meanlog) | is.null(sdlog))){
    stop("Please provide means and sds on log scale for assumed log-normal distributions using arguments meanlog and sdlog")
  }else if(dist.type == "gamma" & (is.null(shape) | is.null(rate))){
    stop("Please provide shapes and rates for assumed gamma distributions using arguments shape and rate")
  }
  
  if(dist.type == "normal" & !is.null(mean) & !is.null(sd)){
    K = n+1; num.obs = length(mean)
    get_xvec <- function(K, mean.i, sd.i){
      integration.vec <- rep(NA, K)
      for(k in 0:(K-1)){
        intergral <- function(x, df, mean.i, sd.i){
          bernstein_orth(tau = x, j = k, n = df)*
            qnorm(p = x, mean = mean.i, sd = sd.i)
        }
        integration.vec[k+1] <- integrate(intergral, 
                                          lower = 0, upper = 1,
                                          df = K-1, 
                                          mean.i = mean.i, sd.i = sd.i)$value
      }
      return(integration.vec)
    }
    X.design <- do.call(rbind, lapply(1:num.obs, function(x) 
      get_xvec(K = K, mean.i = mean[x], sd.i = sd[x])))
    return(X.design)
    
  }else if(dist.type == "lognormal" & 
           !is.null(meanlog) & !is.null(sdlog)){
    K = n+1; num.obs = length(meanlog)
    get_xvec <- function(K, meanlog.i, sdlog.i){
      integration.vec <- rep(NA, K)
      for(k in 0:(K-1)){
        intergral <- function(x, df, meanlog.i, sdlog.i){
          bernstein_orth(tau = x, j = k, n = df)*
            qlnorm(p = x, meanlog = meanlog.i, sdlog = sdlog.i)
        }
        integration.vec[k+1] <- integrate(intergral, 
                                          lower = 0, upper = 1,
                                          df = K-1, 
                                          meanlog.i = meanlog.i, 
                                          sdlog.i = sdlog.i)$value
      }
      return(integration.vec)
    }
    X.design <- do.call(rbind, lapply(1:num.obs, function(x) 
      get_xvec(K = K, meanlog.i = meanlog[x], sdlog.i = sdlog[x])))
    return(X.design)
    
  }else if(dist.type == "gamma" & 
           !is.null(shape) & !is.null(rate)){
    K = n+1; num.obs = length(shape)
    get_xvec <- function(K, shape.i, rate.i){
      integration.vec <- rep(NA, K)
      for(k in 0:(K-1)){
        intergral <- function(x, df, shape.i, rate.i){
          bernstein_orth(tau = x, j = k, n = df)*
            qgamma(p = x, shape = shape.i, rate = rate.i)
        }
        integration.vec[k+1] <- integrate(intergral, 
                                          lower = 0, upper = 1,
                                          df = K-1, 
                                          shape.i = shape.i, 
                                          rate.i = rate.i)$value
      }
      return(integration.vec)
    }
    X.design <- do.call(rbind, lapply(1:num.obs, function(x) 
      get_xvec(K = K, shape.i = shape[x], rate.i = rate[x])))
    return(X.design)
  }
  
}

## A function to get basis functions that used for expanding quantile functions
## Basis functions given in Equations (3) and (4) in Reich (2012) ##
Bl <- function(tau.vec, l, L, basis.fun, shape = 5){
  k.vec = seq(0, 1, length.out = L+1)
  if(basis.fun == "Gau"){
    B.vec <- c()
    if(k.vec[l] < 0.5){
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = qnorm(k.vec[l]) - qnorm(k.vec[l+1]) 
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qnorm(tau) - qnorm(k.vec[l+1]) 
        }else{
          B = 0
        }
        B.vec <- c(B.vec, B)
      }
    }else{
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = 0
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qnorm(tau) - qnorm(k.vec[l]) 
        }else{
          B = qnorm(k.vec[l+1]) - qnorm(k.vec[l])
        }
        B.vec <- c(B.vec, B)
      }
    }
    return(B.vec)
  }
  if(basis.fun == "Gamma"){
    B.vec <- c()
    if(k.vec[l] < 0.5){
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = qgamma(k.vec[l], shape = shape) - qgamma(k.vec[l+1], shape = shape) 
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qgamma(tau, shape = shape) - qgamma(k.vec[l+1], shape = shape) 
        }else{
          B = 0
        }
        B.vec <- c(B.vec, B)
      }
    }else{
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = 0
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qgamma(tau, shape = shape) - qgamma(k.vec[l], shape = shape) 
        }else{
          B = qgamma(k.vec[l+1], shape = shape) - qgamma(k.vec[l], shape = shape)
        }
        B.vec <- c(B.vec, B)
      }
    }
    return(B.vec)
  }
  
}


### A function to compute integration of K(tau)B(tau)^T 
### in Eqn. (7) in our paper ###
get_integration_bernstein <- function(L, n, basis.fun, shape = 5){
  ## INPUTS:
  ## L: number of basis functions to expand exposure quantile functions
  ## n: degrees of orthonormal Bernstein basis functions
  ## basis.fun: name of the basis function used for estimating quan
  ## shape: shape parameter for gamma basis fun used for estimating quan
  ## OUTPUT: (L+1)xK (B(tau)K(.,tau)^T) matrix 
  K = n+1
  integration.mat <- matrix(NA, ncol = K, nrow = L + 1)
  for(ll in 1:(L+1)){
    if(ll == 1){
      for(k in 0:(K-1)){
        intergral <- function(x, df){
          bernstein_orth(tau = x, j = k, n = df)
        }
        integration.mat[ll, k+1] <- integrate(intergral, lower = 0, upper = 1,
                                              df = n)$value
      }
    }else{
      for(k in 0:(K-1)){
        intergral <- function(x, l, L, df){
          Bl(tau.vec = x, l = l, L = L, 
             basis.fun = basis.fun,
             shape = shape)*bernstein_orth(tau = x, j = k, n = df)
        }
        integration.mat[ll, k+1] <- integrate(intergral, lower = 0, upper = 1,
                                              l = ll-1, L = L,
                                              df = n)$value
      }
    }
  }
  return(integration.mat)
}

#' Compute basis coefficients of orthonormal Bernstein polynomials used for
#' approximating a given function
#'
#' This function computes basis coefficients of orthonormal Bernstein polynomials used for
#' approximating a given function following Eqn. (22) provided in Bellucci (2014). 
#' Please refer to Section 5 Function Approcimation in Bellucci (2014) for more details.  
#'
#' @param j a number specifying the index of the basis coefficient, \eqn{j=0,\dots,n}. 
#' @param n a number specifying the degrees of orthonormal Bernstein polynomials. 
#' @param f the name of self-defined R function that to be approximated. 
#' 
#' @references 
#' Bellucci, M. A. (2014). On the explicit representation of orthonormal Bernstein polynomials. 
#' \emph{arXiv preprint} \href{https://arxiv.org/abs/1404.2293}{\emph{arXiv:1404.2293}}.
#' 
#' @examples
#' # Compute basis coefficients of orthonormal Bernstein polynomials used 
#' # for approximating the linear function f(x) = x
#' beta.true <- sapply(0:1, function(x) get_beta_bernstein(x, n = 1, f = function(y){y}))
#' 
#' @export
## A function to compute coefs of orthonormal Bernstein basis functions 
## used for approximating a given function 
## Eqn. (22)
get_beta_bernstein <- function(j, n, f){
  ## INPUTS:
  ## n: degree of orthonormal Bernstein basis functions
  ## f: function to be approximated 
  integral <- function(x){
    bernstein_orth(x, j, n)*f(x)
  }
  re = integrate(integral, lower = 0, upper = 1)$value
  return(re)
}


#' Fit negative-binomial regression using known exposure quantile functions as the covariate 
#'
#' This function fits negative binomial regression models using known 
#' exposure quantile functions as the covariate to analyze 
#' aggregate health outcomes in environmental epidemiological studies with a time-series design.
#'
#' @param y a vector containing observed aggregate health outcomes.
#' @param Z.design a design matrix containing confounders that are controlled for.
#' @param n a number specifying the degrees of orthonormal Bernstein basis functions 
#' that are used for expanding \eqn{\beta(\tau)}.
#' @param dist.type a character specifying the distribution assumptions imposed 
#' for within-unit individual-level exposures. 
#' This function accepts three distribution assumptions 
#' \code{normal} (normal), \code{lognormal} (log-normal distribution), 
#' and \code{gamma} (gamma distribution).
#' @param mean If \code{dist.type} = \code{normal}, 
#' a vector of length that equals the number of units 
#' (e.g., the number of time points in a time-series design). 
#' This vector contains unit-specific means of assumed normal distributions. 
#' @param sd If \code{dist.type} = \code{normal}, 
#' a vector containing unit-specific standard deviations of assumed normal distributions.
#' @param meanlog If \code{dist.type} = \code{lognormal}, 
#' a vector containing unit-specific means on the log scale of assumed log-normal distributions.
#' @param sdlog If \code{dist.type} = \code{lognormal}, 
#' a vector containing unit-specific standard deviations on the log scale of assumed log-normal distributions.
#' @param shape If \code{dist.type} = \code{gamma}, 
#' a vector containing unit-specific shape parameters of assumed gamma distributions.
#' @param rate If \code{dist.type} = \code{gamma}, 
#' a vector containing unit-specific rate parameters of assumed gamma distributions.
#' @param rand.int a logical value indicating whether unit-specific random intercepts are included.
#' @param str a character string specifying the structure of unit-specific random intercepts 
#' when \code{rand.int} = \code{TRUE}. \code{AR-1} fits random walk model of order 1 
#' and \code{exch} specifies exchangeable unit-specific random intercepts.
#' @param niter a number specifying the number of MCMC iterations, the default is 5000.
#' @param burn_in a number specifying the number of iterations that is discarded, the default is 2500. 
#' @param returnpost a logical value indicating whether the returned MCMC samples containing burn-in samples.
#' 
#' @return This function returns a list containing: (1) a matrix containing samples of 
#' basis coefficients of orthonormal Bernstein polynomials and coefficients of controlled confounders, 
#' (2) WAIC, and (3) parameters related to the mean-zero residual process (if \code{rand.int} = \code{TRUE}). 
#'
#' @references 
#' Bellucci, M. A. (2014). On the explicit representation of orthonormal Bernstein polynomials. 
#' \emph{arXiv preprint} \href{https://arxiv.org/abs/1404.2293}{\emph{arXiv:1404.2293}}.
#' 
#' @examples
#' # Simulate time-specific means and SDs of exposure quantile functions
#' mean.vec <- rnorm(1000, mean = 7.2, sd = 1)
#' sd.vec <- rnorm(1000, mean = 1, sd = 0.2)
#' 
#' # Specify the true beta(tau)
#' betatau.true <- function(x) {x}
#' 
#' # Generate aggregate health outcomes 
#' integral <- function(x, betatau.fun, mean, sd){
#'   re = betatau.fun(x)*qnorm(x, mean, sd)
#'   return(re)
#' }
#' 
#' integration.vec <- do.call(rbind,
#' lapply(1:1000, function(y) integrate(integral, lower = 0, upper = 1,
#' betatau.fun = betatau.true, mean = mean.vec[y], sd = sd.vec[y])$value))
#' 
#' y.vec <- apply(cbind(rep(20, 1000), 1/(1+exp(-3.5 + integration.vec))), 1,
#' function(x) rnbinom(1, size = x[1], prob = x[2]))
#' 
#' fit.health.knownquan(y = y.vec, n = 1, dist.type = "normal", 
#' mean = mean.vec, sd = sd.vec, rand.int = FALSE, 
#' niter = 5000, burn_in = 2500)
#' 
#' @export
## A function to fit health model 
## using orthonormal Bernstein polynomials basis to model beta(tau)
### Yt ~ NB(eta_t, xi); 
### E[Yt] = mu_t = xi*exp(eta_t); Var[Yt] = mu_t + 1/xi*mu_t^2
### eta_t = beta0 + \int_Qt(tau)beta(tau)dtau + gammaW
### beta(tau) = \sum Bern(j, n)beta_j; basis functions representation
### orthonormal Bernstein basis polynomials: Eqn. (7)
### Qt(tau) = \sum B(tau)theta_t; basis functions representation
fit.health.knownquan <- function(y, 
                                 Z.design = NULL,
                                 n, dist.type, 
                                 mean = NULL, sd = NULL,
                                 meanlog = NULL, sdlog = NULL,
                                 shape = NULL, rate = NULL,
                                 rand.int = FALSE,
                                 str = NULL,
                                 niter = 5000, 
                                 burn_in = 2500,
                                 returnpost = FALSE, ...){
  ## INPUTs:
  ## y: a vector contains observed aggregate health outcomes 
  ## Z.design: a design matrix containing confounders controlled  
  ## n: the degrees of orthonormal Bernstein basis functions
  ## dist.type: a character specifies the distribution assumptions 
  ## that you are willing to make for the exposure of interest
  ## three options are available:
  ## "normal", "lognormal", "gamma"
  ## mean/sd: a vector of length that equals to the number of time points, 
  ## this vector contains the group-specific mean/standard deviation of 
  ## the assumed normal dist
  ## logmean/logsd: 
  ## this vector contains the group-specific mean/sd on the log scale of 
  ## the assumed log-normal dist
  ## shape/rate:
  ## this vector contains the group-specific shape/rate of the assumed gamma dist

  ## rand.int: 
  ## a logical value indicating whether time-specific intercepts are imposed
  ## str: a character string specifying the structure of the random intercept
  ## "AR-1" (1st order Gaussian Markov random field process), 
  ## "exch" (exchangeable)
  
  ## burn_in: # of burn-in samples (default is 2500)
  ## niter: # of MCMC iterations
  ## returnpost: a logical value indicating whether 
  ## the returned MCMC samples containing burn-in samples
  
  # ----------- PREPARE DATA ------------ #
  ## get design matrix with known quantile functions 
  X.design <- get_Xdesign_knownquan(n = n, dist.type = dist.type, 
                                    mean = mean, sd = sd,
                                    meanlog = meanlog, sdlog = sdlog,
                                    shape = shape, rate = rate)
  X.design <- cbind(1, X.design, Z.design)
  dat.sim = data.frame(y = y, X.design)
  nbeta.input = ncol(X.design)
  colnames(dat.sim) <- c("y", "Int", paste0("x", 1:(nbeta.input-1)))
  
  # ------------ PRIORS ----------- #
  # priors #
  ## beta: basis coefs for expanding beta(tau) 
  c.beta = matrix(rep(0, nbeta.input), ncol = 1)
  C.beta.pre = diag(rep(1/100, nbeta.input))
  
  # tuning parameters for the over-dispersion parameter #
  xi.tun = 0.2
  
  ## initialize data frames to store posterior samples ##
  ntimes = length(y)
  omega.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.input, nrow = niter + 1)
  colnames(beta.mat) <- paste0("beta", 0:(nbeta.input-1))
  beta.mat[1, ] <- rep(0, nbeta.input)
  
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  colnames(accept.mat) <- c("xi")
  
  if(! rand.int){
    parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
    colnames(parm.mat) <- c("xi")
    parm.mat[1, ] <- c(5)
  }else if(rand.int & str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2", "xi")
    parm.mat[1, ] <- c(0.1, 5)
    ## data frame to store random effects
    theta.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
    theta.mat[,1] <- 0
    ## the InvGamma prior for the variance of temporal random effects 
    a = 0.1; b = 0.1
  }else if(rand.int & str == "AR-1"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2", "rho", "xi")
    parm.mat[1, ] <- c(0.1, 0.9, 5)
    ## data frame to store random effects
    theta.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
    theta.mat[,1] <- 0
    ## adjacency matrix ##
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))
    ## the InvGamma prior for variance of temporal random effects 
    a = 0.1; b = 0.1
    ## the discrete prior for rho controlling temporal dependency 
    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)),
                    simplify = TRUE)
  }

  pg.parm1 <- dat.sim$y + parm.mat[1, "xi"]
  pg.parm2 <- as.matrix(dat.sim[,c("Int", paste0("x", 1:(nbeta.input-1)))])%*%
    matrix(beta.mat[1, ], ncol = 1)
  ntotal = ntimes
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    ## compute latent normal variables ##
    z.vec <- (dat.sim$y - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2. update  beta (conjugate)
    ## posterior mean and variance of beta ##
    if(! rand.int){
      M <- solve(crossprod(X.design*sqrt(omega.vec)) + C.beta.pre)
      m <- M%*%(C.beta.pre%*%c.beta + 
                  t(sqrt(omega.vec)*X.design)%*%
                  (sqrt(omega.vec)*matrix(z.vec, ncol = 1)))
      beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    }else{
      res.vec = z.vec - theta.mat[i,]
      M <- solve(crossprod(X.design*sqrt(omega.vec)) + C.beta.pre)
      m <- M%*%(C.beta.pre%*%c.beta + 
                  t(sqrt(omega.vec)*X.design)%*%
                  (sqrt(omega.vec)*matrix(res.vec, ncol = 1)))
      beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    }

    
    # update random effects (conjugate)
    ## posterior mean and variance of random effects ##
    if(rand.int){
      if(str == "AR-1"){
        pre.theta.0 <- (1/parm.mat[i, "tau2"])*(Dw.t - parm.mat[i, "rho"]*W.t)
        pre.theta.post <- Omega + pre.theta.0
        Sigma.theta.post <- solve(pre.theta.post)
        fix.eff <- X.design%*%matrix(beta.mat[i+1, ], ncol = 1)
        mu.theta.post <- Sigma.theta.post%*%(Omega%*%matrix(z.vec - fix.eff, ncol = 1))
        theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.theta.post,
                                                Sigma = Sigma.theta.post) )
        ## update tau2 
        parm.mat[i+1, "tau2"] <- 1/rgamma(1, ntimes/2 + a,
                                          b + as.numeric(sum(theta.mat[i+1, ]^2))/2)
        
      }else if(str == "exch"){
        pre.theta.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2"]), ntimes))
        pre.theta.post <- Omega + pre.theta.0
        Sigma.theta.post <- solve(pre.theta.post)
        fix.eff <- X.design%*%matrix(beta.mat[i+1, ], ncol = 1)
        mu.theta.post <- Sigma.theta.post%*%(Omega%*%matrix(z.vec - fix.eff, ncol = 1))
        theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.theta.post,
                                                Sigma = Sigma.theta.post) )
        ## update tau2 and rho
        theta.unique <- matrix(theta.mat[i+1, ], ncol = 1)
        parm.mat[i+1, "tau2"] <- 1/rgamma(1, ntimes/2 + a1,
                                          b1 + as.numeric(t(theta.unique)%*%
                                                            (Dw.t - parm.mat[i, "rho"]*W.t)%*%theta.unique)/2)
        inter = as.numeric( t(matrix(theta.unique, ncol = 1))%*%W.t%*%
                              matrix(theta.unique, ncol = 1) )
        ll.rho = ll.rho.1 + rho.prior.val/(2*parm.mat[i+1, "tau2"]^2)*inter
        parm.mat[i+1, "rho"] <- sample(x = rho.prior.val, size = 1,
                                       prob = exp(ll.rho - max(ll.rho)))
        
      }
    }
    
    # 3. update xi (the over-dispersion parameter) 
    eta.vec.i1 = X.design%*%matrix(beta.mat[i+1, ], ncol = 1)
    q.nb = 1/(1+exp(eta.vec.i1)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(dat.sim$y, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(dat.sim$y, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    omega.mat[i+1, ] <- rpg(ntotal, dat.sim$y + parm.mat[i+1, "xi"], eta.vec.i1)
    
    ## tuning parameters ##
    if(i <= 5000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
    
  }# END MCMC iterations 
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  compq <- function(x, X.design){
    q = 1/(1+exp(X.design%*%matrix(x[1:(nbeta.input)], ncol = 1)))
  }
  q.vec = apply(beta.mat[-(1:burn_in), ], 1, compq, X.design = X.design)
  
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(dat.sim$y, q.vec), 1, compll, 
                    xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  if(returnpost){
    re <- list(beta = beta.mat[-c(1:burn_in), ], 
               parm = parm.mat[-c(1:burn_in), ], 
               accept.rate = mean(accept.mat[-c(1:burn_in), ]), 
               WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
  }
  if(! returnpost){
    re <- list(beta = beta.mat, parm = parm.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), ]), 
               WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
  }
  ## OUTPUTS
  ## A list containing 
  ## if the argument is true: 
  ## (1) a matrix containing MCMC samples with burn-in samples discarded 
  ##     (basis coefs + coefs of confounders) 
  ## (2) WAIC 
  ## if the argument is false: 
  ## (1) a matrix containing MCMC samples (basis coefs + coefs of confounders) 
  ## (2) WAIC
  return(re)
}





#' Fit negative-binomial regression using estimated exposure quantile functions as the covariate 
#'
#' This function fits negative binomial regression models using estimated 
#' exposure quantile functions as the covariate while accounting for uncertainties associated 
#' with estimating exposure quantile functions by introducing MVN prior for basis coefficients 
#' used for modeling exposure quantile functions to analyze 
#' aggregate health outcomes in environmental epidemiological studies with a time-series design.
#'
#' @param y a vector containing observed aggregate health outcomes.
#' @param Z.design a design matrix containing confounders that are controlled for.
#' @param L a number specifying the number of of basis functions that has been used in the 
#' basis expansion for modeling exposure quantile functions.
#' @param basis.fun a character string specifying which basis functions that have been used 
#' for modeling quantile processes. This function accepts \code{Gau} to specify 
#' piecewise Gaussian functions and \code{Gamma} to specify piecewise Gamma functions.
#' @param basis.shape a number specifying the shape parameter of the 
#' piecewise Gamma functions if \code{basis.fun = Gamma}, the default value is 5.
#' @param n a number specifying the degrees of orthonormal Bernstein basis functions 
#' that are used for expanding \eqn{\beta(\tau)}.
#' @param theta.pri.mu a list of length \eqn{T} (i.e., the number of time points) 
#' with each element including the mean of the multivariate normal (MVN) 
#' prior assumed for basis coefficients that has been used to model exposure quantile functions.
#' @param theta.pri.pre a list with each element including 
#' the precision matrix of the MVN prior. 
#' @param str a character string specifying the structure of unit-specific random intercepts 
#' when \code{rand.int} = \code{TRUE}. \code{AR-1} fits random walk model of order 1 
#' and \code{exch} specifies exchangeable unit-specific random intercepts.
#' @param niter a number specifying the number of MCMC iterations, the default is 5000.
#' @param burn_in a number specifying the number of iterations that is discarded, the default is 2500. 
#' 
#' @return This function returns a list containing: (1) a matrix containing samples of 
#' basis coefficients of orthonormal Bernstein polynomials and coefficients of controlled confounders, 
#' (2) WAIC, and (3) parameters related to the mean-zero residual process (if \code{rand.int} = \code{TRUE}). 
#'
#' @references 
#' Bellucci, M. A. (2014). On the explicit representation of orthonormal Bernstein polynomials. 
#' \emph{arXiv preprint} \href{https://arxiv.org/abs/1404.2293}{\emph{arXiv:1404.2293}}.
#' 
#' @examples
#' # Simulate time-specific means and SDs of exposure quantile functions
#' mean.vec <- rnorm(10, mean = 7.2, sd = 1)
#' sd.vec <- rnorm(10, mean = 1, sd = 0.2)
#' 
#' # Specify the true beta(tau)
#' betatau.true <- function(x) {x}
#' 
#' # Generate aggregate health outcomes 
#' integral <- function(x, betatau.fun, mean, sd){
#'   re = betatau.fun(x)*qnorm(x, mean, sd)
#'   return(re)
#' }
#' 
#' integration.vec <- do.call(rbind,
#' lapply(1:10, function(y) integrate(integral, lower = 0, upper = 1,
#' betatau.fun = betatau.true, mean = mean.vec[y], sd = sd.vec[y])$value))
#' 
#' y.vec <- apply(cbind(rep(20, 10), 1/(1+exp(-3.5 + integration.vec))), 1,
#' function(x) rnbinom(1, size = x[1], prob = x[2]))
#' 
#' # Estimate exposure quantile functions assuming exposure quantile functions 
#' # are independent across time points.
#' re.fit.exp <- NULL
#' for(i in 1:10){
#' re.fit.exp[[i]] <- fit.exposure(x.ind = x.sim.mat[i, ],
#' basis.fun = "Gau", L = 4, niter = 10000, inde = T) }
#' 
#' # Obtain means and precision matrices of the MVN prior.
#' MVNprior.mu.pre.ind <- exposure.prior(re = re.fit.exp, burn_in = 5000, inde = T)
#' 
#' re.fit.quan.errors <- fit.health.quan.errors(y = y.vec, 
#' L = 4, basis.fun = "Gau", n = 1, 
#' theta.pri.mu = MVNprior.mu.pre.ind$theta.pri.mu, 
#' theta.pri.pre = MVNprior.mu.pre.ind$theta.pri.pre,
#' rand.int = FALSE,
#' niter = 5000, burn_in = 2500)
#' 
#' @export
##### update beta as a whole vector #####
##### introduce a "prior" for alpha_t,0 and theta_t,l ######
fit.health.quan.errors <- function(y, 
                                   Z.design = NULL,
                                   L, basis.fun, basis.shape, n, 
                                   theta.pri.mu, theta.pri.pre,
                                   rand.int = FALSE,
                                   str = NULL,
                                   niter = 5000, burn_in = 2500){
  
  ## INPUTs:
  ## y: a vector contains observed aggregate health outcomes 
  ## Z.design: a design matrix containing confounders controlled for
  ## L: the number of basis functions has been used to model quantile functions
  ## basis.fun: a character string specifying what basis functions have been used
  ## for expanding quantile functions 
  ## "Gau" = Gaussian piecewise functions
  ## "Gamma" = Gamma piecewise functions
  ## if basis.fun = "Gamma", basis.shape is the shape parameter specified for
  ## the Gamma piecewise functions; default value is 5
  ## n: the degrees of orthonormal Bernstein basis functions for modeling beta(tau)
  ## theta.pri.mu: 
  ## a list contains means of MVN priors introduced for 
  ## basis coefs used for modeling quantile processes
  ## theta.pri.pre: 
  ## a list contains precision matrix of MVN priors  (alpha0, theta_1, ... ,theta_L) 
  ## rand.int: 
  ## a logical value indicating whether time-specific intercepts are imposed
  ## str: a character string specifying the structure of the random intercept
  ## "AR-1" (1st order Gaussian Markov random field process), 
  ## "exch" (exchangeable)
  ## niter: the number of MCMC iterations 
  ## burn_in: the number of burn-in samples default is 2500
  
  
  # -------- compute BK.mat: int B(tau)phi(tau) dtau (L+1)x(K) ------- #
  BK.mat = get_integration_bernstein(L = L, n = n,
                                     basis.fun = basis.fun,
                                     shape = basis.shape)
  
  # ------------ priors ----------- #
  # specify priors #
  if(is.null(Z.design)){
    num.beta.x = ncol(BK.mat)
    num.beta = num.beta.x + 1
    c.beta.x = matrix(rep(0, num.beta.x), ncol = 1)
    C.beta.pre.x = diag(rep(1/100, num.beta.x))
    c.beta = matrix(rep(0, num.beta), ncol = 1)
    C.beta.pre = diag(rep(1/100, num.beta))
    
    beta.mat <- matrix(NA, ncol = num.beta, nrow = niter + 1)
    colnames(beta.mat) <- c("beta0", 
                            paste0("beta_x", 1:(num.beta.x)))
    beta.mat[1, ] <- rep(0, num.beta)
    
  }else{
    num.beta.z = ncol(Z.design)
    num.beta.x = ncol(BK.mat)
    num.beta = num.beta.z + num.beta.x + 1
    c.beta.z = matrix(rep(0, num.beta.z), ncol = 1)
    C.beta.pre.z = diag(rep(1/100, num.beta.z))
    c.beta.x = matrix(rep(0, num.beta.x), ncol = 1)
    C.beta.pre.x = diag(rep(1/100, num.beta.x))
    c.beta = matrix(rep(0, num.beta), ncol = 1)
    C.beta.pre = diag(rep(1/100, num.beta))
    
    beta.mat <- matrix(NA, ncol = num.beta, nrow = niter + 1)
    colnames(beta.mat) <- c("beta0", 
                            paste0("beta_z", 1:(num.beta.z)),
                            paste0("beta_x", 1:(num.beta.x)))
    beta.mat[1, ] <- rep(0, num.beta)
    
  }
  
  ntimes = length(y)
  theta.pri.mu.mat <- do.call(rbind, theta.pri.mu)
  theta.pri.prod <- lapply(1:ntimes, function(x) theta.pri.pre[[x]]%*%theta.pri.mu[[x]])
  
  # -------- tuning parameters ------- #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = length(y), nrow = niter + 1)
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  colnames(accept.mat) <- c("xi")
  
  ## list of matrices to store MCMC samples of basis coefs used for
  ## modeling quantile functions 
  theta.mat.list <- replicate(1 + niter, 
                              matrix(NA, ncol = L + 1, nrow = ntimes), 
                              simplify = F)
  theta.mat.list[[1]] <- theta.pri.mu.mat
  
  if(! rand.int){
    parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
    colnames(parm.mat) <- c("xi")
    parm.mat[1, ] <- c(5)
  }else if(rand.int & str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2", "xi")
    parm.mat[1, ] <- c(0.1, 5)
    ## data frame to store random effects
    theta.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
    theta.mat[,1] <- 0
    ## the InvGamma prior for the variance of temporal random effects 
    a = 0.1; b = 0.1
  }else if(rand.int & str == "AR-1"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2", "rho", "xi")
    parm.mat[1, ] <- c(0.1, 0.9, 5)
    ## data frame to store random effects
    theta.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
    theta.mat[,1] <- 0
    ## adjacency matrix ##
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))
    ## the InvGamma prior for variance of temporal random effects 
    a = 0.1; b = 0.1
    ## the discrete prior for rho controlling temporal dependency 
    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)),
                    simplify = TRUE)
  }
  
  ## dimension of BK.mat = int B(tau)K(tau)^T dtau  = (L+1)xK
  if(is.null(Z.design)){
    X.design.curr = cbind(1, theta.mat.list[[1]]%*%BK.mat)
    dat.sim = data.frame(y = y, X.design.curr)
    colnames(dat.sim) <- c("y", "Int", paste0("x", 1:(num.beta-1)))
    pg.parm1 <- dat.sim$y + parm.mat[1, "xi"]
    pg.parm2 <- as.matrix(dat.sim[,c("Int", paste0("x", 1:(num.beta-1)))])%*%
      matrix(beta.mat[1, ], ncol = 1)
    ntotal = ntimes
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    X.design.curr = theta.mat.list[[1]]%*%BK.mat
    dat.sim = data.frame(y = y, 1, Z.design, X.design.curr)
    colnames(dat.sim) <- c("y", "Int", paste0("z", 1:num.beta.z),
                           paste0("x", 1:num.beta.x))
    pg.parm1 <- dat.sim$y + parm.mat[1, "xi"]
    index = ! (colnames(dat.sim) %in% "y")
    pg.parm2 <- as.matrix(dat.sim[,index])%*%matrix(beta.mat[1, ], ncol = 1)
    ntotal = ntimes
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    X.design.curr <- theta.mat.list[[i]]%*%BK.mat
    if(is.null(Z.design)){
      X.design.all <- cbind(1, X.design.curr)
    }else{
      X.design.all <- cbind(1, Z.design, X.design.curr)
    }
    
    ## compute latent normal variables ##
    z.vec <- (y - parm.mat[i, "xi"])/(2*omega.vec)
    if(! rand.int){
      M <- solve(crossprod(X.design.all*sqrt(omega.vec)) + C.beta.pre)
      m <- M%*%(C.beta.pre%*%c.beta + 
                  t(sqrt(omega.vec)*X.design.all)%*%
                  (sqrt(omega.vec)*matrix(z.vec, ncol = 1)))
      beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    }else{
      res.vec = z.vec - theta.mat[i,]
      M <- solve(crossprod(X.design.all*sqrt(omega.vec)) + C.beta.pre)
      m <- M%*%(C.beta.pre%*%c.beta + 
                  t(sqrt(omega.vec)*X.design.all)%*%
                  (sqrt(omega.vec)*matrix(res.vec, ncol = 1)))
      beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    }
  
    # 3. update basis coefficients for exposure quantile functions
    if(is.null(Z.design)){
      fix.eff = beta.mat[i+1, "beta0"]
      z.theta.vec <- z.vec - fix.eff
      BK.beta.i <- t(BK.mat%*%matrix(beta.mat[i+1, paste0("beta_x", 1:num.beta.x)], ncol = 1))
      M.list <- get_Mlist(BK_beta = BK.beta.i,
                          omega_vec = omega.vec,
                          theta_pre = theta.pri.pre)
      m.list <- get_mulist(theta_prod = theta.pri.prod,
                           M_list = M.list,
                           omega_vec = omega.vec,
                           z = z.theta.vec,
                           BK_beta = BK.beta.i)
      theta.mat.quan.up <- sample_basiscoef(mu_list = m.list,
                                            sigma_list = M.list,
                                            ntheta = L+1)
      theta.mat.list[[i+1]] <- theta.mat.quan.up
      X.design.up <- theta.mat.quan.up%*%BK.mat
    }else{
      theta.vec.i = theta.mat[i,]
      fix.eff = beta.mat[i+1, "beta0"] + 
        as.numeric(as.matrix(Z.design)%*%
                     matrix(beta.mat[i+1, paste0("beta_z", 1:num.beta.z)], 
                            ncol = 1))
      z.theta.vec <- z.vec - fix.eff - theta.vec.i
      BK.beta.i <- t(BK.mat%*%matrix(beta.mat[i+1, paste0("beta_x", 1:num.beta.x)], ncol = 1))
      M.list <- get_Mlist(BK_beta = BK.beta.i,
                          omega_vec = omega.vec,
                          theta_pre = theta.pri.pre)
      m.list <- get_mulist(theta_prod = theta.pri.prod,
                           M_list = M.list,
                           omega_vec = omega.vec,
                           z = z.theta.vec,
                           BK_beta = BK.beta.i)
      theta.mat.quan.up <- sample_basiscoef(mu_list = m.list,
                                            sigma_list = M.list,
                                            ntheta = L+1)
      theta.mat.list[[i+1]] <- theta.mat.quan.up
      X.design.up <- theta.mat.quan.up%*%BK.mat
    }
    
    if(is.null(Z.design)){
      X.design.all <- cbind(1, X.design.up)
    }else{
      X.design.all <- cbind(1, Z.design, X.design.up)
    }
    
    # update random effects (conjugate)
    ## posterior mean and variance of random effects ##
    if(rand.int){
      if(str == "AR-1"){
        pre.theta.0 <- (1/parm.mat[i, "tau2"])*(Dw.t - parm.mat[i, "rho"]*W.t)
        pre.theta.post <- Omega + pre.theta.0
        Sigma.theta.post <- solve(pre.theta.post)
        fix.eff <- X.design.all%*%matrix(beta.mat[i+1, ], ncol = 1)
        mu.theta.post <- Sigma.theta.post%*%(Omega%*%matrix(z.vec - fix.eff, ncol = 1))
        theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.theta.post,
                                                Sigma = Sigma.theta.post) )
        ## update tau2 
        parm.mat[i+1, "tau2"] <- 1/rgamma(1, ntimes/2 + a,
                                          b + as.numeric(sum(theta.mat[i+1, ]^2))/2)
        
      }else if(str == "exch"){
        pre.theta.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2"]), ntimes))
        pre.theta.post <- Omega + pre.theta.0
        Sigma.theta.post <- solve(pre.theta.post)
        fix.eff <- X.design.all%*%matrix(beta.mat[i+1, ], ncol = 1)
        mu.theta.post <- Sigma.theta.post%*%(Omega%*%matrix(z.vec - fix.eff, ncol = 1))
        theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.theta.post,
                                                Sigma = Sigma.theta.post) )
        ## update tau2 and rho
        theta.unique <- matrix(theta.mat[i+1, ], ncol = 1)
        parm.mat[i+1, "tau2"] <- 1/rgamma(1, ntimes/2 + a1,
                                          b1 + as.numeric(t(theta.unique)%*%
                                                            (Dw.t - parm.mat[i, "rho"]*W.t)%*%theta.unique)/2)
        inter = as.numeric( t(matrix(theta.unique, ncol = 1))%*%W.t%*%
                              matrix(theta.unique, ncol = 1) )
        ll.rho = ll.rho.1 + rho.prior.val/(2*parm.mat[i+1, "tau2"]^2)*inter
        parm.mat[i+1, "rho"] <- sample(x = rho.prior.val, size = 1,
                                       prob = exp(ll.rho - max(ll.rho)))
        
      }
    }
  
    # 4. update xi (over dispersion) 
    eta.vec.i1 = X.design.all%*%matrix(beta.mat[i+1, ], ncol = 1)
    q.nb = 1/(1+exp(eta.vec.i1)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(dat.sim$y, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(dat.sim$y, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    omega.mat[i+1, ] <- rpg(ntotal, dat.sim$y + parm.mat[i+1, "xi"], eta.vec.i1)
    
    ## tuning parameters ##
    if(i <= 5000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
    
  }# END MCMC iterations 

  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  if(! rand.int){
    beta.post <- beta.mat[-(1:burn_in), ]
    X.quan.post <- lapply(theta.mat.list[-c(1:burn_in)], function(x) x%*%BK.mat)
    if(is.null(Z.design)){
      q.vec <- compq_cpp_XZonly(beta_x_post = beta.post[,paste0("beta_x", 1:num.beta.x)],
                                beta_z_post = as.matrix(beta.post[,"beta0"], ncol = 1),
                                Xdesign_quan_post = X.quan.post,
                                Zdesign = as.matrix(dat.sim[,"Int"], ncol = 1),
                                ntotal = ntotal, npost = nrow(beta.post))
    }else{
      q.vec <- compq_cpp_XZonly(beta_x_post = beta.post[,paste0("beta_x", 1:num.beta.x)],
                                beta_z_post = beta.post[,c("beta0",paste0("beta_z", 1:num.beta.x))],
                                Xdesign_quan_post = X.quan.post,
                                Zdesign = cbind(1, Z.design),
                                ntotal = ntotal, npost = nrow(beta.post))
    }
    log.pd = t( apply(cbind(y, q.vec), 1, compll,
                      xi.post = parm.mat[-c(1:burn_in),"xi"]) )
    lppd = sum(log(apply(exp(log.pd), 1, mean)))
    pWAIC = sum(apply(log.pd, 1, var))
    
    results <- list(beta = beta.mat,
                    parm = parm.mat,
                    basis.coef.quan = theta.mat.list,
                    WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
    
  }else{
    beta.post <- beta.mat[-(1:burn_in), ]
    rand.post <- theta.mat[-c(1:burn_in), ]
    X.quan.post <- lapply(theta.quan.list[-c(1:burn_in)], function(x) x%*%BK.mat)
    if(is.null(Z.design)){
      q.vec <- compq_cpp(beta_x_post = beta.post[,paste0("beta_x", 1:num.beta.x)],
                         beta_z_post = as.matrix(beta.post[,"beta0"], ncol = 1),
                         rand_post = rand.post,
                         Xdesign_quan_post = X.quan.post,
                         Zdesign = as.matrix(dat.sim[,"Int"], ncol = 1),
                         ntotal = ntotal, npost = nrow(beta.post))
    }else{
      q.vec <- compq_cpp(beta_x_post = beta.post[,paste0("beta_x", 1:num.beta.x)],
                         beta_z_post = beta.post[,c("beta0",paste0("beta_z", 1:num.beta.x))],
                         rand_post = rand.post,
                         Xdesign_quan_post = X.quan.post,
                         Zdesign = cbind(1, Z.design),
                         ntotal = ntotal, npost = nrow(beta.post))
    }
    log.pd = t( apply(cbind(y, q.vec), 1, compll,
                      xi.post = parm.mat[-c(1:burn_in),"xi"]) )
    lppd = sum(log(apply(exp(log.pd), 1, mean)))
    pWAIC = sum(apply(log.pd, 1, var))
    
    results <- list(beta = beta.mat,
                    parm = parm.mat,
                    basis.coef.quan = theta.mat.list,
                    rand.intercept = theta.mat,
                    WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
    
  }

  return(results)
  
  ## OUTPUTS
  ## A list containing 
  ## if the argument rand.ind is FALSE: 
  ## (1) a matrix containing MCMC samples for 
  ## (basis coefs for modeling beta(tau);
  ##  coefs of confounders;
  ##  over-dispersion parameter;
  ##  basis coefs for modeling exposure quantile functions) 
  ## (2) WAIC 
  ## if the argument rand.ind is TRUE: 
  ## (1) a matrix containing MCMC samples for
  ## (basis coefs for modeling beta(tau);
  ##  coefs of confounders;
  ##  basis coefs for modeling exposure quantile functions;
  ##  random intercepts;
  ##  over-dispersion parameter;
  ##  parameters related to the mean-zero Gaussian process specified for modeling residuals) 
  ## (2) WAIC
  return(re)
}





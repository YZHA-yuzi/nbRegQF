############################################################################
### Functions to estimate quantile functions using individual-level exposure data
### yuzi.zhang@emory.edu
############################################################################

### A function to find initial values of basis coefs using quadratic programming 
find.initial <- function(y.obs, Bl.mat, L, tau.vec){
  Dmat <- t(Bl.mat) %*% Bl.mat
  dvec <- t(Bl.mat) %*% matrix(y.obs, ncol = 1)
  Amat <- diag(rep(1, L+1))
  bvec <- matrix(0, ncol = L + 1)
  theta.initial <- solve.QP(Dmat = Dmat, dvec = dvec,
                            Amat = Amat, bvec = bvec)$solution
  return(theta.initial)
}

ff <- function(x, epsilon){
  x[x <= epsilon] <- epsilon
  return(x)
}

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




#' Quantile function specified via the basis expansion 
#'
#' This function evaluates the quantile function specified via the basis expansion using 
#' piecewise Gaussian functions or piecewise Gamma functions.   
#'
#' @param tau a vector of probabilities (i.e., percentile levels). 
#' @param L a number specifying the number of of basis functions used in the 
#' basis expansion for defining the quantile function. 
#' @param theta.vec a vector of lenght \code{L} containing basis coefficients determing 
#' the shape of the distribution.
#' @param alpha a numeric value specifying the median of the distribution. 
#' @param basis.fun a character string specifying which basis functions are used 
#' for modeling quantile processes. This function accepts \code{Gau} to specify 
#' piecewise Gaussian functions and \code{Gamma} to specify piecewise Gamma functions.
#' @param shape a number specifying the shape parameter of the 
#' piecewise Gamma functions if \code{basis.fun = Gamma}, the default value is 5.
#' 
#' @examples
#' # Evalute the quantile function of a random variable following a normal distribution with 
#' # mean 0 and variance 1, but specify its quantile function via the basis expansion using 
#' # 4 piecewise Gaussian functions. 
#' quan.fun(tau = seq(0, 1, 0.2), L = 4, theta.vec = rep(1, 4), alpha = 0, basis.fun = "Gau")
#' # This function should return the same value as applying the function
#' qnorm(p = seq(0, 1, 0.2), mean = 0, sd = 1)
#' 
#' @export
## define quantile function based on those basis functions ##
quan.fun <- function(tau, L, theta.vec, alpha, basis.fun, shape = 5){
  if(basis.fun == "Gau"){
    B.mat <- matrix(NA, nrow = L, ncol = length(tau))
    for(i in 1:L){
      B.mat[i, ] <- sapply(tau, Bl, l = i, L = L, 
                           basis.fun = basis.fun,
                           shape = shape,
                           simplify = T)
    }
    theta.vec <- matrix(theta.vec, nrow = 1)
    Q.sim <- (alpha + theta.vec%*%B.mat)
    return(Q.sim)
  }
  if(basis.fun == "Gamma"){
    B.mat <- matrix(NA, nrow = L, ncol = length(tau))
    for(i in 1:L){
      B.mat[i, ] <- sapply(tau, Bl, l = i, L = L, 
                           basis.fun = basis.fun,
                           shape = shape,
                           simplify = T)
    }
    theta.vec <- matrix(theta.vec, nrow = 1)
    Q.sim <- (alpha + theta.vec%*%B.mat)
    return(Q.sim)
  }
}


#' Density function specified via the basis expansion 
#'
#' This function evaluates the density function specified via the basis expansion using 
#' piecewise Gaussian functions or piecewise Gamma functions.   
#'
#' @param x a vector of quantiles. 
#' @param L a number specifying the number of of basis functions used in the 
#' basis expansion for defining the quantile function. 
#' @param theta.vec a vector of lenght \code{L} containing basis coefficients determing 
#' the shape of the distribution.
#' @param alpha a numeric value specifying the median of the distribution. 
#' @param basis.fun a character string specifying which basis functions are used 
#' for modeling quantile processes. This function accepts \code{Gau} to specify 
#' piecewise Gaussian functions and \code{Gamma} to specify piecewise Gamma functions.
#' @param shape a number specifying the shape parameter of the 
#' piecewise Gamma functions if \code{basis.fun = Gamma}, the default value is 5.
#' 
#' @examples
#' # Evalute the density function of a random variable following a normal distribution with 
#' # mean 0 and variance 1, but specify its quantile function via the basis expansion using 
#' # 4 piecewise Gaussian functions. 
#' den.fun(x = 0.1, L = 4, theta.vec = rep(1, 4), alpha = 0, basis.fun = "Gau")
#' # This function should return the same value as applying the function
#' dnorm(x = 0.1, mean = 0, sd = 1)
#' 
#' @export
## A function to evaluate density function ##
den.fun <- function(x, L, theta.vec, alpha, basis.fun, shape = 5){
  k.vec = seq(0, 1, length.out = L+1)
  if(basis.fun == "Gau"){
    q.vec <- quan.fun(tau = k.vec, L = L, 
                      theta.vec = theta.vec,
                      alpha = alpha,
                      basis.fun = basis.fun,
                      shape = shape)
    interval.mat <- cbind(q.vec[1:L], q.vec[2:(L+1)])
    a.vec <- get_basis_a(k.vec, q.vec, theta.vec)
    den.vec <- get_den(x, interval.mat, avec = a.vec,
                       thetavec = theta.vec)
    return(den.vec)
  }
  if(basis.fun == "Gamma"){
    q.vec <- quan.fun(tau = k.vec, L = L, 
                      theta.vec = theta.vec,
                      alpha = alpha, 
                      basis.fun = basis.fun,
                      shape = shape)
    interval.mat <- cbind(q.vec[1:L], q.vec[2:(L+1)])
    a.vec <- get_basis_a_gamma(k.vec, q.vec, theta.vec, shape)
    den.vec <- get_den_gamma(x, interval.mat, avec = a.vec,
                             thetavec = theta.vec, shape = shape)
    return(den.vec)
  }
}

# A function to estimate quantile function at a time point 
fit.exposure.inde <- function(x.ind, 
                              basis.fun, basis.shape = 5, L = 4,
                              niter, ...){
  ## INPUTs:
  ## x.ind: a vector containing individual-level exposures 
  ## basis.fun: a character string specifying what basis functions are used
  ## for expanding quantile functions 
  ## "Gau" = Gaussian piecewise functions
  ## "Gamma" = Gamma piecewise functions
  ## if basis.fun = "Gamma", basis.shape is the shape parameter specified for
  ## the Gamma piecewise functions; default value is 5
  ## L: the number of basis functions used in the basis expansion  
  ## niter: the number of MCMC iterations 
  
  ## find good initial values for basis coefs ##
  tau.vec.int = seq(0.01, 0.99, length.out = 50)
  Bl.mat <- cbind(1, do.call(cbind, lapply(1:L, function(x) 
    Bl(tau.vec = tau.vec.int, l = x, L = L, 
       basis.fun = basis.fun, shape = basis.shape))))
  Bl.mat <- as.matrix(Bl.mat)
  
  y.obs <- as.numeric(quantile(x.ind, tau.vec.int, na.rm = T))
  theta.initial <- find.initial(y.obs = y.obs,
                                Bl.mat = Bl.mat, L = L, tau.vec = tau.vec.int)
  
  # -------- PRIORS ------- #
  ## specify priors for basis coefs ##
  c.prior = 0; C.prior = 100
  
  # -------- tuning parameters -------- #
  alph.tun = 0.01
  theta.star.tun = rep(0.01, L)
  
  # -------- initialize parameters --------- #
  parm.mat <- matrix(NA, ncol = L + 1, nrow = niter + 1)
  colnames(parm.mat) <- c("alpha", paste0("theta", 1:L))
  theta.star.mat <- matrix(NA, ncol = L, nrow = niter + 1)
  parm.mat[1, ] <- theta.initial
  theta.star.mat[1, ] <- theta.initial[-1]
  accept.mat <- matrix(NA, ncol = 2, nrow = niter+1)
  epsilon = 0.01
  
  # ----- BEGIN MCMC -------- #
  for(i in 1:niter){
    
    # 1. update alpha_t,0 (ntimes by 1) M-H
    ## current value of alpha_t,0 and theta_t,l
    theta.i <- ff(theta.star.mat[i, ], epsilon = epsilon)
    alph.i <- parm.mat[i, "alpha"]
    
    ## current log-likelihood
    ll.curr.0.vec = log(den.fun(x = na.omit(x.ind),
                                L = L,
                                theta.vec = theta.i,
                                alpha = alph.i, 
                                basis.fun = basis.fun, 
                                shape = basis.shape))
    ll.curr.0.vec[!is.finite(ll.curr.0.vec)] <- log(1e-20)
    ll.curr.0 <- sum(ll.curr.0.vec)

    ## propose a candidate alpha_t0*
    alph.star <- rnorm(1, alph.i, sd = alph.tun)
    ## likelihood for the proposed alpha
    ll.star.0.vec = log(den.fun(x = na.omit(x.ind),
                                L = L,
                                theta.vec = theta.i,
                                alpha = alph.star, 
                                basis.fun = basis.fun, 
                                shape = basis.shape))
    ll.star.0.vec[!is.finite(ll.star.0.vec)] <- log(1e-20)
    ll.star.0 <- sum(ll.star.0.vec)

    ll.curr <- ll.curr.0 + dnorm(x = alph.star,
                                 mean = c.prior,
                                 sd = C.prior, log = T)
    ll.star <- ll.star.0 + dnorm(x = alph.i,
                                 mean = c.prior,
                                 sd = C.prior, log = T)
    
    ratio = min(1, exp(ll.star - ll.curr))
    if(ratio >= runif(1)){
      parm.mat[i+1, "alpha"] <- alph.star
      accept.mat[i+1, 1] <- 1
      ll.curr.0 <- ll.star.0
    }else{
      parm.mat[i+1, "alpha"] <- alph.i
      accept.mat[i+1, 1] <- 0
    }
    
    # 2. update theta_tl* t = itime, l = 1, ..., L
    ## update theta_t as a Lx1 vector
    theta.star.i <- theta.star.mat[i, ]
    theta.star.star <- rnorm(L, theta.star.i, theta.star.tun)
    
    theta.i <- ff(theta.star.i, epsilon = epsilon)
    theta.star <- ff(theta.star.star, epsilon = epsilon)
    
    ll.star.0.vec = log(den.fun(x = na.omit(x.ind),
                                L = L,
                                theta.vec = theta.star,
                                alpha = parm.mat[i+1, "alpha"],
                                basis.fun = basis.fun, 
                                shape = basis.shape))
    ll.star.0.vec[!is.finite(ll.star.0.vec)] <- log(1e-20)
    ll.star.0 = sum(ll.star.0.vec)

    ll.curr <- ll.curr.0 + sum(dnorm(x = theta.star.star,
                                     mean = c.prior,
                                     sd = C.prior, log = T))
    ll.star <- ll.star.0 + sum(dnorm(x = theta.star.i,
                                     mean = c.prior,
                                     sd = C.prior, log = T))
    
    ratio = min(1, exp(ll.star - ll.curr))
    if(ratio >= runif(1)){
      theta.star.mat[i+1, ] <- theta.star.star
      parm.mat[i+1, 2:(L+1)] <- theta.star
      accept.mat[i+1, 2] <- 1
      ll.curr.0 <- ll.star.0
    }else{
      theta.star.mat[i+1, ] <- theta.star.i
      parm.mat[i+1, 2:(L+1)] <- theta.i
      accept.mat[i+1, 2] <- 0
    }
    
    ## tuning parameters ##
    if(i <= 5000 & i%%1000 == 0){
      accept.rate.all <- colMeans(accept.mat[2:i, ])
      
      accept.rate = accept.rate.all[1]
      if(accept.rate > 0.45){alph.tun = 1.2*alph.tun}
      if(accept.rate < 0.25){alph.tun = 0.8*alph.tun}
      
      accept.rate = accept.rate.all[2]
      if(accept.rate > 0.55){
        theta.star.tun = 1.2*theta.star.tun
      }
      if(accept.rate < 0.25){
        theta.star.tun = 0.8*theta.star.tun
      }
    } # end tuning
    
  } # end MCMC procedure
  # print(paste0((proc.time()[3] - start)/60, "min"))
  
  re  <- list(parm = parm.mat,
              theta.star = theta.star.mat,
              accept = accept.mat)
  return(re)
  ## OUTPUT:
  ## a list containing MCMC samples of basis coefs (latent variables), and
  ## a matrix containing whether the MCMC sample was accepted in MH algorithms  
}

# A function to estimate quantile function assuming quantile processes are
# temporal correlated by assuming basis coefs follow a first-order Gaussian 
# process  
fit.exposure.AR1 <- function(x.ind, 
                             basis.fun, basis.shape = 5, L = 4, 
                             niter, ...){
  ## INPUTs:
  ## x.ind: a matrix containing individual-level exposures
  ## each row contains individual-level exposures collected from a time point  
  ## basis.fun: a character string specifying what basis functions are used
  ## for expanding quantile functions 
  ## "Gau" = Gaussian piecewise functions
  ## "Gamma" = Gamma piecewise functions
  ## if basis.fun = "Gamma", basis.shape is the shape parameter specified for
  ## the Gamma piecewise functions; default value is 5
  ## L: the number of basis functions used in the basis expansion  
  ## niter: the number of MCMC iterations 
  
  ntimes = nrow(x.ind)
  
  #### find good initial values ###
  tau.vec.int <- seq(0.01, 0.99, length.out = 50)
  Bl.mat <- cbind(1, do.call(cbind, lapply(1:L, function(x) 
    Bl(tau.vec = tau.vec.int, l = x, L = L, 
       basis.fun = basis.fun, shape = basis.shape))))
  Bl.mat <- as.matrix(Bl.mat)
  
  theta.init.mat <- matrix(NA, ncol = L + 1, nrow = ntimes)
  for(itime in 1:ntimes){
    y.obs <- as.numeric(quantile(x.ind[itime, ], tau.vec.int, na.rm = T))
    theta.initial <- find.initial(y.obs = y.obs,
                                  Bl.mat = Bl.mat, L = L, tau.vec = tau.vec.int)
    theta.init.mat[itime, ] <- theta.initial
  }
  theta.bar.init <- find.initial(y.obs = 
                                   as.numeric(quantile(as.numeric(x.ind), 
                                                       tau.vec.int, na.rm = T)),
                                 Bl.mat = Bl.mat, L = L, tau.vec = tau.vec.int)
  
  ## adjacency matrix ##
  adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
  W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
  Dw.t <- Diagonal(x = apply(W.t, 1, sum))
  
  # ---------- specify priors --------- #
  ## alpha0, theta_l_bar
  c.prior = 0; C.prior = 100
  ## tau12, tau22 
  a = 0.1; b = 0.1
  ## rho1, rho2 
  # compute quantities beforehand to facilitate the updating of rho1, and rho2 
  lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
  rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
  ll.rho1.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)), 
                     simplify = TRUE)
  
  rho2.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
  ll.rho2.1 = sapply(rho2.prior.val, function(x) 0.5*L*sum(log(1-x*lambda.t)), 
                     simplify = TRUE)
  
  # --------- tuning parameters -------- #
  alph.tun = rep(0.1, ntimes)
  theta.star.tun = matrix(rep(0.1, ntimes*L),
                          nrow = L, ncol = ntimes)
  alph <- matrix(NA, ncol = ntimes, nrow = niter+1); 
  theta.star.list <- replicate(L, 
                               matrix(NA, ncol = ntimes, nrow = niter+1), 
                               simplify = F)
  theta.bar.mat <- matrix(NA, ncol = L, nrow = niter + 1)
  colnames(theta.bar.mat) <- paste0("theta", 1:L)
  parm.mat <- matrix(NA, ncol = 5, nrow = niter + 1)
  colnames(parm.mat) <- c("alpha", "tau12", "tau22", 
                          "rho1", "rho2")
  
  ### assign initial values ###
  for(ll in 1:L){
    theta.star.list[[ll]][1, ] <- theta.init.mat[,ll+1] 
  }
  alph[1, ] <- theta.init.mat[,1] 
  theta.bar.mat[1, ] <- colMeans(theta.init.mat)[-1]
  parm.mat[1,  ] <- c(colMeans(theta.init.mat)[1], rep(0.5, ncol(parm.mat)-1))
  accept.mat <- matrix(NA, ncol = ntimes + ntimes, nrow = niter+1)
  epsilon = 0.01
  
  # start = proc.time()[3]
  for(i in 1:niter){
    
    # 1. update alpha_t,0 (ntimes by 1) M-H
    theta.start.mat.i <- matrix(NA, nrow = L, ncol = ntimes)
    for(ll in 1:L){
      theta.start.mat.i[ll, ] <- theta.star.list[[ll]][i, ]
    }
    theta.mat.i <- apply(theta.start.mat.i, 2, ff, epsilon = epsilon)
    Omega1.i <- as.matrix(1/parm.mat[i, "tau12"]*(Dw.t - parm.mat[i,"rho1"]*W.t))
    
    ## current value of alpha = (alpha_10, ..., alpha_T0) (T x 1)
    alph.i <- alph[i, ]
    alph.val.i.vec <- rep(parm.mat[i, "alpha"], ntimes)
    
    ll.curr.vec <- comp_ll_alpha(alph_c = alph.i,
                                 x_sim_mat = x.ind,
                                 theta_mat_i = theta.mat.i,
                                 f = den.fun, L = L, 
                                 shape = basis.shape,
                                 ntimes = ntimes, 
                                 basis_fun = basis.fun)
    ll.curr.vec[!is.finite(ll.curr.vec)] <- log(1e-200)
    
    ### pre-compute the conditional mean and variance for alpha_t|alpha_-t ###
    alph.cond.mu <- alph.cond.sd <- rep(NA, ntimes)
    alph.bar.i = parm.mat[i, "alpha"]
    for(a.ii.cond in 1:ntimes){
      index2 = (1:ntimes)[-a.ii.cond]
      ### using precision matrix ###
      omega12.ii = Omega1.i[a.ii.cond, index2]
      alph.cond.mu[a.ii.cond] = alph.bar.i - (1/Omega1.i[a.ii.cond, a.ii.cond])*
        sum(omega12.ii*(alph.i[-a.ii.cond] - alph.bar.i))
      alph.cond.sd[a.ii.cond] = sqrt(1/Omega1.i[a.ii.cond, a.ii.cond])
    }
    
    for(a.ii in 1:ntimes){
      ## propose a candidate alpha_t0*
      alph.star.ii <- rnorm(1, alph.i[a.ii], sd = alph.tun[a.ii])
      alph.star <- alph.i 
      alph.star[a.ii] <- alph.star.ii
      ## evaluate likelihood 
      ## likelihood for the proposed alpha
      ll.star.t <- sum(log(den.fun(x = na.omit(x.ind[a.ii, ]), L = L,
                                   theta.vec = theta.mat.i[ ,a.ii],
                                   alpha = alph.star[a.ii],
                                   shape = basis.shape, 
                                   basis.fun = basis.fun)))
      if(!is.finite(ll.star.t)){ll.star.t <- log(1e-200)}
      
      ll.star <- ll.star.t + dnorm(alph.star.ii, 
                                   alph.cond.mu[a.ii], 
                                   alph.cond.sd[a.ii], log = T)
      
      ll.curr <- ll.curr.vec[a.ii] + dnorm(alph.i[a.ii], 
                                           alph.cond.mu[a.ii], 
                                           alph.cond.sd[a.ii], log = T)
      
      ratio = min(1, exp(ll.star - ll.curr))
      if(ratio >= runif(1)){
        alph.i[a.ii] <- alph.star.ii
        ll.curr.vec[a.ii] <- ll.star.t
        accept.mat[i+1, a.ii] <- 1
      }else{
        accept.mat[i+1, a.ii] <- 0
      }
      # cat(a.ii,",")
    } # end loop over alpha_t0, t = 1, ..., T
    alph[i+1, ] <- alph.i
    
    # 2. update theta_tl* t = 1, ..., ntimes, l = 1, ..., L
    ## update theta_t as a Lx1 vector 
    Omega2.i.sparse <- 1/parm.mat[i, "tau22"]*(Dw.t - parm.mat[i,"rho2"]*W.t)
    Omega2.i <- as.matrix(Omega2.i.sparse)
    theta.start.mat.i <- matrix(NA, nrow = L, ncol = ntimes)
    for(ll in 1:L){
      theta.start.mat.i[ll, ] <- theta.star.list[[ll]][i, ]
    }
    
    ### pre-compute the conditional variance for theta_tl|theta_-tl ###
    theta.cond.sd <- sqrt(1/diag(Omega2.i))
    theta.bar.i = theta.bar.mat[i, ]
    for(the.ii in 1:ntimes){
      ## propose a new theta_lt*|current values 
      theta.star.star.ii <- rnorm(L, theta.start.mat.i[,the.ii],
                                  theta.star.tun[,the.ii])
      theta.star.ii <- ff(theta.star.star.ii, epsilon = epsilon)
      ll.star.t <- sum(log(den.fun(x = na.omit(x.ind[the.ii, ]), L = L,
                                   theta.vec = theta.star.ii,
                                   alpha = alph.i[the.ii],
                                   shape = basis.shape,
                                   basis.fun = basis.fun)))
      if(!is.finite(ll.star.t)){ll.star.t <- log(1e-200)}
      
      ### compute conditional mean of theta_tl|theta_-tl using precision matrix ###
      ### compute it within loop, 
      ### since the conditional mean depends on current values of theta_-tl
      index2 = (1:ntimes)[-the.ii]
      omega12.ii = Omega2.i[the.ii, index2]
      H.inv = 1/Omega2.i[the.ii, the.ii]
      theta.cond.mu.t <- rep(NA, L)
      for(lll in 1:L){
        theta.cond.mu.t[lll] = theta.bar.i[lll] - 
          H.inv*sum(omega12.ii*(theta.start.mat.i[lll, -the.ii] - theta.bar.i[lll]))
      }
      
      p.star = sum(dnorm(theta.star.star.ii, theta.cond.mu.t,
                         rep(theta.cond.sd[the.ii], L), log = T))
      p.curr = sum(dnorm(theta.start.mat.i[,the.ii], theta.cond.mu.t,
                         rep(theta.cond.sd[the.ii], L), log = T))
      ll.star <- ll.star.t + p.star
      ll.curr <- ll.curr.vec[the.ii] + p.curr
      ratio = min(1, exp(ll.star - ll.curr))
      if(ratio >= runif(1)){
        theta.start.mat.i[,the.ii] <- theta.star.star.ii
        accept.mat[i+1, ntimes + the.ii] <- 1
        ll.curr.vec[the.ii] <- ll.star.t
      }else{
        accept.mat[i+1, ntimes + the.ii] <- 0
      }
      # cat(the.ii, ",")
    } # end loop over time for theta_lt*
    for(ll in 1:L){
      theta.star.list[[ll]][i+1, ] <- theta.start.mat.i[ll, ] 
    }
    
    # 3. update theta_bar_l, l = 1,..., L conjugate 
    X.design.theta = matrix(rep(1, ntimes), ncol = 1)
    pre.comp = t(X.design.theta) %*% Omega2.i
    for(ll in 1:L){
      M = 1/(pre.comp %*% X.design.theta + 1/C.prior)
      m = M*(pre.comp%*%matrix(theta.start.mat.i[ll, ], ncol = 1))
      theta.bar.mat[i+1, ll] <- rnorm(1, m, sqrt(M))
    }
    
    # 4. update alpha_bar conjugate 
    X.design.alph = matrix(rep(1, ntimes), ncol = 1)
    pre.comp = t(X.design.alph) %*% Omega1.i
    M = 1/(pre.comp %*% X.design.alph + 1/C.prior)
    m = M*(pre.comp%*%matrix(alph[i+1, ], ncol = 1))
    parm.mat[i+1, "alpha"] <- rnorm(1, m, sqrt(M))
    
    # 5. update tau12, conjugate, inverse gamma 
    # variance parameter of the intercept 
    res.alph.i <- alph[i+1, ] - parm.mat[i+1, "alpha"]
    parm.mat[i+1, "tau12"] <- 1/rgamma(1, ntimes/2 + a,
                                       b + as.numeric(t(matrix(res.alph.i, ncol = 1))%*%
                                                        (Dw.t - parm.mat[i, "rho1"]*W.t)%*%
                                                        matrix(res.alph.i, ncol = 1))/2)
    
    # 6. update tau22, conjugate, inverse gamma 
    # variance parameter of coefs of basis functions
    innter.mat = (Dw.t - parm.mat[i, "rho2"]*W.t)
    sum.ig = 0
    for(ll in 1:L){
      res.theta.ll = theta.star.list[[ll]][i+1, ] - theta.bar.mat[i+1, ll]
      re = as.numeric(t(matrix(res.theta.ll, ncol = 1))%*%
                        innter.mat%*%
                        matrix(res.theta.ll, ncol = 1))
      sum.ig = sum.ig + re
    }
    parm.mat[i+1, "tau22"] <- 1/rgamma(1, L*ntimes/2 + a,
                                       b + sum.ig/2)
    
    # 7. update rho1 (M-H)
    inter = matrix(alph[i+1, ] - parm.mat[i+1, "alpha"], ncol = 1)
    inter1 = as.numeric( t(inter)%*%W.t%*%inter )
    ll.rho1 = ll.rho1.1 + rho1.prior.val/(2*parm.mat[i+1, "tau12"])*inter1
    parm.mat[i+1, "rho1"] <- sample(x = rho1.prior.val, size = 1, 
                                    prob = exp(ll.rho1 - max(ll.rho1)))
    
    # 8. update rho2 (M-H)
    inter2 = 0
    for(ll in 1:L){
      inter = matrix(theta.star.list[[ll]][i+1, ] - 
                       theta.bar.mat[i+1, ll], ncol = 1)
      inter1 = as.numeric( t(inter)%*%W.t%*%inter )
      inter2 = inter2 + inter1
    }
    ll.rho2 = ll.rho2.1 + rho2.prior.val/(2*parm.mat[i+1, "tau22"])*inter2
    parm.mat[i+1, "rho2"] <- sample(x = rho2.prior.val, size = 1, 
                                    prob = exp(ll.rho2 - max(ll.rho2)))
    
    ## tuning parameters ##
    if(i <= 1000 & i%%100 == 0){
      accept.rate.all <- colMeans(accept.mat[2:i, ])
      for(acc.ii in 1:ntimes){
        accept.rate = accept.rate.all[acc.ii]
        if(accept.rate > 0.45){alph.tun[acc.ii] = 1.2*alph.tun[acc.ii]}
        if(accept.rate < 0.25){alph.tun[acc.ii] = 0.8*alph.tun[acc.ii]}
        acc.the.ii = acc.ii+ntimes
        accept.rate.the = accept.rate.all[acc.the.ii]
        if(accept.rate.the > 0.45){
          theta.star.tun[,acc.ii] = 1.2*theta.star.tun[,acc.ii]
        }
        if(accept.rate.the < 0.25){
          theta.star.tun[,acc.ii] = 0.8*theta.star.tun[,acc.ii]
        }
      }
    } # end tunning 
    
  } # end MCMC procedure
  
  theta.overall = cbind(parm.mat[,"alpha"], theta.bar.mat)
  colnames(theta.overall) <- paste0("theta", 0:L, "bar")
  re  <- list(theta0.time = alph,
              theta.overall = theta.overall,
              parm = parm.mat[,-1],
              theta.star.time = theta.star.list,
              accept = accept.mat)
  return(re)
  ## OUTPUT:
  ## a list containing MCMC samples of 
  ## overall and time-specific basis coefs (latent variables), 
  ## "theta0.time": time-specific basis coef, i.e., theta_0i
  ## "theta.overall": overall basis coefs, i.e., theta_l_bar, l=0,...,L
  ## "parm": parameters controlling temporal dependency of quantile processes
  ## "theta.star.time": time-specific latent basis coefs, i.e., theta*_li, l=1,...L
  ## "accept": a matrix containing whether the MCMC sample was accepted in MH algorithms
}


#' Estimate exposure quantile functions using semiparametric Bayesian approaches 
#'
#' This function estimates basis coefficients of basis functions used for expanding 
#' exposure quantile functions using 
#' semiparametric Bayesian approaches based on individual-level exposures.
#'
#' @param x.ind a matrix with each row contains individual-level exposures.
#' @param basis.fun a character string specifying which basis functions are used 
#' for modeling quantile processes. This function accepts \code{Gau} to specify 
#' piecewise Gaussian functions and \code{Gamma} to specify piecewise Gamma functions.
#' @param basis.shape a number specifying the shape parameter of the 
#' piecewise Gamma functions if \code{basis.fun = Gamma}, the default value is 5.
#' @param L a number specifying the number of of basis functions used in the 
#' basis expansion for modeling exposure quantile functions.   
#' @param niter a number specifying the number of MCMC iterations.
#' @param inde a logical value indicating whether exposure quantile processes are 
#' assumed to be independent across time points. 
#' \code{inde = TRUE}, this function estimates basis coefficients assuming 
#' exposure quantile processes are independent across time points; 
#' \code{inde = FALSE} this function estimates basis coefficients assuming 
#' exposure quantile processes are temporally correlated by introducing 1-st order
#' Gaussian Markov random process for the unconstrained latent basis coefficients. 
#' 
#' @return If \code{inde = TRUE}, this function returns a list contains: (1) a matrix containing samples of 
#' basis coefficients (\code{parm}), (2) a matrix containing samples of unconstrained latent basis coefficients 
#' (\code{theta.star}), and (3) a matrix containing whether the MCMC sample was accepted in MH algorithms (\code{accept}), 
#' and (3) a matrix containing samples of unconstrained latent variables introduced for basis coefficients; 
#' 
#' if \code{inde = FALSE}, this function returns a list contains: (1) a matrix containing 
#' samples of overall basis coefficients (\code{theta.overall}), (2) a matrix containing 
#' samples of parameters related to introduce temporal dependency of quantile processes (\code{parm}), 
#' (3) a list of length \code{L} with each element includes a matrix containing samples of time-specific 
#' \code{l}-th unconstrained latent basis coefficients (\code{theta.star.time}), 
#' (4) a matrix containing time-specific basis coefficients (i.e., intercept included in the basis expansion)
#' (\code{theta0.time}).
#' 
#' @references 
#' Reich, B. J. (2012). Spatiotemporal quantile regression for detecting 
#' distributional changes in environmental processes. 
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics), 61(4)}, 535-553.
#' 
#' @examples
#' # Simulate time-specific means and SDs of exposure quantile functions
#' mean.vec <- rnorm(1000, mean = 7.2, sd = 1)
#' sd.vec <- rnorm(1000, mean = 1, sd = 0.2)
#' 
#' # Simulate individual-level exposures 
#' x.sim.mat <- do.call(rbind, lapply(1:1000, 
#' function(x) rnorm(num.size.time, mean.vec[x], sd.vec[x])))
#' 
#' # Estimate exposure quantile functions at the time point 1 
#'   re.fit.exp.ind <- fit.exposure(x.ind = x.sim.mat[1, ],
#'   basis.fun = "Gau", L = 4, niter = 10000, inde = T)
#' 
#' @export
fit.exposure <- function(x.ind, 
                         basis.fun, basis.shape = 5, L = 4, 
                         niter,
                         inde = TRUE, ...){
  
  ## INPUTs:
  ## inde: a logical value indicating whether quantile processes are correlated
  ## x.ind: 
  ## if inde = FALSE a matrix containing individual-level exposures
  ## each row contains individual-level exposures collected from a time point
  ## if inde = TRUE a vector containing individual-level exposures
  ## basis.fun: a character string specifying what basis functions are used
  ## for expanding quantile functions 
  ## "Gau" = Gaussian piecewise functions
  ## "Gamma" = Gamma piecewise functions
  ## if basis.fun = "Gamma", basis.shape is the shape parameter specified for
  ## the Gamma piecewise functions; default value is 5
  ## L: the number of basis functions used in the basis expansion  
  ## niter: the number of MCMC iterations 

  if(inde){
    re <- fit.exposure.inde(x.ind = x.ind,
                            basis.fun = basis.fun, 
                            basis.shape = basis.shape, 
                            L = L, niter = niter)

  }else{
    re <- fit.exposure.AR1(x.ind = x.ind,
                           basis.fun = basis.fun, 
                           basis.shape = basis.shape, 
                           L = L, niter = niter)
  }
  return(re)
}


### prepare MVN prior of basis coefs for fitting health models ###
exposure.prior.inde <- function(re, burn_in = 5000, ...){
  ## INPUTs:
  ## re: a list of length T (# of time points) containing
  ## outputs from applying self-defined function fit.exposure to 
  ## estimate exposure quantile functions for each time point
  ## burn_in: # of burn-in samples (default is 5000)
  
  theta.pri.mu.list <- list()
  theta.pri.pre.list <- list()
 
  ntimes = length(re)
  for(itime in 1:ntimes){
    parm.post.t <- re[[itime]]$parm[-c(1:burn_in), ]
    inter.mu = as.numeric(colMeans(parm.post.t))
    inter.cov = as.matrix(cov(parm.post.t))
    colnames(inter.cov) <- rownames(inter.cov) <- NULL
    inter.cov[inter.cov == 0] <- 1e-8
    theta.pri.mu.list[[itime]] <- inter.mu
    
    if(det(inter.cov) == 0 | det(inter.cov) < 1e-20){
      inter.diag <- diag(inter.cov)
      inter.diag[inter.diag < 1e-10] <- 1e-8
      inter.cov <- diag(inter.diag)
      inter.inv <- solve(inter.cov)
      theta.pri.pre.list[[itime]] <- inter.inv
    }else{
      inter.inv <- solve(inter.cov)
      theta.pri.pre.list[[itime]] <- inter.inv
    }
  }
  return(list(theta.pri.mu = theta.pri.mu.list,
              theta.pri.pre = theta.pri.pre.list))
  ## OUTPUTs:
  ## a list containing two lists:
  ## theta.pri.mu: a list of length T (# of time points) contains of
  ## means of the MVN prior
  ## theta.pri.pre: a list of length T (# of time points) contains of
  ## the precision matrix of the MVN prior 
}


exposure.prior.AR1 <- function(re, burn_in, ...){
  ## INPUTs:
  ## re: a list containing the output from applying self-defined function 
  ## fit.exposure with argument equals FALSE
  ## to estimate exposure quantile functions 
  ## burn_in: # of burn-in samples (default is 5000)
  
  get_pointest_t <- function(re, burn_in){
    ntimes = ncol(re$theta0.time)
    alph.t <- re$theta0.time[-c(1:burn_in), ]
    theta.star.list <- re$theta.star.time
    theta.star.post <- lapply(theta.star.list, function(x) x[-c(1:burn_in), ])
    theta.t.post <- lapply(theta.star.post, ff, epsilon = 0.01)
    theta.t.hat <- sapply(theta.t.post, colMeans)
    inter <- data.frame(alphat = colMeans(alph.t), theta.t.hat)
    colnames(inter)[-1] = paste0("theta_", 1:(ncol(theta.t.hat)))
    ### var-cov matrix of alpha_tl, theta_tl, ... at a given zipcode 
    cov.i.0 <- lapply(1:ntimes, function(y) 
      cbind(alph.t[,y], 
            do.call(cbind, lapply(theta.t.post, function(x) x[,y]))))
    cov.i <- lapply(cov.i.0, cov)
    return(list(mu = inter, cov = cov.i))
  }
  
  inter.mu.cov <- get_pointest_t(re = re, burn_in = burn_in)
  inter.mu.mat <- apply(inter.mu.cov$mu, 1, function(x) as.numeric(x, nrow = 1))
  theta.pri.mu.list <- split(inter.mu.mat, col(inter.mu.mat))
  
  theta.pri.pre.list <- list()
  ntimes = ncol(re$theta0.time)
  for(itime in 1:ntimes){
    inter.cov = inter.mu.cov$cov[[itime]]
    inter.cov[inter.cov == 0] <- 1e-8
    if(det(inter.cov) == 0 | det(inter.cov) < 1e-20){
      inter.diag <- diag(inter.cov)
      inter.diag[inter.diag < 1e-10] <- 1e-8
      inter.cov <- diag(inter.diag)
      inter.inv <- solve(inter.cov)
      theta.pri.pre.list[[itime]] <- inter.inv
    }else{
      inter.inv <- solve(inter.cov)
      theta.pri.pre.list[[itime]] <- inter.inv
    }
  }
  return(list(theta.pri.mu = theta.pri.mu.list,
              theta.pri.pre = theta.pri.pre.list))
}


#' Construct design matrices for estimated exposure quantile functions
#'
#' This function constructs design matrices corresponding to estimated 
#' exposure quantile functions, which will be used subsequently for 
#' fitting negative binomial regression models using the function \code{fit.health.quan.errors}.
#'
#' @param inde a logical value indicating whether exposure quantile processes are 
#' assumed to be independent across time points. 
#' @param re the output from applying the function \code{\link{fit.exposure}}.
#' @param burn_in a number specifying the number of iterations that is discarded. 
#' 
#' @return This function returns a list containing: 
#' (1) a list include vectors containing the mean of the multivariate normal (MVN) 
#' prior assumed for basis coefficients used for modeling exposure quantile functions (\code{theta.pri.mu}), 
#' (2) a list include matrices containing the precision matrix of the MVN prior (\code{theta.pri.pre}). 
#' 
#' @examples
#' # Simulate time-specific means and SDs of exposure quantile functions
#' mean.vec <- rnorm(10, mean = 7.2, sd = 1)
#' sd.vec <- rnorm(10, mean = 1, sd = 0.2)
#' 
#' # Simulate individual-level exposures for all 10 time points
#' x.sim.mat <- do.call(rbind, lapply(1:10, 
#' function(x) rnorm(num.size.time, mean.vec[x], sd.vec[x])))
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
#' @export
exposure.prior <- function(inde, re, burn_in, ...){
  
  ## INPUTs:
  ## inde: a logical value indicating 
  ## whether exposure quantile processes are assumed to be independent 
  ## across time points
  ## if YES, re is a list containing
  ## the output from applying self-defined function 
  ## fit.exposure with argument equals FALSE
  ## to estimate exposure quantile functions 
  ## if NO, a list of length T (# of time points) containing
  ## outputs from applying self-defined function fit.exposure 
  ## with argument inde equals TRUE to 
  ## estimate exposure quantile functions for each time point
  ## burn_in: # of burn-in samples (default is 5000)
  
  if(inde){
    pri.mu.pre <- exposure.prior.inde(re = re, burn_in = burn_in)
  }else{
    pri.mu.pre <- exposure.prior.AR1(re = re, burn_in = burn_in)
  }
  return(pri.mu.pre)
}


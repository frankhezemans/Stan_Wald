######################## Wald ----

# Code below simplified from the Wald-SSEXG version, since we are not estimating the
# parameter A (upper bound of uniform start point distribution).
# Parameterised with drift rate v, starting point z, threshold b (>0),
# and hence B = b-z (>=0) as threshold gap.
# The mean of the Wald distribution is B/v, and shape is B^2.
# Assumes noise of Brownian motion sigma = 1.

# Copyright notice below to acknowledge the original code (ignore comments about parameterisation):

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with: 
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to 
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.
# Following functions use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence 
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A/2 


# random function for single accumulator
rWald <- function(n, v, B) {
  
  # check input and initialise output
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  out <- numeric(n)
  
  # round extremely small but positive drift rates down to zero
  tiny <- v > 0 & v <= 1e-6
  v[tiny] <- 0
  
  # assume that runner with negative drift rate will never terminate, hence RT = Inf
  ok <- v >= 0
  out[!ok] <- Inf
  
  # compute mean and shape parameters of Wald distribution
  mu <- B/v
  lambda <- B^2 
  
  # generate random samples
  out[ok] <- statmod::rinvgauss(n = sum(ok), mean = mu[ok], shape = lambda[ok])
  return(out)
  
}

# density function for single accumulator
dWald <- function(t, v, B) {
  
  # check input and initialise output
  if (length(v)!=length(t)) v <- rep(v,length.out=length(t))
  if (length(B)!=length(t)) B <- rep(B,length.out=length(t))
  out <- numeric(length(t))
  
  # round extremely small but positive drift rates down to zero
  tiny <- v > 0 & v <= 1e-6
  v[tiny] <- 0
  
  # assume that runner with negative drift rate will never terminate, hence density = 0
  ok <- v >= 0
  out[!ok] <- 0
  
  # custom density function - significantly faster than statmod::dinvgauss
  dig <- function(t, v, B) {
    
    lambda <- B^2
    v_zero <- v == 0
    e <- numeric(length(t))
    
    if (any(!v_zero)) {
      mu <- B[!v_zero] / v[!v_zero]
      e[!v_zero] <- -(lambda[!v_zero] / (2*t[!v_zero])) * (t[!v_zero]^2 / mu^2 - 2*t[!v_zero]/mu  + 1)
    }
    if ( any(v_zero) )  e[v_zero] <- -.5*lambda[v_zero] / t[v_zero]
    
    x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
    x[t <= 0] <- 0
    x[x < 0 | is.nan(x)] <- 0
    return(x)
    
  }
  
  # compute densities
  out[ok] <- dig(t = t[ok], v = v[ok], B = B[ok])
  return(out)
  
}

# cumulative density for single accumulator
pWald <- function(t, v, B) {
  
  # check input and initialise output
  if (length(v)!=length(t)) v <- rep(v,length.out=length(t))
  if (length(B)!=length(t)) B <- rep(B,length.out=length(t))
  out <- numeric(length(t))
  
  # round extremely small but positive drift rates down to zero
  tiny <- v > 0 & v <= 1e-6
  v[tiny] <- 0
  
  # assume that runner with negative drift rate will never terminate, hence cumulative density = 0
  ok <- v >= 0
  out[!ok] <- 0
  
  # custom cumulative density function - significantly faster than statmod::pinvgauss
  pig <- function(t, v, B) {
    
    mu <- B/v
    lambda <- B^2
    
    e <- exp(log(2*lambda) - log(mu))
    add <- sqrt(lambda/t) * (1 + t/mu)
    sub <- sqrt(lambda/t) * (1 - t/mu)
    
    p.1 <- 1 - pnorm(add)
    p.2 <- 1 - pnorm(sub)
    x <- exp(e + log(p.1)) + p.2
    
    x[t < 0] <- 0
    x[x < 0 | is.nan(x)] <- 0
    return(x)
    
  }
  
  # compute cumulative densities
  out[ok] <- pig(t = t[ok], v = v[ok], B = B[ok])
  return(out)
  
}

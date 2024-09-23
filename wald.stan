functions {

  /**
    * Convenience function to test if the input value is at or below the lower limit of the
    * inverse Gaussian (a.k.a. Wald) distribution, parameterised in terms of the mean mu and shape
    * lambda.
    *
    * @param x         The input value whose probability needs to be evaluated.
    * @param mu        The mean of the inverse Gaussian distribution.
    * @param lambda    The shape of the inverse Gaussian distribution.
    *
    * @return A boolean value indicating whether input value x is at or below the lower limit of
    *         the distribution (1) or not (0).
    *
    * Acknowledgement: The code below is based on the function `dinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L1-L60
    *
    */
  int is_lower_limit_wald(real x, real mu, real lambda) {
    int is_lower = (
      x < 0 ||
      (x == 0 && mu > 0 && lambda > 0) ||
      (x < mu && lambda == positive_infinity())
    );
    return is_lower;
  }

  /**
    * Convenience function to test if the input value is at or above the upper limit of the
    * inverse Gaussian (a.k.a. Wald) distribution, parameterised in terms of the mean mu and shape
    * lambda.
    *
    * @param x         The input value whose probability needs to be evaluated.
    * @param mu        The mean of the inverse Gaussian distribution.
    * @param lambda    The shape of the inverse Gaussian distribution.
    *
    * @return A boolean value indicating whether input value x is at or above the upper limit of
    *         the distribution (1) or not (0).
    *
    * Acknowledgement: The code below is based on the function `dinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L1-L60
    *
    */
  int is_upper_limit_wald(real x, real mu, real lambda) {
    int is_upper = (
      x == positive_infinity() ||
      (x > 0 && lambda == 0) ||
      (x > mu && (mu == 0 || lambda == positive_infinity()))
    );
    return is_upper;
  }

  /**
    * Convenience function to test if the input value is at a "spike" of the inverse Gaussian
    * (a.k.a. Wald) distribution, parameterised in terms of the mean mu and shape lambda.
    *
    * @param x         The input value whose probability needs to be evaluated.
    * @param mu        The mean of the inverse Gaussian distribution.
    * @param lambda    The shape of the inverse Gaussian distribution.
    *
    * @return A boolean value indicating whether input value x is at a spike (1) or not (0).
    *
    * Acknowledgement: The code below is based on the function `dinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L1-L60
    *
    */
  int is_spike_wald(real x, real mu, real lambda) {
    int is_spike = (
      (x == 0 && lambda == 0) ||
      (x == mu && (mu == 0 || lambda == positive_infinity()))
    );
    return is_spike;
  }

  /**
    * Log probability density function of the inverse Gaussian (a.k.a. Wald) distribution,
    * parameterised in terms of the drift coefficient, diffusion coefficient, and decision
    * threshold.
    * 
    * @param x         The input value whose probability needs to be evaluated.
    * @param drift     The drift coefficient, i.e. the mean rate of evidence accumulation.
    * @param diffusion The diffusion coefficient, i.e. the moment-to-moment variability in
    *                  evidence accumulation.
    * @param boundary  The decision threshold, i.e. the distance from the starting point of 
    *                  evidence accumulation (assumed to be zero) to the absorbing boundary.
    *
    * @return The log probability that random variable X is equal to input value x, given the
    *         parameters of the Wald distribution.
    *
    * @throws If any of the input arguments are NaN, or if the parameters drift, diffusion, and
    *         boundary result in a negative mean or shape of the Wald distribution.
    *
    * Acknowledgement: The code below is based on the function `dinv_gaussian` from the R package
    * `brms`, as well as the function `dinvgauss` from the R package `statmod`.
    * Permanent links to relevant code snippets:
    * https://github.com/paul-buerkner/brms/blob/0c188149495ba488336025e7a741f6b78cc694a6/R/distributions.R#L825-L859
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L1-L79
    *
    */
  real wald_lpdf(real x, real drift, real diffusion, real boundary) {
      
    // initialise local variables
    real drift_scaled; // drift coefficient, scaled by the diffusion coefficient
    real boundary_scaled; // decision threshold, scaled by the diffusion coefficient
    real mu; // mean of the inverse Gaussian distribution
    real lambda; // shape of the inverse Gaussian distribution
    real term1; // constant term of probability density function (does not depend on input x)
    real term2; // term of probability density function that depends on input x
    real lprob; // log probability of input x, given model parameters
        
    // calculate the parameters of the inverse Gaussian distribution
    drift_scaled = drift / diffusion;
    boundary_scaled = boundary / diffusion;

    mu = boundary_scaled / drift_scaled;
    lambda = square(boundary_scaled);

    // check for invalid cases which cannot be evaluated
    if (is_nan(x) || is_nan(mu) || is_nan(lambda) || mu < 0 || lambda < 0) {
      reject("wald_lpdf: unable to evaluate function given input: ",
             "x = ", x, "; drift = ", drift, "; diffusion = ", diffusion,
             "; boundary = ", boundary);
    }

    // check for cases leading to probability of 0, i.e., negative infinity log probability
    if (
      is_lower_limit_wald(x, mu, lambda) ||
      is_upper_limit_wald(x, mu, lambda)
    ) {
      return negative_infinity();
    }

    // check for cases leading to probability of 1, i.e., positive infinity log probability
    if (is_spike_wald(x, mu, lambda)) {
      return positive_infinity();
    }

    // if mu is positive infinity (which occurs e.g. when drift_scaled is zero), the log
    // probability is given by the Levy distribution with parameters location = 0 and
    // scale = lambda
    if (mu == positive_infinity()) {
      lprob = -0.5 * (
        (lambda / x) - log(lambda) + log(2 * pi()) + lmultiply(3, x)
      );
      return lprob;
    }

    // if none of the above cases apply, then calculate the log probability of input x as planned
    term1 = lmultiply(0.5, lambda / (2 * pi()));
    term2 = square(x - mu) / (x * square(mu));

    lprob = term1 - lmultiply(1.5, x) - 0.5 * lambda * term2;

    return lprob;
      
  }

  /**
    * Log cumulative distribution function of the inverse Gaussian (a.k.a. Wald) distribution,
    * parameterised in terms of the drift coefficient, diffusion coefficient, and decision
    * threshold.
    * 
    * @param x         The input value whose probability needs to be evaluated.
    * @param drift     The drift coefficient, i.e. the mean rate of evidence accumulation.
    * @param diffusion The diffusion coefficient, i.e. the moment-to-moment variability in 
    *                  evidence accumulation.
    * @param boundary  The decision threshold, i.e. the distance from the starting point of 
    *                  evidence accumulation (assumed to be zero) to the absorbing boundary.
    *
    * @return The log probability that random variable X will take a value smaller than or equal to
    *         input value x, given the parameters of the Wald distribution.
    *
    * @throws If any of the input arguments are NaN, or if the parameters drift, diffusion, and
    *         boundary result in a negative mean or shape of the Wald distribution.
    *
    * Acknowledgement: The code below is based on the function `pinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L81-L179
    *
    */
  real wald_lcdf(real x, real drift, real diffusion, real boundary) {
      
    // initialise local variables
    real drift_scaled; // drift coefficient, scaled by the diffusion coefficient
    real boundary_scaled; // decision threshold, scaled by the diffusion coefficient
    real mu; // mean of the inverse Gaussian distribution
    real lambda; // shape of the inverse Gaussian distribution
    // terms following Giner & Smyth (2016) page 342
    real x_mu;
    real phi_mu;
    real r;
    real a;
    real b;
    real lcprob; // log cumulative probability
        
    // calculate the parameters of the inverse Gaussian distribution
    drift_scaled = drift / diffusion;
    boundary_scaled = boundary / diffusion;

    mu = boundary_scaled / drift_scaled;
    lambda = square(boundary_scaled);

    // check for invalid cases which cannot be evaluated
    if (is_nan(x) || is_nan(mu) || is_nan(lambda) || mu < 0 || lambda < 0) {
        reject("wald_lcdf: unable to evaluate function given input: ",
                "x = ", x, "; drift = ", drift, "; diffusion = ", diffusion,
                "; boundary = ", boundary);
    }

    // check for cases leading to cumulative probability of zero, i.e., negative infinity log
    // cumulative probability
    if (is_lower_limit_wald(x, mu, lambda)) {
      return negative_infinity();
    }

    // check for cases leading to cumulative probabiliyt of one, i.e., zero log cumulative
    // probability
    if (
      is_upper_limit_wald(x, mu, lambda) ||
      is_spike_wald(x, mu, lambda)
    ) {
      return 0.0;
    }
    
    // if mu is positive infinity (which occurs e.g. when drift_scaled is zero), the log cumulative
    // probability is given by the log complementary cumulative probability of the Chi-Square
    // distribution with parameter degrees of freedom = 1
    if (mu == positive_infinity()) {
        return chi_square_lccdf((lambda / x) | 1);
    }

    // the coefficient of variation of the inverse Gaussian distribution is equal to the square
    // root of mu over lambda. if mu over lambda is tiny, the log cumulative probability is given
    // by the log cumulative probability of the Gamma distribution with parameters
    // shape = lambda / mu and rate (inverse scale) lambda / mu^2
    if ((mu / lambda) < 1e-14) {
        return gamma_lcdf(x | (lambda / mu), (lambda / square(mu)));
    }

    // if none of the above cases apply, then calculate the log cumulative probability of input x
    // as planned
    x_mu = x / mu;
    phi_mu = mu / lambda;
    r = sqrt(x_mu * phi_mu);
    // NB using std_normal_lcdf() instead of Phi() or Phi_approx(), as the former seems to be much
    // more accurate (contrary to the Stan function reference, which still recommends Phi_approx())
    // https://github.com/stan-dev/math/issues/2470
    a = std_normal_lcdf((x_mu - 1) / r);
    b = 2 / phi_mu + std_normal_lcdf(-(x_mu + 1) / r);
    lcprob = a + log1p_exp(b - a);

    return lcprob;

  }

  /**
    * Log complementary cumulative distribution function (a.k.a. survival function) of the inverse
    * Gaussian (a.k.a. Wald) distribution, parameterised in terms of the drift coefficient,
    * diffusion coefficient, and decision threshold.
    * 
    * @param x         The input value whose probability needs to be evaluated.
    * @param drift     The drift coefficient, i.e. the mean rate of evidence accumulation.
    * @param diffusion The diffusion coefficient, i.e. the moment-to-moment variability in 
    *                  evidence accumulation.
    * @param boundary  The decision threshold, i.e. the distance from the starting point of 
    *                  evidence accumulation (assumed to be zero) to the absorbing boundary.
    *
    * @return The log probability that random variable X will take a value greater than input value
    *         x, given the parameters of the Wald distribution.
    *
    * @throws If any of the input arguments are NaN, or if the parameters drift, diffusion, and
    *         boundary result in a negative mean or shape of the Wald distribution.
    *
    * Acknowledgement: The code below is based on the function `pinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L81-L179
    *
    */
  real wald_lccdf(real x, real drift, real diffusion, real boundary) {
      
    // initialise local variables
    real drift_scaled; // drift coefficient, scaled by the diffusion coefficient
    real boundary_scaled; // decision threshold, scaled by the diffusion coefficient
    real mu; // mean of the inverse Gaussian distribution
    real lambda; // shape of the inverse Gaussian distribution
    // terms following Giner & Smyth (2016) page 342
    real x_mu;
    real phi_mu;
    real x_phi_mu;
    real r;
    real a;
    real b;
    real lccprob; // log complementary cumulative probability
        
    // calculate the parameters of the inverse Gaussian distribution
    drift_scaled = drift / diffusion;
    boundary_scaled = boundary / diffusion;

    mu = boundary_scaled / drift_scaled;
    lambda = square(boundary_scaled);

    // check for invalid cases which cannot be evaluated
    if (is_nan(x) || is_nan(mu) || is_nan(lambda) || mu < 0 || lambda < 0) {
        reject("wald_lccdf: unable to evaluate function given input: ",
                "x = ", x, "; drift = ", drift, "; diffusion = ", diffusion,
                "; boundary = ", boundary);
    }

    // check for cases leading to complementary cumulative probability of zero, i.e., negative
    // infinity log complementary cumulative probability
    if (is_upper_limit_wald(x, mu, lambda)) {
      return negative_infinity();
    }

    // check for cases leading to complementary cumulative probabiliyt of one, i.e., zero log
    // complementary cumulative probability
    if (
      is_lower_limit_wald(x, mu, lambda) ||
      is_spike_wald(x, mu, lambda)
    ) {
      return 0.0;
    }
    
    // if mu is positive infinity (which occurs e.g. when drift_scaled is zero), the log
    // complementary cumulative probability is given by the log cumulative probability of the
    // Chi-Square distribution with parameter degrees of freedom = 1
    if (mu == positive_infinity()) {
        return chi_square_lcdf((lambda / x) | 1);
    }

    // the coefficient of variation of the inverse Gaussian distribution is equal to the square
    // root of mu over lambda. if mu over lambda is tiny, the log complementary cumulative
    // probability is given by the log complementary cumulative probability of the Gamma distribution
    // with parameters shape = lambda / mu and rate (inverse scale) lambda / mu^2
    if ((mu / lambda) < 1e-14) {
        return gamma_lccdf(x | (lambda / mu), (lambda / square(mu)));
    }

    // if none of the above cases apply, then calculate the log complementary cumulative probability
    // of input x as planned
    x_mu = x / mu;
    phi_mu = mu / lambda;
    x_phi_mu = x_mu / (2 * phi_mu);
    // if input value x is extremely large, use asymptotic expression given by Giner & Smyth (2016)
    // otherwise, use conventional calculation
    if (x_mu > 1e6 || x_phi_mu > 5e5) {
      lccprob = (1 / phi_mu) - lmultiply(0.5, pi()) - log(2 * phi_mu) -
        1.5 * log1p(x_phi_mu) - x_phi_mu;
    } else {
      r = sqrt(x_mu * phi_mu);
      // NB using std_normal_lcdf(-((x_mu - 1) / r)) instead of std_normal_lccdf((x_mu - 1) / r) to
      // compute term a. Theoretically these should be identical, but in practice the former might
      // be more accurate:
      // https://github.com/stan-dev/math/issues/1985
      a = std_normal_lcdf(-((x_mu - 1) / r));
      b = 2 / phi_mu + std_normal_lcdf(-(x_mu + 1) / r);
      lccprob = a + log1m_exp(b - a);
    }

    return lccprob;

  }

  /**
    * Generate a random variate from the inverse Gaussian (a.k.a. Wald) distribution, parameterised
    * in terms of the drift coefficient, diffusion coefficient, and decision threshold.
    *
    * @param drift     The drift coefficient, i.e. the mean rate of evidence accumulation.
    * @param diffusion The diffusion coefficient, i.e. the moment-to-moment variability in 
    *                  evidence accumulation.
    * @param boundary  The decision threshold, i.e. the distance from the starting point of 
    *                  evidence accumulation (assumed to be zero) to the absorbing boundary.
    *
    * @return A random variate from the inverse Gaussian distribution parameterised with the given
    *         values.
    *
    * @throws If any of the input arguments are NaN, or if the parameters drift, diffusion, and
    *         boundary result in a negative mean or shape of the Wald distribution.
    *
    * Acknowledgement: The code below is based on the function `rinvgauss` from the R package
    * `statmod`.
    * Permanent link to relevant code snippet:
    * https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R#L181-L249
    *
    */
  real wald_rng(real drift, real diffusion, real boundary) {

    // initialise local variables
    real drift_scaled; // drift coefficient, scaled by the diffusion coefficient
    real boundary_scaled; // decision threshold, scaled by the diffusion coefficient
    real mu; // mean of the inverse Gaussian distribution
    real lambda; // shape of the inverse Gaussian distribution
    // some terms following statmod::rinvgauss
    real y;
    real y_phi;
    real x_1;
    real unif_sim;
    real x_sim; // output: random variate from inverse Gaussian distribution

    // calculate the parameters of the inverse Gaussian distribution
    drift_scaled = drift / diffusion;
    boundary_scaled = boundary / diffusion;

    mu = boundary_scaled / drift_scaled;
    lambda = square(boundary_scaled);

    // check for invalid cases which cannot be evaluated
    if (is_nan(mu) || is_nan(lambda) || mu < 0 || lambda < 0) {
        reject("wald_rng: unable to evaluate function given input: ",
                "; drift = ", drift, "; diffusion = ", diffusion,
                "; boundary = ", boundary);
    }

    if (mu == 0 || lambda == positive_infinity()) {
      return not_a_number();
    }

    if (lambda == 0) {
      return 0.0;
    }

    // special case: mu positive infinity but lambda okay
    // NB we already established that lambda < Inf with previous if-statement
    if (mu == positive_infinity() && lambda > 0) {
      return ((1 / square(std_normal_rng())) * lambda);
    }

    // if the above cases don't apply, then we know that mu and lambda are
    // strictly positive finite values, so we can continue as planned
    y = square(std_normal_rng());
    y_phi = (y / lambda) * mu;

    if (y_phi > 5e5) {
      x_1 = 1 / y_phi;
    } else {
      x_1 = 1 + (y_phi / 2) * (1 - sqrt(1 + (4 / y_phi)));
    }

    unif_sim = uniform_rng(0, 1);
    if (unif_sim < (1 / (1 + x_1))) {
      x_sim = mu * x_1;
    } else {
      x_sim = mu / x_1;
    }

    return x_sim;

  }

}

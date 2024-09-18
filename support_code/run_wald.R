# define function for evaluating all the different Wald probability functions
run_wald <- function(x, v, s, B) {

  # get parameters of inverse Gaussian distribution, based on the parameters
  # of the diffusion process
  v_scaled <- v / s
  B_scaled <- B / s
  mu <- B_scaled / v_scaled
  lambda <- B_scaled^2

  # log probability density function
  out_lpdf <- data.frame(
    x = x,
    fun = "lpdf",
    statmod = statmod::dinvgauss(x = x, mean = mu, shape = lambda, log = TRUE),
    Stan_Wald = wald_lpdf(x = x, drift = v, diffusion = s, boundary = B),
    brms = brms::dinv_gaussian(x = x, mu = mu, shape = lambda, log = TRUE),
    DMC = log(dWald(t = x, v = v_scaled, B = B_scaled))
  )

  # log cumulative distribution function
  out_lcdf <- data.frame(
    x = x,
    fun = "lcdf",
    statmod = statmod::pinvgauss(
      q = x, mean = mu, shape = lambda, lower.tail = TRUE, log.p = TRUE
    ),
    Stan_Wald = wald_lcdf(x = x, drift = v, diffusion = s, boundary = B),
    brms = brms::pinv_gaussian(
      q = x, mu = mu, shape = lambda, lower.tail = TRUE, log.p = TRUE
    ),
    DMC = log(pWald(t = x, v = v_scaled, B = B_scaled))
  )

  # log complementary cumulative distribution function
  out_lccdf <- data.frame(
    x = x,
    fun = "lccdf",
    statmod = statmod::pinvgauss(
      q = x, mean = mu, shape = lambda, lower.tail = FALSE, log.p = TRUE
    ),
    Stan_Wald = wald_lccdf(x = x, drift = v, diffusion = s, boundary = B),
    brms = brms::pinv_gaussian(
      q = x, mu = mu, shape = lambda, lower.tail = FALSE, log.p = TRUE
    ),
    DMC = log1p(-pWald(t = x, v = v_scaled, B = B_scaled))
  )

  # bind output together, and put in long format with respect to the
  # implementation used
  out <- dplyr::bind_rows(out_lpdf, out_lcdf, out_lccdf) |>
    tidyr::pivot_longer(
      cols = tidyr::all_of(c(
        "statmod", "Stan_Wald", "brms", "DMC"
      )),
      names_to = "method",
      values_to = "value"
    ) |>
    dplyr::select(-x)

  return(out)

}

# Comparing implementations of Wald probability functions
Frank H. Hezemans
2024-09-18

This notebook compares my Stan implementation of the Wald distribution
(`Stan_Wald`) with several R implementations that were released with the
following packages / toolboxes: `statmod` (Giner & Smyth, 2016), `brms`
(Bürkner, 2017), and Dynamic Models of Choice (`DMC`; Heathcote et al.,
2019). I will treat the `statmod` implementation as the reference point;
I expect my Stan implementation to perform identically to the `statmod`
implementation, and I expect the `brms` and `DMC` implementations to
potentially be less accurate, especially for the log (complementary)
cumulative distribution functions.

First, I’ll load the required R packages as well as some support code
that I’ve provided in the folder `support_code/`, and then I’ll compile
my Stan code and load (“expose”) the functions into the R workspace.

``` r
# load library for various data curating and plotting tools
library(tidyverse)
# load "reference point" of probability functions for Wald distribution
library(statmod)
# load brms package, which also contains an implementation of the Wald
# distribution
library(brms)
# load various convenience functions, as well as Wald probability functions
# that were released with the DMC toolbox
purrr::walk(
  .x = list.files(path = "support_code", full.names = TRUE),
  .f = source 
)
# use the rstan package to compile my Stan code and "expose" the functions into
# the R workspace
library(rstan)
rstan::expose_stan_functions(
  stanmodel = rstan::stanc(file = "wald.stan")
)
```

Next, I’ll define the parameter values for which we want to evaluate the
Wald probability functions, as well as a grid of input values. I’m using
an equally-spaced grid from 0 to 6, and I’m deliberately including one
negative value (-0.1) as an extra test.

I’ve chosen the following parameter sets (with the resulting mu and
lambda values provided for your reference):

    # A tibble: 3 × 5
          v     s     B    mu lambda
      <dbl> <dbl> <dbl> <dbl>  <dbl>
    1   4       1     1  0.25      1
    2   2.5     1     1  0.4       1
    3   2.5     1     3  1.2       9

I then use my convenience function `support_code/run_wald` to apply the
log probability density function, log cumulative distribution function,
and log complementary cumulative distribution function to each of the
parameter sets, for each of the four implementations under
consideration.

For each Figure, the input values are plotted along the x-axis; the
computed log-probabilities are plotted along the y-axis; the different
implementations are presented in separate panels; and the parameter sets
(v, s, B) are represented by different colour hues.

``` r
# declare grid of parameter values for the Wald distribution, as well as input
# values for which we want to evaluate the probabilities.
wald_params <- tibble::as_tibble(expand.grid(
  x = c(          # input values
    -0.1, seq(from = 0, to = 6, by = 0.1)
  ),
  v = c(2.5, 4),   # drift coefficient
  s = 1,          # diffusion coefficients
  B = c(1, 3)     # boundary
)) |>
  # remove one parameter set to avoid overplotting
  dplyr::filter(!(v == 4 & B == 3))

# now calculate the various log probabilities for each of the above parameter
# values, using my convenience function `run_wald`
wald_result <- wald_params |>
  dplyr::mutate(
    out = purrr::pmap(
      .l = list(.data[["x"]], .data[["v"]], .data[["s"]], .data[["B"]]),
      .f = run_wald
    )
  ) |>
  tidyr::unnest(out) |>
  # turn into one nested data frame for each probability function
  # (lpdf / lcdf / lccdf)
  tidyr::nest(data = -fun)
```

Next, I’ll use my convenience function `support_code/plot_wald` to plot
the results. I’ll discuss the results spearately for each of the three
probability functions, comparing the accuracy of the results across
parameter sets and implementations.

``` r
# use my convenience function `plot_wald` to plot the results for each
# probability function (lpdf / lcdf / lccdf)
wald_result_plots <- wald_result |>
  dplyr::mutate(
    plots = purrr::map2(
      .f = plot_wald,
      .x = data,
      .y = fun
    )
  ) |>
  dplyr::pull(plots)
```

## Log probability density function

``` r
wald_result_plots[[1]]
```

![](comparing_wald_files/figure-commonmark/show_lpdf_plot-1.png)

The results of the LPDF look fine for all implementations and parameter
sets. The only thing to note is that the `brms` implementation returns
`NaN` for negative input values, whereas the other implementations
return negative infinity.

## Log cumulative distribution function

``` r
wald_result_plots[[2]]
```

![](comparing_wald_files/figure-commonmark/show_lcdf_plot-1.png)

At first glance, the results of the LCDF also look fine for all
implementations, but in fact the `brms` and `DMC` implementations become
inaccurate for relatively extreme input values. This is best illustrated
by focusing on the parameter set with a drift rate of 4, and zooming in
on relatively high x values:

``` r
wald_lcdf_zoomed_data <- wald_result |>
  dplyr::filter(fun == "lcdf") |>
  dplyr::select(data) |>
  tidyr::unnest(data) |>
  dplyr::filter(
    x > 4.2 & v == 4 & B == 1
  )

wald_lcdf_zoomed_plot <- plot_wald(
  data = wald_lcdf_zoomed_data,
  fun = "lcdf",
  facet_scales = "free_y",
  xlim = c(4, 6.5)
) +
  ggplot2::scale_y_continuous(
    labels = scales::label_scientific(digits = 3)
  ) +
  ggplot2::theme(
    legend.position = "null"
  )

wald_lcdf_zoomed_plot
```

![](comparing_wald_files/figure-commonmark/show_lcdf_plot_zoomed-1.png)

Here, we can start to appreciate that the `brms` and `DMC`
implementations begin to overflow to a probability of 1 (i.e., log
probability of 0) for input values greater than approximately 4.5,
whereas the `statmod` and `Stan_Wald` implementations still return
accurate (non-zero) log probabilities. We can highlight this further by
zooming in even more:

``` r
wald_lcdf_zoomed_more_plot <- plot_wald(
  data = wald_lcdf_zoomed_data |>
    dplyr::filter(x > 5),
  fun = "lcdf",
  facet_scales = "free_y",
  xlim = c(5, 6.5)
) +
  ggplot2::scale_y_continuous(
    labels = scales::label_scientific(digits = 3)
  ) +
  ggplot2::theme(
    legend.position = "null"
  )

wald_lcdf_zoomed_more_plot
```

![](comparing_wald_files/figure-commonmark/show_lcdf_plot_zoomed_alt-1.png)

Again, the `brms` and `DMC` implementations are unable to evaluate the
LCDF in the far tails of the distribution, whereas the `statmod` and
`Stan_Wald` implementations remain accurate.

## Log complementary cumulative distribution function

``` r
wald_result_plots[[3]]
```

![](comparing_wald_files/figure-commonmark/show_lccdf_plot-1.png)

For the LCCDF, we can clearly see that the `brms` and `DMC`
implementations underflow to zero probability (i.e., negative infinity
log probability), specifically for the parameter set including a drift
rate of 4 (darkest colour hues).

## References

Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using
Stan. *Journal of Statistical Software*, *80*(1), 1-28.

Giner, G., & Smyth, G. K. (2016). statmod: Probability Calculations for
the Inverse Gaussian Distribution. *The R Journal*, *8*(1), 339-351.

Heathcote, A., Lin, Y. S., Reynolds, A., Strickland, L., Gretton, M., &
Matzke, D. (2019). Dynamic models of choice. *Behavior research
methods*, *51*, 961-985.

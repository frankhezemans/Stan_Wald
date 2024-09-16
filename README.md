# Stan_Wald: Stan probability functions for the Wald distribution

<!-- badges: start -->
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

The **Wald distribution** (also known as the [inverse Gaussian distribution](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)) is a [positively skewed](https://en.wikipedia.org/wiki/Skewness) [unimodal](https://en.wikipedia.org/wiki/Unimodality) distribution commonly used to model [response time data](https://en.wikipedia.org/wiki/Mental_chronometry) (Anders, Alario, & Van Maanen, 2016; Schwarz, 2001). Parameterised as a [single-boundary diffusion process](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#Relationship_with_Brownian_motion), the Wald distribution allows for inference on the latent psychological mechanisms underlying observed data (Heathcote & Matzke, 2022; Matzke & Wagenmakers, 2009; Tillman, Van Zandt, & Logan, 2020).

**Stan** is a state-of-the-art platform for statistical modeling, but it lacks native support for the Wald distribution. This repository provides **Stan implementations of the Wald distribution's probability functions**, including:
- Log probability density function (`wald_lpdf`)
- Log cumulative distribution function (`wald_lcdf`)
- Log complementary cumulative distribution function (`wald_lccdf`)
- Random deviate generator (`wald_rng`)

These functions can be incorporated into any Stan model as [user-defined probability functions](https://mc-stan.org/docs/stan-users-guide/user-functions.html#user-defined-probability-functions). The functions comply with recent Stan language updates (September 2023[^1]) and are [clearly documented](https://mc-stan.org/docs/stan-users-guide/user-functions.html#documenting-functions.section).

[^1]: For details, see the [list of features](https://mc-stan.org/docs/reference-manual/removals.html) that were removed in Stan version 2.33.

## Reference implementation

These Stan functions are, for the most part, direct translations of R code released with the `statmod` package (Giner & Smyth, 2016), specifically the file [`invgauss.R`](https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R). Any changes relative to the `statmod` code were made for enhanced readability and to comply with the Stan language syntax.

The `statmod` implementation was chosen as the reference point as it has been shown to be highly accurate for an extremely wide range of parameter values, and it robustly handles edge cases of parameter values that would produce undefined or non-finite probabilities (Giner & Smyth, 2016). In particular, `statmod`'s implementations of the Wald LCDF and LCCDF are significantly more accurate than other R implementations, which suffer from numerical problems such as [arithmetic under- or overflow](https://en.wikipedia.org/wiki/Integer_overflow) and [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation).

## Usage considerations

### Parameterisation

The Wald distribution represents the first passage time distribution of a Wiener diffusion process with a single absorbing boundary. In the context of psychological experiments involving speeded decision-making, we interpret the diffusion process as a process of accumulating noisy partial information ("evidence") over time. Following this interpretation, the Wald probability functions consist of the following parameters:
- drift coefficient, _v_, representing the mean gain of evidence per unit time
- diffusion coefficient, _s_, representing the variance in the gain of evidence per unit time
- boundary, _B_, representing the total amount of evidence needed to terminate the evidence accumulation process (i.e., reach a decision)

These parameters are used to determine the mean _mu_ and shape _lambda_ of the inverse Gaussian distribution, which, in turn, are used to compute the requested log probability:

```
v_scaled = v / s
B_scaled = B / s
mu = B_scaled / v_scaled
lambda = (B_scaled)^2
```

From the above identities, it becomes evident that one of _v_, _s_, or _B_ must be fixed to a constant value so that the two remaining parameters are identifiable (Donkin, Brown & Heathcote, 2009; Van Maanen & Miletić, 2021). Conventionally, the diffusion coefficient _s_ is fixed to 1 (Schwarz, 2001; Heathcote, 2004), but ultimately this is the user's choice. 

### Vectorisation

Unlike the probability functions in R, these Stan functions are not vectorised: they return a single probability value for a single given input value[^2]. Thus, a for-loop must be used to apply these functions to multiple data points. This should not affect performance significantly, as loops are optimised in Stan's C++ back-end, such that for-loops may actually perform similarly to vectorised code in practice (at least when _user-defined_ as opposed to _built-in_ probability functions are involved).

[^2]: In Stan, user-defined functions may only be used as probability functions in distribution statements if they return a real (scalar) value as output. Thus, a vectorised probability function in Stan would return a single value (namely, the _summed_ log probability) for a given vector of input values.

## Table of contents

- `wald.stan` contains Stan functions corresponding to the Wald probability functions (`_lpdf`, `_lcdf`, `_lccdf`, and `_rng`), as well as internal functions used to check for edge cases of extreme parameter values.

## Citation

If you use this code in your work, please cite the repository:

Hezemans (2024). TEXT HERE

A BibTeX entry for LaTeX users is as follows:

```
TEXT HERE
```

Please also consider citing the R package `statmod`:

Giner, G., & Smyth, G. K. (2016). statmod: Probability Calculations for the Inverse Gaussian Distribution. _The R Journal_, _8_(1), 339-351.

```
@article{statmod2016,
  author = {G{\"o}knur Giner and Gordon K. Smyth},
  title = {{statmod: Probability Calculations for the Inverse Gaussian Distribution}},
  year = {2016},
  journal = {{The R Journal}},
  doi = {10.32614/RJ-2016-024},
  url = {https://doi.org/10.32614/RJ-2016-024},
  pages = {339--351},
  volume = {8},
  number = {1}
}
```

You should also [cite Stan](https://mc-stan.org/users/citations/) according to the specific Stan packages used.

## License

This repository is licensed under [**GPL-2**](https://choosealicense.com/licenses/gpl-2.0/) due to its derivation from the `statmod` R package, which was licensed under GPL-2 | [GPL-3](https://choosealicense.com/licenses/gpl-3.0/).

## References

Anders, R., Alario, F., & Van Maanen, L. (2016). The shifted wald distribution for response time data analysis. _Psychological Methods_, _21_(3), 309.

Donkin, C., Brown, S. D., & Heathcote, A. (2009). The overconstraint of response time models: Rethinking the scaling problem. _Psychonomic Bulletin & Review_, _16_(6), 1129–1135.

Giner, G., & Smyth, G. K. (2016). statmod: Probability Calculations for the Inverse Gaussian Distribution. _The R Journal_, _8_(1), 339-351.

Heathcote, A. (2004). Fitting wald and ex-wald distributions to response time data: An example using functions for the s-plus package. _Behavior Research Methods, Instruments, & Computers_, _36_(4), 678–694.

Heathcote, A., & Matzke, D. (2022). Winner takes all! What are race models, and why and how should psychologists use them?. _Current Directions in Psychological Science_, _31_(5), 383-394.

Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. _Psychonomic Bulletin & Review_, _16_, 798-817.

Schwarz, W. (2001). The ex-wald distribution as a descriptive model of response times. _Behavior Research Methods, Instruments, & Computers_, _33_(4), 457–469.

Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: The racing diffusion model of speeded decision making. _Psychonomic Bulletin & Review_, _27_(5), 911-936.
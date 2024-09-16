# Stan_Wald: Stan probability functions for the Wald distribution

<!-- badges: start -->
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

The Wald distribution (also known as [inverse Gaussian distribution](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)) is a [positively skewed](https://en.wikipedia.org/wiki/Skewness) [unimodal](https://en.wikipedia.org/wiki/Unimodality) distribution that is commonly used to model [response time (RT) data](https://en.wikipedia.org/wiki/Mental_chronometry) (REFS). It can be parameterised in terms of a [diffusion process with a single boundary](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#Relationship_with_Brownian_motion), which enables inference on the latent psychological mechanisms underlying observed RT data. Stan is a state-of-the-art platform for statistical modeling and high-performance statistical computation, but the Wald distribution is not readily available in the [Stan programming language](https://mc-stan.org/docs/reference-manual/).

This repository provides probability functions for the Wald distribution in Stan code. These functions can be incorporated into any `.stan` model file as [user-defined probability functions](https://mc-stan.org/docs/stan-users-guide/user-functions.html#user-defined-probability-functions). The functions have been written to comply with recent (September 2023) [updates](https://mc-stan.org/docs/reference-manual/removals.html) in the Stan language syntax, and are [clearly documented](https://mc-stan.org/docs/stan-users-guide/user-functions.html#documenting-functions.section).

These Stan functions are essentially translations of R code released with the `statmod` package, specifically [the file `invgauss.R`](https://github.com/cran/statmod/blob/f85e32011346fb75d2b967cf2aff1f2e01a10ba8/R/invgauss.R). Any changes compared to the `statmod` code were made for enhanced readability and to comply with the Stan language syntax.

The `statmod` implementation was chosen as the reference point as it has been shown to be highly accurate for an extremely wide range of parameter values, and carefully handles special cases of parameter values that are known to produce non-finite probabilities. In particular, `statmod`'s implementation of the Wald (complementary) cumulative distribution function is significantly more accurate than other R implementations, which suffer from numerical problems such as [arithmetic under- or overflow](https://en.wikipedia.org/wiki/Integer_overflow) and [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation).

Note that these Stan functions are not vectorised --- they return a single probability value for a single given input value. This is in contrast to the behaviour of probability functions in R (including those in the `statmod` package), which return a _vector_ of probability values for a corresponding vector of input values. The reason for this is that in the Stan language, user-defined functions may only be used as probability functions in distribution statements if they return a single value as output (i.e., a `real` return type). Thus, a "vectorised" Stan probability function would return a single value --- namely, the _summed_ log probability --- for a given vector of input values. However, when working with user-defined probability functions, a for-loop may in practice actually perform similarly to vectorised code in Stan, because loops are heavily optimised by Stan's C++ back-end.

## Table of contents

- `wald.stan` contains Stan functions corresponding to the following probability functions of the Wald distribution:
  * log probability density function (LPDF)
  * log cumulative distribution function (LCDF)
  * log complementary cumulative distribution function (LCCDF)
  
  In addition, a function is provided to generate a random deviate (RNG) from the Wald distribution.

  Lastly, several convenience functions are included to test for edge cases involving extreme parameter values that are known to result in undefined or non-finite probabilities. These functions are used internally by the probability functions to prevent calculating log probabilities for these invalid cases, which could otherwise result in undefined (`NaN`) output.

## Citation

If you use this code in your work, please cite the repository as follows:

Hezemans (2024). TEXT HERE

A BibTeX entry for LaTeX users is as follows:

```
TEXT HERE
```

Please also consider citing the R package `statmod`, which served as the basis for the current repository's code:

Giner, G., & Smyth, G. K. (2016). statmod: Probability Calculations for the Inverse Gaussian Distribution. _The R Journal_, 8(1), 339-351.

A BibTeX entry for LaTeX users is as follows:

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

Naturally, you should also [cite Stan](https://mc-stan.org/users/citations/) according to the specific packages used.

## License

The key file in this repository, `wald.stan` was largely based on code from the R package `statmod`, which was released under two versions of the GNU General Public License ([GPL-2](https://choosealicense.com/licenses/gpl-2.0/) | [GPL-3](https://choosealicense.com/licenses/gpl-3.0/)). These licenses have "copyleft" provisions requiring derivative work to be distributed under a GPL-compatible license. Therefore, the code in this repository is released under a GPL-2 license.
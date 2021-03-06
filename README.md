[![Build Status](https://github.com/PabRod/sleepR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/PabRod/sleepR/actions)
[![codecov](https://codecov.io/gh/PabRod/sleepR/branch/master/graph/badge.svg?token=tVmhqzqNuM)](https://codecov.io/gh/PabRod/sleepR)
[![codecov](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# SleepR
Simulation of different sleep-wake deterministic models. See references below.

## Installing

### Latest stable version
Type `devtools::install_github("PabRod/sleepR", ref = "master")` in your `R` command console.

### Latest version
Type `devtools::install_github("PabRod/sleepR", ref = "develop")` in your `R` command console.

### Running the tests
The integrity of this package can be checked by running the battery of tests available at `./tests`.

## Author
Pablo Rodríguez-Sánchez [![](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-2855-940X)

## Examples of usage
Examples available at the vignettes for:

- [Phillips and Robinson's model](inst/doc/philrob.pdf)
- [Kronauer's model](inst/doc/kronauer.pdf)
- [Strogatz's model](inst/doc/strogatz.pdf)

## References
- Strogatz, S. H. (1987). Human sleep and circadian rhythms: a simple model based on two coupled oscillators. Journal of Mathematical Biology, 25(3), 327–347. http://doi.org/10.1007/BF00276440
- Jewett ME, Kronauer RE. (1998). Refinement of Limit Cycle Oscillator Model of the Effects of Light on the Human Circadian Pacemaker. J Theor Biol. 1998 Jun 21;192(4):455–65. https://www.sciencedirect.com/science/article/pii/S0022519398906671
- Phillips AJK, Robinson PA. (2007). A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System. J Biol Rhythms. 2007 Apr 29;22(2):167–79. http://journals.sagepub.com/doi/10.1177/0748730406297512
- Skeldon AC, Dijk D-J, Derks G. (2014) Mathematical Models for Sleep-Wake Dynamics: Comparison of the Two-Process Model and a Mutual Inhibition Neuronal Model. PLoS One. 2014 Aug 1;9(8):e103877. Available from: http://dx.plos.org/10.1371/journal.pone.0103877

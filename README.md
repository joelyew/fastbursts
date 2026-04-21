# fastbursts

The `fastbursts` package provides a performance-optimized implementation of Kleinberg's (2002) burst detection algorithm, utilizing a linear-time dynamic programming approach for efficient analysis of large event sequences.


## Why use fastbursts?

The standard `bursts` package in R has an $O(n^2)$ time and memory complexity because it stores the entire state path at each step. This makes it unusable for large event streams (e.g., 50,000+ events).

`fastbursts` implements the same algorithm using a back-pointer (Viterbi) approach, reducing the complexity to $O(nk)$. In benchmarks, `fastbursts` is **30x-100x faster** for moderate datasets and maintains 100% output parity with the original package.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("joelyew/fastbursts")
```

## Usage

```r
library(fastbursts)

# Example event offsets
offsets <- c(seq(0, 400, 100), seq(410, 450, 5), seq(451, 470, 2),
             seq(480, 600, 5), 700, seq(710, 800, 5), 900, 1000)

# Run burst detection
result <- kleinberg(offsets)

# Plot results
plot(result)
```

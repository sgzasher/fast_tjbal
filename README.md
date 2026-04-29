# fast.tjbal

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**R** package for trajectory balancing — a kernel-based reweighting method for causal inference with time-series cross-sectional (panel) data.

This is a fork of [`tjbal`](https://github.com/xuyiqing/tjbal) (Hazlett & Xu) that replaces the `kbal` dependency with a self-contained solver and adds a Nystrom approximation for scalability to large panels.

## What's different from `tjbal`

The original `tjbal` delegates kernel balancing to the [`kbal`](https://github.com/csterbenz1/KBAL) package, which solves an exact balancing problem using a truncated eigendecomposition of the full N x N Gram matrix.

`fast.tjbal` reformulates the balancing objective as a **kernel maximum mean discrepancy (MMD)** minimization problem with entropy regularization:

```
w* = argmin   w' K_CC w  -  2 k_TC' w  +  eta * sum(w_i log w_i)
      s.t.    sum(w_i) = 1,  w_i >= 0
```

where `K_CC` is the kernel Gram matrix among controls and `k_TC` is the kernel mean map from treated units to controls. The entropy term ensures strict convexity (unique solution), keeps all weights positive, and provides dimension-free weight stability via Pinsker's inequality.

The KKT conditions reduce to a **k-dimensional fixed-point system** `G(z) = 0`, solved by Newton's method with backtracking line search. This replaces the N-dimensional primal problem with a system whose dimension equals the rank of the kernel approximation.

Two computational paths are available:

- **Exact path** (`nystrom = FALSE`): eigendecompose the full Gram matrix, automatically select the rank `k` that minimizes a bias bound over the residual eigenspace.
- **Nystrom path** (`nystrom = TRUE`, default for N > 1000): sample `m` landmark points, form a Nystrom feature matrix `Phi` such that `K ~ Phi Phi'`, and solve the same fixed-point system without ever forming the N x N kernel matrix. Cost is O(N m^2) instead of O(N^3).

## Installation

```r
install.packages("devtools", repos = "http://cran.us.r-project.org") # if not already installed
devtools::install_github("sgzasher/fast_tjbal")
```

No external kernel-balancing package is required — all balancing is handled internally.

## Usage

The interface is the same as `tjbal`:

```r
library(fast.tjbal)

out <- tjbal(Y ~ D, data = panel, index = c("unit", "time"))

# Or with covariates and explicit options:
out <- tjbal(Y ~ D + X1 + X2, data = panel, index = c("unit", "time"),
             estimator = "kernel",  # "kernel" (default) or "mean"
             eta = 0.001,           # entropy regularization strength
             nystrom = NULL,        # NULL = auto (TRUE if N > 1000)
             nystrom_m = NULL,      # number of landmarks (default: 3*sqrt(N))
             vce = "jackknife")     # "jackknife", "bootstrap", or "fixed.weights"

print(out)
plot(out, type = "gap")
plot(out, type = "counterfactual")
```

Staggered adoption designs are supported — the estimator loops over treatment-timing subgroups and aggregates ATT estimates in event time.

## Reference

Asher, Sam, Chad Hazlett, and Yiqing Xu. "Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data." Available at SSRN: <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3214231>.

## Authors

Chad Hazlett, Yiqing Xu, Samuel Asher

## License

MIT


<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatbiome

A collection of functions for spatial analysis of microbiome data.
Currently contains ARShift, a random shift test for bivariate point
patterns.

## Installation

You can install the development version of spatbiome from
[GitHub](https://github.com/) with:

``` r

# Install a package from GitHub
pak::pak("asmita112358/spatbiome")
#> ℹ Loading metadata database✔ Loading metadata database ... done
#>  
#> → Package library at '/Library/Frameworks/R.framework/Versions/4.6/Resources/library'.
#> ℹ No downloads are needed
#> ✔ 1 pkg + 40 deps: kept 40 [4.2s]
library(spatbiome)
```

## Data generating functions and simulation scripts

The package includes the data-generating functions in scripts/datagen.R.
To locate this script and the results_in_paper directory in your
installed copy of spatbiome, run

``` r
system.file("results_in_paper", package = "spatbiome")
#> [1] "/Library/Frameworks/R.framework/Versions/4.6/Resources/library/spatbiome/results_in_paper"
system.file("scripts", package = "spatbiome")
#> [1] "/Library/Frameworks/R.framework/Versions/4.6/Resources/library/spatbiome/scripts"
```

These commands return the full paths to the files and directory. Open
datagen.R in RStudio with:

``` r
file.edit(system.file("scripts", "datagen.R", package = "spatbiome"))
```

An example data generating function is provided below. The function
`generate_LGCP` generates two Log Gaussian Cox Processes (LGCPs) in a
specified window. The user can choose to generate independent point
patterns or dependent point patterns where points of type 2 are placed
close to points of type 1.

``` r
generate_LGCP <- function(win, corrfun = "exponential",mu = 5, scale = 0.1, dependant = FALSE){

  data1 <- spatstat.random::rLGCP(corrfun, mu = mu, win = win, scale = scale, var = 1)

  if(dependant == FALSE){
    data2 <- spatstat.random::rLGCP(corrfun, mu = mu, win = win, scale = scale, var = 1)
  } else {
    data2 <- data1
    runif_sh <- spatstat.random::runifdisc(data1$n, radius = 0.1)
    data2$x <- data1$x + runif_sh$x
    data2$y <- data1$y + runif_sh$y

    keep <- spatstat.geom::inside.owin(data2$x, data2$y, win)
    data2 <- data2[keep]
  }

  spatstat.geom::marks(data1) <- factor(rep(1, spatstat.geom::npoints(data1)), levels = c(1, 2))
  spatstat.geom::marks(data2) <- factor(rep(2, spatstat.geom::npoints(data2)), levels = c(1, 2))

  data_all <- spatstat.geom::superimpose(data1, data2, W = win)
  spatstat.geom::marks(data_all) <- factor(spatstat.geom::marks(data_all), levels = c(1, 2))

  return(data_all)
}
```

## Example 1. Null hypothesis: independent point patterns

This is a basic example which demonstrates the implementation of ARShift
for two independent Log Gaussian point processes. The null hypothesis is
that the two point patterns are independent.

``` r
library(ggplot2)
X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = FALSE)
plot_X <- data.frame(X)
ggplot(plot_X, aes(x = x, y = y, color = marks)) + 
  geom_point() +
  labs(title = "Simulated independent LGCP", x = "X coordinate", y = "Y coordinate", color = "Point Type") +
  theme_minimal()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="" width="70%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.764
#> 
#> $K
#>  [1] 0.000000e+00 2.829335e-05 8.488004e-05 2.829335e-04 4.977460e-04
#>  [6] 7.240927e-04 1.054846e-03 1.352364e-03 1.833351e-03 2.288022e-03
#> [11] 2.719890e-03 3.333585e-03 3.794787e-03 4.497280e-03 5.190479e-03
#> [16] 5.892194e-03 6.767064e-03 7.779780e-03 8.551561e-03 9.393780e-03
#> [21] 1.047453e-02 1.193514e-02 1.318313e-02 1.439533e-02 1.557807e-02
#> [26] 1.653092e-02 1.795903e-02 1.904657e-02 2.053686e-02 2.214044e-02
#> [31] 2.391906e-02 2.548134e-02 2.735324e-02 2.903961e-02 3.091004e-02
#> [36] 3.282906e-02 3.465309e-02 3.649202e-02 3.869864e-02 4.093342e-02
#> [41] 4.299097e-02 4.483886e-02 4.692194e-02 4.910274e-02 5.150021e-02
#> [46] 5.375536e-02 5.554869e-02 5.767832e-02 6.035106e-02 6.312486e-02
#> 
#> $r
#>  [1] 0.000000000 0.003061224 0.006122449 0.009183673 0.012244898 0.015306122
#>  [7] 0.018367347 0.021428571 0.024489796 0.027551020 0.030612245 0.033673469
#> [13] 0.036734694 0.039795918 0.042857143 0.045918367 0.048979592 0.052040816
#> [19] 0.055102041 0.058163265 0.061224490 0.064285714 0.067346939 0.070408163
#> [25] 0.073469388 0.076530612 0.079591837 0.082653061 0.085714286 0.088775510
#> [31] 0.091836735 0.094897959 0.097959184 0.101020408 0.104081633 0.107142857
#> [37] 0.110204082 0.113265306 0.116326531 0.119387755 0.122448980 0.125510204
#> [43] 0.128571429 0.131632653 0.134693878 0.137755102 0.140816327 0.143877551
#> [49] 0.146938776 0.150000000
#> 
#> $n_attempts
#> [1] 549
#> 
#> $n_valid
#> [1] 549
#> 
#> $n_perm
#> [1] 499
```

## Example 2. Alternative hypothesis: co-clustering point patterns

Here, points of type 2 are placed close to points of type 1 to mimic
association between point types.

``` r
X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = TRUE)
plot_X <- data.frame(X)
ggplot(plot_X, aes(x = x, y = y, color = marks)) + 
  geom_point() +
  labs(title = "Simulated dependent LGCP", x = "X coordinate", y = "Y coordinate", color = "Point Type") +
  theme_minimal()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" alt="" width="70%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.008
#> 
#> $K
#>  [1] 0.000000e+00 2.097755e-05 1.887980e-04 3.146633e-04 5.663940e-04
#>  [6] 1.032365e-03 1.649023e-03 2.425192e-03 3.080502e-03 4.074369e-03
#> [11] 5.343922e-03 6.546276e-03 7.363707e-03 8.369406e-03 9.601120e-03
#> [16] 1.104050e-02 1.245895e-02 1.437307e-02 1.658198e-02 1.852601e-02
#> [21] 2.068220e-02 2.318496e-02 2.582038e-02 2.828967e-02 3.094387e-02
#> [26] 3.340607e-02 3.577716e-02 3.828891e-02 4.097915e-02 4.365514e-02
#> [31] 4.683477e-02 4.974718e-02 5.284573e-02 5.595556e-02 5.852402e-02
#> [36] 6.100328e-02 6.392414e-02 6.701225e-02 6.966800e-02 7.286243e-02
#> [41] 7.613812e-02 7.907558e-02 8.232717e-02 8.511120e-02 8.834164e-02
#> [46] 9.173156e-02 9.505769e-02 9.791920e-02 1.009262e-01 1.043200e-01
#> 
#> $r
#>  [1] 0.000000000 0.003061224 0.006122449 0.009183673 0.012244898 0.015306122
#>  [7] 0.018367347 0.021428571 0.024489796 0.027551020 0.030612245 0.033673469
#> [13] 0.036734694 0.039795918 0.042857143 0.045918367 0.048979592 0.052040816
#> [19] 0.055102041 0.058163265 0.061224490 0.064285714 0.067346939 0.070408163
#> [25] 0.073469388 0.076530612 0.079591837 0.082653061 0.085714286 0.088775510
#> [31] 0.091836735 0.094897959 0.097959184 0.101020408 0.104081633 0.107142857
#> [37] 0.110204082 0.113265306 0.116326531 0.119387755 0.122448980 0.125510204
#> [43] 0.128571429 0.131632653 0.134693878 0.137755102 0.140816327 0.143877551
#> [49] 0.146938776 0.150000000
#> 
#> $n_attempts
#> [1] 549
#> 
#> $n_valid
#> [1] 549
#> 
#> $n_perm
#> [1] 499
```

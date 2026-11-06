
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
#> ✔ 1 pkg + 40 deps: kept 40 [4.5s]
```

## Generate data

For more data generating functions, please refer to the `datagen.R`
script in the `scripts` folder.

``` r
library(spatbiome)
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

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.734
#> 
#> $K
#>  [1] 0.0000000000 0.0000000000 0.0001598465 0.0002397698 0.0003196931
#>  [6] 0.0005861040 0.0009057971 0.0012431177 0.0018178701 0.0022263173
#> [11] 0.0028287516 0.0033405012 0.0037607092 0.0044685066 0.0051111127
#> [16] 0.0058117311 0.0067692398 0.0076091569 0.0085567557 0.0098494532
#> [21] 0.0107741500 0.0119637650 0.0133247823 0.0148594433 0.0163486824
#> [26] 0.0177488920 0.0190152092 0.0204699008 0.0219333402 0.0238452746
#> [31] 0.0256291224 0.0278174352 0.0299507035 0.0319338691 0.0338163745
#> [36] 0.0355390238 0.0376603552 0.0395596056 0.0418380610 0.0449112985
#> [41] 0.0472506869 0.0500893856 0.0524900363 0.0550074434 0.0578975107
#> [46] 0.0614293714 0.0642734868 0.0672939882 0.0701635457 0.0732117117
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

## Example 2. Null hypothesis: dependent point patterns

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

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.008
#> 
#> $K
#>  [1] 0.000000e+00 7.410701e-05 2.964280e-04 3.705351e-04 5.558026e-04
#>  [6] 8.522306e-04 1.259819e-03 1.963836e-03 2.651382e-03 3.402544e-03
#> [11] 4.342162e-03 5.564305e-03 6.753510e-03 8.169117e-03 9.314459e-03
#> [16] 1.040607e-02 1.187087e-02 1.332239e-02 1.473571e-02 1.569910e-02
#> [21] 1.761700e-02 1.946405e-02 2.117340e-02 2.264205e-02 2.442902e-02
#> [26] 2.666466e-02 2.882016e-02 3.055601e-02 3.269732e-02 3.498960e-02
#> [31] 3.749778e-02 3.993395e-02 4.265239e-02 4.499247e-02 4.724654e-02
#> [36] 4.922195e-02 5.144335e-02 5.376509e-02 5.688665e-02 5.951856e-02
#> [41] 6.239257e-02 6.501051e-02 6.788052e-02 7.026818e-02 7.312592e-02
#> [46] 7.728439e-02 8.015372e-02 8.350443e-02 8.667021e-02 8.987049e-02
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

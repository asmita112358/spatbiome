
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatbiome

A collection of functions for spatial analysis of microbiome data.
Currently contains ARShift, a random shift test for bivariate point
patterns.

## Installation

You can install the development version of spatbiome from
[GitHub](https://github.com/) with:

``` r
pak::pak("asmita112358/spatbiome", dependencies = TRUE)
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
Change the number of workers according to the number of cores available
in your machine.

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.808
#> 
#> $K
#>  [1] 0.000000e+00 0.000000e+00 4.075976e-05 1.902122e-04 3.532513e-04
#>  [6] 6.385696e-04 1.018994e-03 1.420570e-03 1.787408e-03 2.344458e-03
#> [11] 3.064548e-03 3.675944e-03 4.314514e-03 5.117531e-03 5.692513e-03
#> [16] 6.608430e-03 7.524357e-03 8.529535e-03 9.668238e-03 1.070242e-02
#> [21] 1.190695e-02 1.303161e-02 1.451770e-02 1.587737e-02 1.709487e-02
#> [26] 1.845848e-02 1.973588e-02 2.121829e-02 2.281568e-02 2.440832e-02
#> [31] 2.610007e-02 2.769921e-02 2.954450e-02 3.147971e-02 3.331678e-02
#> [36] 3.514919e-02 3.688439e-02 3.900727e-02 4.119179e-02 4.320210e-02
#> [41] 4.577714e-02 4.806735e-02 5.046895e-02 5.265650e-02 5.451654e-02
#> [46] 5.711060e-02 5.982784e-02 6.221078e-02 6.474052e-02 6.766795e-02
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
#> [1] 0.002
#> 
#> $K
#>  [1] 0.000000e+00 8.865248e-05 3.071488e-04 5.878816e-04 8.733823e-04
#>  [6] 1.301869e-03 2.001825e-03 2.724094e-03 3.284242e-03 4.064773e-03
#> [11] 4.974974e-03 5.791050e-03 7.121126e-03 8.152503e-03 9.385716e-03
#> [16] 1.054538e-02 1.208790e-02 1.355098e-02 1.534196e-02 1.702904e-02
#> [21] 1.889649e-02 2.119165e-02 2.309122e-02 2.508585e-02 2.749688e-02
#> [26] 2.985181e-02 3.239926e-02 3.453657e-02 3.678365e-02 3.950171e-02
#> [31] 4.196018e-02 4.485399e-02 4.765082e-02 4.985079e-02 5.270706e-02
#> [36] 5.508196e-02 5.808329e-02 6.109218e-02 6.345930e-02 6.602463e-02
#> [41] 6.920328e-02 7.166897e-02 7.474927e-02 7.796014e-02 8.059412e-02
#> [46] 8.310472e-02 8.566928e-02 8.854268e-02 9.096850e-02 9.380421e-02
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

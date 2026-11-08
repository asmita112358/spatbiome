
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
#> ✔ 1 pkg + 40 deps: kept 40 [4.3s]
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

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.518
#> 
#> $K
#>  [1] 0.000000e+00 3.702607e-05 1.295912e-04 2.406694e-04 3.948234e-04
#>  [6] 7.280580e-04 9.872405e-04 1.283449e-03 1.686349e-03 2.014552e-03
#> [11] 2.580149e-03 3.171922e-03 3.846875e-03 4.525279e-03 5.435653e-03
#> [16] 6.266032e-03 7.386533e-03 8.271505e-03 9.435031e-03 1.041306e-02
#> [21] 1.165767e-02 1.297108e-02 1.443229e-02 1.628515e-02 1.758996e-02
#> [26] 1.905852e-02 2.018625e-02 2.154306e-02 2.306018e-02 2.487371e-02
#> [31] 2.664522e-02 2.838744e-02 3.004947e-02 3.207292e-02 3.429941e-02
#> [36] 3.613092e-02 3.794093e-02 4.025927e-02 4.247522e-02 4.464237e-02
#> [41] 4.681722e-02 4.883006e-02 5.094687e-02 5.339222e-02 5.589468e-02
#> [46] 5.851065e-02 6.136705e-02 6.432498e-02 6.668385e-02 6.963216e-02
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

<img src="man/figures/README-unnamed-chunk-7-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.002
#> 
#> $K
#>  [1] 0.000000e+00 6.426322e-05 2.731187e-04 5.944348e-04 1.060343e-03
#>  [6] 1.654778e-03 2.178408e-03 2.887907e-03 3.806873e-03 4.761770e-03
#> [11] 5.673299e-03 6.792481e-03 7.901873e-03 9.467146e-03 1.117503e-02
#> [16] 1.278648e-02 1.468525e-02 1.643385e-02 1.834201e-02 2.038201e-02
#> [21] 2.248643e-02 2.449750e-02 2.665492e-02 2.887610e-02 3.156259e-02
#> [26] 3.427074e-02 3.691579e-02 3.986143e-02 4.290812e-02 4.556909e-02
#> [31] 4.840722e-02 5.152933e-02 5.479912e-02 5.797748e-02 6.090388e-02
#> [36] 6.405428e-02 6.740900e-02 7.049741e-02 7.381396e-02 7.680292e-02
#> [41] 7.999821e-02 8.342107e-02 8.716612e-02 9.049001e-02 9.401770e-02
#> [46] 9.756456e-02 1.013892e-01 1.051079e-01 1.085411e-01 1.124066e-01
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

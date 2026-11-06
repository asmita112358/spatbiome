
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
#> ✔ Updated metadata database: 7.33 MB in 14 files.
#> ℹ Updating metadata database✔ Updating metadata database ... done
#>  
#> → Package library at '/Library/Frameworks/R.framework/Versions/4.6/Resources/library'.
#> ℹ No downloads are needed
#> ✔ 1 pkg + 40 deps: kept 40 [9.8s]
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
    runif_sh <- spatstat.random::runifdisc(data1$n, radius = 0.2)
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
X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = FALSE)
plot(X, main = "Simulated independent LGCP")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.978
#> 
#> $K
#>  [1] 0.000000e+00 0.000000e+00 4.638219e-05 1.391466e-04 3.555968e-04
#>  [6] 6.388793e-04 8.576736e-04 1.107047e-03 1.693899e-03 2.107979e-03
#> [11] 2.477456e-03 2.902251e-03 3.619701e-03 4.614773e-03 5.423976e-03
#> [16] 6.197026e-03 7.171567e-03 8.081193e-03 8.800024e-03 9.950374e-03
#> [21] 1.094189e-02 1.208581e-02 1.328100e-02 1.446550e-02 1.622188e-02
#> [26] 1.743779e-02 1.906610e-02 2.040570e-02 2.194648e-02 2.384951e-02
#> [31] 2.574690e-02 2.765557e-02 2.969256e-02 3.150141e-02 3.313605e-02
#> [36] 3.513753e-02 3.738270e-02 3.940355e-02 4.144544e-02 4.332919e-02
#> [41] 4.520679e-02 4.774279e-02 4.974814e-02 5.257372e-02 5.472340e-02
#> [46] 5.754131e-02 6.007532e-02 6.245702e-02 6.507961e-02 6.764358e-02
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
plot(X, main = "Simulated dependent LGCP")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="" width="100%" />

``` r
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")
res
#> $p
#> [1] 0.002
#> 
#> $K
#>  [1] 0.0000000000 0.0001747030 0.0003105831 0.0006211662 0.0009511608
#>  [6] 0.0014752698 0.0020724307 0.0026935970 0.0033924090 0.0043004212
#> [11] 0.0053680507 0.0066880290 0.0080801185 0.0093029635 0.0104311789
#> [16] 0.0119153516 0.0131553170 0.0151246590 0.0166368232 0.0183234268
#> [21] 0.0200989223 0.0222146912 0.0244851611 0.0265865302 0.0287634592
#> [26] 0.0310283010 0.0338795159 0.0367619803 0.0393546887 0.0417828531
#> [31] 0.0447244058 0.0479815226 0.0507638352 0.0535794006 0.0571314807
#> [36] 0.0603820265 0.0635348447 0.0672794404 0.0705643397 0.0741141128
#> [41] 0.0775584261 0.0813193541 0.0852693537 0.0887232375 0.0926399377
#> [46] 0.0965988340 0.1009346834 0.1051919224 0.1094183366 0.1134878228
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


<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatbiome

A collection of functions for spatial analysis of microbiome data.
Currently contains ARShift, a random shift test for bivariate point
patterns.

## Installation

You can install the development version of spatbiome from
[GitHub](https://github.com/) with:

``` r
pak::pak("asmita112358/spatbiome")
#> ℹ Loading metadata database✔ Loading metadata database ... done
#>  
#> → Package library at '/Library/Frameworks/R.framework/Versions/4.6/Resources/library'.
#> → Will update 1 package.
#> → Will download 1 package with unknown size.
#> + spatbiome 0.1.0 → 0.1.0 👷🏽‍♂️🔧 ⬇ (GitHub: aefd986)
#> ℹ Getting 1 pkg with unknown size
#> ✔ Cached copy of spatbiome 0.1.0 (source) is the latest build
#> ✔ Installed spatbiome 0.1.0 (github::asmita112358/spatbiome@aefd986) (19ms)
#> ✔ 1 pkg + 40 deps: kept 39, upd 1 [4.2s]
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
#> [1] 0.086
#> 
#> $K
#>  [1] 0.0000000000 0.0001095770 0.0001095770 0.0003686007 0.0005000932
#>  [6] 0.0007411626 0.0009384013 0.0012890478 0.0016148341 0.0020905580
#> [11] 0.0028408995 0.0033643323 0.0041413840 0.0048351346 0.0057168404
#> [16] 0.0065494905 0.0075813417 0.0085437241 0.0096419515 0.0107469930
#> [21] 0.0118015558 0.0128513016 0.0140879084 0.0152223303 0.0166095903
#> [26] 0.0182126128 0.0195767597 0.0209743239 0.0225258143 0.0243431732
#> [31] 0.0261985763 0.0276961127 0.0297173322 0.0314113176 0.0332157374
#> [36] 0.0354443485 0.0375964775 0.0396608913 0.0424924407 0.0445122195
#> [41] 0.0465131321 0.0487287100 0.0510605550 0.0533563527 0.0557592833
#> [46] 0.0582804378 0.0612650564 0.0639103498 0.0663434553 0.0689686948
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
#>  [1] 0.0000000000 0.0000621118 0.0001656315 0.0004347826 0.0007742486
#>  [6] 0.0014763089 0.0018903876 0.0025943214 0.0036452689 0.0045163171
#> [11] 0.0055804955 0.0069288131 0.0082921706 0.0099247913 0.0114246953
#> [16] 0.0132108965 0.0146735830 0.0165915883 0.0185039934 0.0205321945
#> [21] 0.0230823406 0.0250996366 0.0274502526 0.0299356302 0.0323112293
#> [26] 0.0347864281 0.0372355356 0.0402265483 0.0430968283 0.0462454903
#> [31] 0.0493210775 0.0519871211 0.0551998449 0.0584118817 0.0613018047
#> [36] 0.0640923354 0.0667637549 0.0698772523 0.0728391855 0.0756903210
#> [41] 0.0790827907 0.0824124876 0.0850701617 0.0894418388 0.0930944468
#> [46] 0.0969539691 0.1007869496 0.1044797229 0.1089896882 0.1136519465
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

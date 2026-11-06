#load functions in datagen.R here
source("scripts/datagen.R", echo = TRUE)

## Simulate p-values under the null hypothesis
p <- c()
for(sim in 1:500){
  X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = FALSE)
  res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 10, type = "hom")
  p[sim] <- res$p
}

mean(p <= 0.05)
qqplot(ppoints(500), p, main = "QQ plot of p-values", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")


##Simulate p-values under the alternative hypothesis (co-clustering)
p <- c()
for(sim in 1:500){
  X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = TRUE)
  res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 10, type = "hom")
  p[sim] <- res$p
}

mean(p <= 0.05)

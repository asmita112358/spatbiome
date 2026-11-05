X <- generate_LGCP(win = spatstat.geom::owin(c(0, 1), c(0, 1)), corrfun = "exponential", mu = 5, scale = 0.1, dependant = TRUE)
res <- ARShift(data = X, i = "1", j = "2", n_perm = 499, workers = 5, type = "hom")

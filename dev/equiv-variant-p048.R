# p048 control: the univariate t2() model on the same data, read at L4
spec$data$y <- spec$data$y1
spec$bf <- "bf(y ~ t2(x, z))"
spec$l4 <- TRUE
spec$no_points <- TRUE

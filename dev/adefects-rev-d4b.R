source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(brms) })

cat("== update(newdata =) must carry the NEW frame\n")
set.seed(20260917)
d1 <- data.frame(x = rnorm(40)); d1$y <- rnorm(40, 1 + d1$x, 1)
d2 <- data.frame(x = rnorm(17)); d2$y <- rnorm(17, 1 + d2$x, 1)
f1 <- frm(frmtmb::bf(y ~ x), family = gaussian(), data = d1)
f2 <- update(f1, newdata = d2)
cat("  nrow(f1$data) =", nrow(f1$data), " nrow(f2$data) =", nrow(f2$data),
    " nrow(model.frame(f2)) =", nrow(model.frame(f2)), "\n")
cat("  f2$data is the new frame:",
    isTRUE(all.equal(f2$data$y, d2$y)), "\n")
cat("  attr(f2$data, 'data_name') =",
    format(attr(f2$data, "data_name")), "\n")

cat("\n== NA rows: fit$data is the model frame, not the input data\n")
d3 <- d1; d3$x[1:3] <- NA
f3 <- suppressWarnings(frm(frmtmb::bf(y ~ x), family = gaussian(),
                           data = d3))
cat("  nrow(input) =", nrow(d3), " nrow(fit$data) =", nrow(f3$data), "\n")

cat("\n== column names: frmtmb's frame against brms's data element\n")
bf1 <- get("brmsfit_example1", envir = asNamespace("brms"))
cat("  brms  fit1$data:", paste(names(bf1$data), collapse = ", "), "\n")
cat("  brms  dim:", paste(dim(bf1$data), collapse = "x"), "\n")
cat("  brms  attrs:", paste(setdiff(names(attributes(bf1$data)),
                                    "row.names"), collapse = ", "), "\n")
dd <- data.frame(count = rpois(40, 3), Trt = factor(rep(0:1, 20)),
                 Age = rnorm(40), volume = runif(40))
fo <- frm(frmtmb::bf(count ~ Trt + Age + offset(Age)),
          family = poisson(), data = dd)
cat("  frmtmb fit$data:", paste(names(fo$data), collapse = ", "), "\n")
cat("  frmtmb dim:", paste(dim(fo$data), collapse = "x"), "\n")
cat("  frmtmb attrs:", paste(setdiff(names(attributes(fo$data)),
                                     "row.names"), collapse = ", "), "\n")

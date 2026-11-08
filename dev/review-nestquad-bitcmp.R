# Compare the bitwise RDS of the two arms.
o <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-out"
a <- readRDS(file.path(o, "bitwise-lane-1.rds")); b <- readRDS(file.path(o, "bitwise-base-1.rds"))
for (nm in names(a)) {
  x <- a[[nm]]; y <- b[[nm]]
  same <- vapply(c("par", "obj", "ll", "vc", "b"), function(k) identical(x[[k]], y[[k]]), TRUE)
  md <- if (!is.null(x$par)) max(abs(x$par - y$par)) else NA
  cat(sprintf("BITWISE %-10s identical %s  max|dpar| %.3g %s\n", nm,
              paste(names(same), same, sep = "=", collapse = " "), md,
              paste(c(x$err, y$err, ""), collapse = " ")))
}

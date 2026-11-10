# test-bcm-mpt.R's hierarchical MPT fit: what each build says.
#   CODE8_ARM=before|after Rscript dev/code8-mpt.R
source("dev/code8-common.R")
options(width = 200)
env <- new.env(parent = asNamespace("frmtmb"))
invisible(utils::capture.output(
  testthat::source_test_helpers("tests/testthat", env = env)))
# the test file's own definitions (its top-level assignments)
for (ex in parse("tests/testthat/test-bcm-mpt.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bcm_mpt_data()
r <- rfit(env$bcm_mpt_hier_formula(), family = env$bcm_mpt_pairs(), data = d)
f <- r$f
cat("objective", format(f$opt$objective, digits = 12), "code",
    f$opt$convergence, f$opt$message, "\n")
cat(gsub(" || ", "\n", cond_txt(r), fixed = TRUE), "\n")
g <- drop(f$obj$gr(f$opt$par))
cat("max|grad|", format(max(abs(g))), "\n")
str(fns$grad_verdict(f, g, f$control))
room <- fns$conv_headroom(f, g)
cat("conv_headroom", format(room), "\n")
invisible(utils::capture.output(dg <- diagnose(f)))
cat("bad_se:", dg$bad_se, "\n")
sdr <- quiet(fns$sdr_of(f))
cat("se_lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
print(data.frame(name = fns$outer_par_names(f), est = f$opt$par, g = g,
                 se = sqrt(diag(sdr$cov.fixed))), digits = 4)
# the probe's own evaluations, and the objective along -g
p <- f$opt$par
H <- fns$fit_outer_hessian(f)$H
cat("diag H:", format(diag(H), digits = 3), "\n")
f0 <- as.numeric(f$obj$fn(p))
gn <- -g / sqrt(sum(g^2))
for (t in 10^(-12:-3)) {
  cat(sprintf("t %.0e  f(p + t dir) - f(p) = %.6g\n", t,
              as.numeric(f$obj$fn(p + t * gn)) - f0))
}
# the objective along each coordinate: a kink grows linearly with |t|
nm <- fns$outer_par_names(f)
for (j in seq_along(p)) {
  out <- sapply(10^(-9:-5), function(t) {
    e <- numeric(length(p)); e[j] <- t
    c(as.numeric(f$obj$fn(p + e)) - f0, as.numeric(f$obj$fn(p - e)) - f0)
  })
  cat(sprintf("%-14s +: %s\n%-14s -: %s\n", nm[j],
              paste(format(out[1, ], digits = 3), collapse = " "), "",
              paste(format(out[2, ], digits = 3), collapse = " ")))
}

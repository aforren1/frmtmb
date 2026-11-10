# Review of lane code8: the poly(x, 2) tight-sd fit loses 8 of 10 SEs on
# the lane. Reference: the Richardson Hessian at the stop, its unit-
# diagonal spectrum, and the SEs of the fixed effects with the theta
# block held (the Schur complement is not needed for that question:
# the fixed effects' SEs from the full inverse and from the pseudo-
# inverse past the flat directions).
#   RV_ARM=after Rscript dev/review-code8-tightsd2.R
source("dev/review-code8-common.R")
options(width = 180)
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22", group = "g")
f <- quiet(frm(form, data = d, prior = p))
x <- f$opt$par
nm <- ns$outer_par_names(f)
obj <- f$obj
J <- numDeriv::jacobian(function(q) drop(obj$gr(q)), x, method = "Richardson")
H <- (J + t(J)) / 2
cat("asym", max(abs(J - t(J))), "\n")
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("unit-diag eigenvalues:", signif(e$values, 3), "\n")
for (k in which(e$values < 1e-3 * max(e$values))) {
  v <- e$vectors[, k]; o <- order(-abs(v))[1:4]
  cat(" ev", signif(e$values[k], 3), ":",
      paste(sprintf("%s %.3f", nm[o], v[o]), collapse = ", "), "\n")
  dn <- v / D
  f0 <- obj$fn(x)
  for (t in c(0.1, 1, 3)) cat(sprintf("    t %g: %.3e %.3e\n", t,
                                    obj$fn(x + t * dn) - f0,
                                    obj$fn(x - t * dn) - f0))
}
ok <- e$values > 1e-9 * max(e$values)
V <- (e$vectors[, ok] %*% (t(e$vectors[, ok]) / e$values[ok])) / outer(D, D)
s <- se_of(f)
print(data.frame(par = nm, est = signif(x, 4), lane_se = signif(s$se, 4),
                 lane_lost = s$lost[nm] %||% NA,
                 ref_se = signif(sqrt(diag(V)), 4)))
cat("DONE\n")

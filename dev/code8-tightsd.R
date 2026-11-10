# Review M1: the poly(x, 2) tight-sd fit of test-brms-names-specials.R.
# Where the finite-difference Hessian's noise sits, what the RE block's
# correlations are, and what se_edge_cor() and se_at_edge() say.
#   CODE8_ARM=before|after Rscript dev/code8-tightsd.R
source("dev/code8-common.R")
options(width = 200)
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22",
               group = "g")
f <- quiet(frm(form, data = d, prior = p))
nm <- fns$outer_par_names(f)
pp <- f$opt$par
print(setNames(signif(pp, 5), nm))
print(quiet(VarCorr(f))$g)
h <- fns$fit_outer_hessian(f)
cat("diag H:", signif(diag(h$H), 3), "\n")
cat("row max of E:", signif(apply(h$E, 1, max), 3), "\n")
print(signif(h$E, 2))
for (j in grep("theta", nm)) {
  cat(nm[j], "edge -:", fns$se_at_edge(f, nm[j], -1), " edge +:",
      fns$se_at_edge(f, nm[j], 1), "\n")
}
# the objective along each theta, others held, to see where it is flat
f0 <- as.numeric(f$obj$fn(pp))
for (j in grep("theta", nm)) {
  ch <- sapply(c(-10, -1, -0.1, 0.1, 1, 10), function(t) {
    q <- pp; q[j] <- q[j] + t; as.numeric(f$obj$fn(q)) - f0
  })
  cat(sprintf("%-8s %s\n", nm[j], paste(format(ch, digits = 3),
                                         collapse = " ")))
}
# tier 3's eigen step over the rows it keeps (theta_5, theta_6 emptied)
free <- 1:8
H <- h$H; E <- h$E
Hf <- H[free, free]; Df <- sqrt(abs(diag(Hf)))
e <- eigen(Hf / outer(Df, Df), symmetric = TRUE)
Es <- E[free, free] / outer(Df, Df)
nk <- sqrt(colSums((Es %*% e$vectors)^2))
cat("ev:", signif(e$values, 3), "\nnk:", signif(nk, 3), "\nthr:",
    signif(pmax(1e-9 * max(abs(e$values)), 3 * nk), 3), "\n")
print(round(e$vectors, 3))
an <- quiet(fns$se_fit_analysis(f, force = TRUE))
cat("analysis lost:", paste(names(an$lost), an$lost), "\n")

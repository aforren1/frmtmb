# Review M1 prototype: the 3-by-3 block's third Cholesky row (theta_5,
# theta_6 at 808 and 1814) scaled back along itself; the objective along
# that scaling, and the Hessian's noise and tier 3's verdict at the
# pulled-back point.
#   CODE8_ARM=after Rscript dev/code8-tightsd-pull.R
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
row <- match(c("theta_5", "theta_6"), nm)
f0 <- as.numeric(f$obj$fn(pp))
tc <- fns$se_cor_clamp_theta()
cat("clamp theta", tc, "\n")
for (cc in c(tc / max(abs(pp[row])), 0.1, 0.5, 2, 10, 100)) {
  q <- pp; q[row] <- q[row] * cc
  cat(sprintf("row scale %.4g: f - f0 = %.4g\n", cc,
              as.numeric(f$obj$fn(q)) - f0))
}
q <- pp; q[row] <- q[row] * tc / max(abs(pp[row]))
h <- fns$fd_hessian(q, function(z) drop(f$obj$gr(z)),
                    eps = fns$outer_fd_steps(f$frame, names(pp)))
cat("asym max at the pulled point:", max(h$E), " row max of E:",
    signif(apply(h$E, 1, max), 3), "\n")
keep <- setdiff(seq_along(pp), row)
Hk <- h$H[keep, keep]
V <- solve(Hk)
cat("SEs of the others at the pulled point:\n")
print(setNames(signif(sqrt(diag(V)), 4), nm[keep]))
an <- fns$se_tier3(f, h$H, h$E, q, exact = FALSE,
                   held = logical(length(pp)),
                   edge = seq_along(pp) %in% row)
cat("tier 3 at the pulled point, row held:", paste(names(an$lost), an$lost),
    "\n")

# Review of lane code8: the tight sd prior of test-brms-names-specials.R,
# whose gradient warning the lane's probe removed. Is the stop the
# maximum? Polish it with nlminb and optim(BFGS) from the stop on the
# fit's own objective and report the gain, with the warnings each build
# gives.
#   RV_ARM=before|after Rscript dev/review-code8-tightsd.R
source("dev/review-code8-common.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
for (cs in list(c("y ~ I(x^2) + (1 + I(x^2) | g)", "I(x^2)", "IxE2"),
                c("y ~ poly(x, 2) + (1 + poly(x, 2) | g)", "poly(x, 2)2",
                  "polyx22"))) {
  form <- bf(stats::as.formula(cs[1]))
  for (spell in cs[2:3]) {
    p <- set_prior("normal(2, 0.001)", class = "sd", coef = spell, group = "g")
    r <- rfit(form, data = d, prior = p)
    f <- r$f
    obj <- f$obj
    x <- f$opt$par
    f0 <- obj$fn(x)
    o1 <- nlminb(x, obj$fn, obj$gr, control = list(iter.max = 2000,
                                                    eval.max = 4000))
    o2 <- optim(x, obj$fn, obj$gr, method = "BFGS",
                control = list(maxit = 5000, reltol = 1e-14))
    g <- drop(obj$gr(x))
    cat(sprintf("%s | %s: code %d gmax %.3g; polish gain nlminb %.3g, BFGS %.3g\n",
                cs[1], spell, f$opt$convergence, max(abs(g)),
                f0 - o1$objective, f0 - o2$value))
    cat("   conds:", cond_txt(r), "\n")
    s <- se_of(f)
    cat("   lost:", paste(names(s$lost), s$lost), "\n")
  }
}
cat("DONE\n")

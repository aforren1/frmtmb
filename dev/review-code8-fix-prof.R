# Fix-round check (M1): the poly tight-sd fit on the fix build. Its SEs
# and lost set, and an independent profile of theta_1, theta_2, theta_4
# at +-0.5 and +-1 of the fix build's SEs (rest re-optimized from the
# stop, nlminb and then BFGS, best kept).
#   RV_ARM=after Rscript dev/review-code8-fix-prof.R
source("dev/review-code8-common.R")
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22", group = "g")
r <- rfit(form, data = d, prior = p)
f <- r$f
cat("conds:", cond_txt(r), "\n")
s <- se_of(f)
print(signif(s$se, 4)); print(s$lost)
x <- f$opt$par; nm <- names(s$se); obj <- f$obj; f0 <- obj$fn(x)
for (pn in c("theta_1", "theta_2", "theta_4", "(Intercept)")) {
  i <- match(pn, nm)
  sc <- s$se[[pn]]
  out <- sapply(c(-1, -0.5, 0.5, 1), function(m) {
    st <- x; st[i] <- x[i] + m * sc
    oth <- setdiff(seq_along(x), i)
    fo <- function(q) { z <- st; z[oth] <- q; obj$fn(z) }
    go <- function(q) { z <- st; z[oth] <- q; drop(obj$gr(z))[oth] }
    o1 <- nlminb(x[oth], fo, go, control = list(iter.max = 2000,
                                                 eval.max = 4000))
    o2 <- optim(o1$par, fo, go, method = "BFGS",
                control = list(maxit = 5000, reltol = 1e-14))
    min(o1$objective, o2$value) - f0
  })
  cat(sprintf("%-12s SE %.4g  loss at -1,-0.5,0.5,1 SE: %s\n", pn, sc,
              paste(signif(out, 3), collapse = " ")))
}
cat("DONE\n")

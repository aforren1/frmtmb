# Attack 4, profile: why confint(method = "profile") on a nested fit
# excludes its own estimate. Same data as review-nestquad-postfit.R.
#   Rscript dev/review-nestquad-profile.R
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
set.seed(20261107)
ns_ <- 14
d <- do.call(rbind, lapply(seq_len(ns_), function(s) {
  do.call(rbind, lapply(seq_len(sample(1:4, 1)), function(c) {
    data.frame(school = s, cls = c, x = rnorm(sample(2:5, 1)))
  }))
}))
d$school <- factor(d$school)
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
ua <- rnorm(ns_, 0, 0.9); ui <- rnorm(nlevels(d$cls_uni), 0, 1)
d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.7 * d$x + ua[d$school] +
                                   ui[d$cls_uni]))
fo <- bf(y ~ x + (1 | school/cls_uni)) + bernoulli()
fq <- frm(fo, data = d, quadrature = TRUE)
fl <- frm(fo, data = d)
par <- fq$opt$par
cat("par:", format(par, digits = 6), "\n")
cat("last.par.best:", format(fq$obj$env$last.par.best, digits = 6), "\n")
cat("obj$par:", format(fq$obj$par, digits = 6), "\n")
pr <- TMB::tmbprofile(fq$obj, name = 2, trace = FALSE)
print(head(pr[order(pr[[1]]), ], 60))
cat("profile min value", min(pr$value), "at", pr[[1]][which.min(pr$value)],
    " fit objective", fq$opt$objective, "\n")
print(confint(pr))
# the frozen tape against the rule re-anchored at each point, others fixed
for (dx in c(-0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6)) {
  q <- par; q[2] <- q[2] + dx
  cat(sprintf("x %+.2f frozen %.6f reanchored %.6f laplace %.6f\n", dx,
              -fq$obj$fn(q), nest_at(fq, q, c(11, 15)), -fl$obj$fn(q)))
}
# Laplace and one-block quadrature fits, for the convention
cat("laplace profile CI:\n"); print(confint(fl, method = "profile", parm = "x"))
cat("laplace wald CI:\n"); print(confint(fl)[2, ])
f1 <- frm(bf(y ~ x + (1 | school)) + bernoulli(), data = d, quadrature = TRUE)
cat("one-block quad profile CI:\n")
print(try(confint(f1, method = "profile", parm = "x")))
cat("one-block quad wald CI:\n"); print(confint(f1)[2, ])
cat("DONE\n")

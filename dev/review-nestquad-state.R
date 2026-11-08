# Where the frozen nested tape goes below its own optimum: fix beta_x,
# minimize the frozen tape over the rest, then re-anchor there.
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
fq <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d,
          quadrature = TRUE)
par <- fq$opt$par
lp <- function(eta, r) sum(dbinom(d$y[r], 1, plogis(eta), log = TRUE))
cat("Wald fresh:", confint(fq)[2, ], "\n")
fq2 <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d, quadrature = TRUE)
cat("last.par.best before:", format(fq2$obj$env$last.par.best, digits = 6), "\n")
cu <- confint(fq2, method = "uniroot", parm = "x"); print(cu)
cat("last.par.best after uniroot:", format(fq2$obj$env$last.par.best, digits = 6), "\n")
cat("Wald after uniroot:", confint(fq2)[2, ], "\n")
cat("logLik after uniroot:", as.numeric(logLik(fq2)), "\n")
fq3 <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d, quadrature = TRUE)
cp <- confint(fq3, method = "profile", parm = "x")
cat("last.par.best after profile:", format(fq3$obj$env$last.par.best, digits = 6), "\n")
cat("Wald after profile:", confint(fq3)[2, ], "\n")

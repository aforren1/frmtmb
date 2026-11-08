# What a correct profile interval for x costs on a nested quadrature
# fit: one nested refit per profile point, with x fixed through an
# offset. The review's design (117 Bernoulli rows, seed 20261107).
#   Rscript dev/nestquad-profcost.R
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
t0 <- proc.time()[[3]]
fq <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d,
          quadrature = TRUE)
t_fit <- proc.time()[[3]] - t0
ll0 <- as.numeric(logLik(fq))
bx <- unname(fq$estimates$beta[2])
cat(sprintf("PROF fit %.2fs logLik %.8f x %.6f settled %s rounds %d\n", t_fit,
            ll0, bx, fq$quad_nested$settled, fq$quad_nested$rounds))
fl <- frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d)
t0 <- proc.time()[[3]]
pr <- TMB::tmbprofile(fl$obj, name = 2, trace = FALSE)
cat(sprintf("PROF laplace tmbprofile points %d time %.2fs\n", nrow(pr),
            proc.time()[[3]] - t0))
n_fit <- 0
t_all <- 0
drop2 <- function(c) {
  d$oc <- c * d$x
  t0 <- proc.time()[[3]]
  f <- suppressWarnings(frm(bf(y ~ 1 + offset(oc) + (1 | school/cls_uni)) +
                              bernoulli(), data = d, quadrature = TRUE))
  t_all <<- t_all + proc.time()[[3]] - t0
  n_fit <<- n_fit + 1
  2 * (ll0 - as.numeric(logLik(f)))
}
w <- confint(fq)[2, ]
lo <- uniroot(function(c) drop2(c) - 3.84, c(w[1] - 1, bx), tol = 1e-4)$root
hi <- uniroot(function(c) drop2(c) - 3.84, c(bx, w[2] + 1), tol = 1e-4)$root
cat(sprintf("PROF refit interval (%.4f, %.4f) Wald (%.4f, %.4f)\n", lo, hi,
            w[1], w[2]))
cat(sprintf("PROF refits %d total %.1fs mean %.2fs per refit\n", n_fit, t_all,
            t_all / n_fit))
cat(sprintf("PROF 2*drop at the endpoints %.4f %.4f\n", drop2(lo), drop2(hi)))

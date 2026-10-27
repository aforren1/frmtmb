# Lane loglik: the grouped student() cases of dev/loglik-brms-probe.R
# again, with the reference taken group by group. brms's time-series
# student likelihood is one multivariate t PER GROUP (its
# student_t_time_*_lpdf loops over groups), so the first probe's
# reference, one t over all rows, was a different model and its 0.1
# differences say nothing about brms. Seed 1, data as in the first
# probe. Output: dev/loglik-brms-probe2.txt.
source("C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/loglik-brms-lib.R")
set.seed(1)
sizes <- c(6, 5, 4, 6, 3)
d <- do.call(rbind, lapply(seq_along(sizes), function(g) {
  data.frame(g = g, t = seq_len(sizes[g]))
}))
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d))
d <- d[sample(nrow(d)), ]
rownames(d) <- NULL
n <- nrow(d)

grouped_cond <- function(dd, m, S, nu) {
  out <- numeric(nrow(dd))
  for (g in unique(dd$g)) {
    r <- which(dd$g == g)
    out[r] <- ll_dense_cond(dd$y[r], m[r], S[r, r, drop = FALSE], nu)
  }
  out
}
grp_cov <- function(dd, corfun) {
  S <- matrix(0, nrow(dd), nrow(dd))
  for (g in unique(dd$g)) {
    r <- which(dd$g == g)
    r <- r[order(dd$t[r])]
    S[r, r] <- corfun(length(r), dd$t[r])
  }
  S
}
ar1 <- function(a) function(k, tt) a^abs(outer(1:k, 1:k, "-")) / (1 - a^2)
ma1 <- function(m) function(k, tt) {
  M <- diag(1 + m^2, k)
  M[abs(outer(1:k, 1:k, "-")) == 1] <- m
  M
}
arma1 <- function(a, m) function(k, tt) {
  g0 <- (1 + 2 * a * m + m^2) / (1 - a^2)
  g1 <- (1 + a * m) * (a + m) / (1 - a^2)
  L <- abs(outer(1:k, 1:k, "-"))
  out <- g1 * a^pmax(L - 1, 0)
  out[L == 0] <- g0
  out
}
cosy1 <- function(r) function(k, tt) {
  M <- matrix(r, k, k)
  diag(M) <- 1
  M
}
mu_of <- function(p, dd) p$Intercept + p$b[1] * (dd$x - mean(dd$x))
report <- function(nm, bfit, ref) {
  ll <- log_lik(bfit)
  cat(sprintf("%-26s draws %d  rows %d  max|brms - dense| / max|dense| = %.3e\n",
              nm, nrow(ll), ncol(ll), max(abs(ll - ref)) / max(abs(ref))))
}
base <- function(extra) {
  list(c(list(b = array(0.4, 1), Intercept = 1.1, sigma = 0.8, nu = 5),
         extra[[1]]),
       c(list(b = array(0.6, 1), Intercept = 0.9, sigma = 1.2, nu = 12),
         extra[[2]]))
}
cases <- list(
  ar = list(f = bf(y ~ x + ar(t, g, cov = TRUE)),
            extra = list(list(ar = array(0.5, 1)), list(ar = array(-0.3, 1))),
            cov = function(p) grp_cov(d, ar1(p$ar[1]))),
  ma = list(f = bf(y ~ x + ma(t, g, cov = TRUE)),
            extra = list(list(ma = array(0.4, 1)), list(ma = array(-0.6, 1))),
            cov = function(p) grp_cov(d, ma1(p$ma[1]))),
  arma = list(f = bf(y ~ x + arma(t, g, cov = TRUE)),
              extra = list(list(ar = array(0.5, 1), ma = array(0.3, 1)),
                           list(ar = array(-0.2, 1), ma = array(0.5, 1))),
              cov = function(p) grp_cov(d, arma1(p$ar[1], p$ma[1]))),
  cosy = list(f = bf(y ~ x + cosy(t, g)),
              extra = list(list(cosy = 0.4), list(cosy = 0.1)),
              cov = function(p) grp_cov(d, cosy1(p$cosy))))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  inits <- base(cs$extra)
  bfit <- ll_fixed_fit(cs$f, student(), d, inits)
  ref <- t(vapply(inits, function(p) {
    grouped_cond(d, mu_of(p, d), p$sigma^2 * cs$cov(p), p$nu)
  }, numeric(n)))
  report(paste(nm, "student"), bfit, ref)
}
inits <- list(list(b = array(0.4, 1), Intercept = 1.1, ar = array(0.5, 1),
                   b_sigma = array(0.3, 1), Intercept_sigma = -0.2, nu = 4),
              list(b = array(0.6, 1), Intercept = 0.9, ar = array(-0.4, 1),
                   b_sigma = array(-0.2, 1), Intercept_sigma = 0.1, nu = 9))
bfit <- ll_fixed_fit(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), student(),
                     d, inits)
ref <- t(vapply(inits, function(p) {
  s <- exp(p$Intercept_sigma + p$b_sigma[1] * (d$x - mean(d$x)))
  grouped_cond(d, mu_of(p, d), outer(s, s) * grp_cov(d, ar1(p$ar[1])), p$nu)
}, numeric(n)))
report("ar sigma~x student", bfit, ref)
du <- d[!(d$g == 1 & d$t == 3) & !(d$g == 4 & d$t == 2), ]
rownames(du) <- NULL
C <- matrix(0.3, 6, 6) + diag(0.7, 6)
C[1, 2] <- C[2, 1] <- 0.6
Lc <- t(chol(C))
inits <- list(list(b = array(0.4, 1), Intercept = 1.1, sigma = 0.8,
                   Lcortime = Lc, nu = 5),
              list(b = array(0.6, 1), Intercept = 0.9, sigma = 1.2,
                   Lcortime = diag(6), nu = 12))
bfit <- ll_fixed_fit(bf(y ~ x + unstr(t, g)), student(), du, inits)
ref <- t(vapply(inits, function(p) {
  CC <- tcrossprod(p$Lcortime)
  grouped_cond(du, mu_of(p, du),
               p$sigma^2 * grp_cov(du, function(k, tt) CC[tt, tt]), p$nu)
}, numeric(nrow(du))))
report("unstr gap student", bfit, ref)

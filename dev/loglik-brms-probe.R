# Lane loglik: what brms 2.23.0's log_lik() returns for every residual
# covariance term, against the leave-one-out conditional computed
# directly from the joint covariance (ll_dense_cond(), no precision
# matrix), at two chosen parameter vectors per case (Fixed_param, one
# chain each). Establishes brms's definition before frmtmb copies it.
# Seed 1. Output: dev/loglik-brms-probe.txt.
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
V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.2, n)
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
dg <- 1 / sqrt(rowSums(P))
W <- P * outer(dg, dg)
d2 <- list(V = V, W = W)

# brms's joint covariance (or scale matrix) of each group, by position
# within the group in time order, at brms's own parameters
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
unstr1 <- function(C) function(k, tt) C[tt, tt]

mu_of <- function(p, dd) p$Intercept + p$b[1] * (dd$x - mean(dd$x))
report <- function(nm, bfit, ref) {
  ll <- log_lik(bfit)
  dif <- max(abs(ll - ref)) / max(abs(ref))
  cat(sprintf("%-26s draws %d  rows %d  max|brms - dense| / max|dense| = %.3e\n",
              nm, nrow(ll), ncol(ll), dif))
  # the row sum is not the joint density
  invisible(ll)
}
pars2 <- function(extra) {
  list(c(list(b = array(0.4, 1), Intercept = 1.1, sigma = 0.8), extra[[1]]),
       c(list(b = array(0.6, 1), Intercept = 0.9, sigma = 1.2), extra[[2]]))
}
joint_ref <- function(dd, inits, covfun, nu = NULL, mean_fun = NULL) {
  t(vapply(inits, function(p) {
    m <- if (is.null(mean_fun)) mu_of(p, dd) else mean_fun(p, dd)
    S <- p$sigma^2 * covfun(p, dd)
    ll_dense_cond(dd$y, m, S, if (is.null(nu)) NULL else p$nu)
  }, numeric(nrow(dd))))
}

cases <- list(
  ar = list(f = bf(y ~ x + ar(t, g, cov = TRUE)),
            extra = list(list(ar = array(0.5, 1)), list(ar = array(-0.3, 1))),
            cov = function(p, dd) grp_cov(dd, ar1(p$ar[1]))),
  ma = list(f = bf(y ~ x + ma(t, g, cov = TRUE)),
            extra = list(list(ma = array(0.4, 1)), list(ma = array(-0.6, 1))),
            cov = function(p, dd) grp_cov(dd, ma1(p$ma[1]))),
  arma = list(f = bf(y ~ x + arma(t, g, cov = TRUE)),
              extra = list(list(ar = array(0.5, 1), ma = array(0.3, 1)),
                           list(ar = array(-0.2, 1), ma = array(0.5, 1))),
              cov = function(p, dd) grp_cov(dd, arma1(p$ar[1], p$ma[1]))),
  cosy = list(f = bf(y ~ x + cosy(t, g)),
              extra = list(list(cosy = 0.4), list(cosy = 0.1)),
              cov = function(p, dd) grp_cov(dd, cosy1(p$cosy))),
  fcor = list(f = bf(y ~ x + fcor(V)),
              extra = list(list(), list()),
              cov = function(p, dd) V),
  lag = list(f = bf(y ~ x + sar(W)),
             extra = list(list(lagsar = 0.4), list(lagsar = -0.3)),
             cov = function(p, dd) {
               A <- diag(n) - p$lagsar * W
               solve(crossprod(A))
             },
             mean = function(p, dd) {
               as.vector(solve(diag(n) - p$lagsar * W, mu_of(p, dd)))
             }),
  err = list(f = bf(y ~ x + sar(W, type = "error")),
             extra = list(list(errorsar = 0.4), list(errorsar = -0.3)),
             cov = function(p, dd) {
               A <- diag(n) - p$errorsar * W
               solve(crossprod(A))
             }))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  for (fam in c("gaussian", "student")) {
    inits <- pars2(cs$extra)
    if (fam == "student") {
      inits[[1]]$nu <- 5
      inits[[2]]$nu <- 12
    }
    bfit <- tryCatch(ll_fixed_fit(cs$f, brmsfamily(fam), d, inits,
                                  data2 = d2),
                     error = function(e) e)
    if (inherits(bfit, "error")) {
      cat(nm, fam, "FIT ERROR", conditionMessage(bfit), "\n")
      next
    }
    ref <- joint_ref(d, inits, cs$cov,
                     nu = if (fam == "student") TRUE, mean_fun = cs$mean)
    ll <- tryCatch(report(paste(nm, fam), bfit, ref),
                   error = function(e) {
                     cat(sprintf("%-26s brms log_lik ERROR: %s\n",
                                 paste(nm, fam), conditionMessage(e)))
                     NULL
                   })
  }
}

# sigma ~ x with ar(cov = TRUE): brms's normal_time_het
inits <- list(list(b = array(0.4, 1), Intercept = 1.1, ar = array(0.5, 1),
                   b_sigma = array(0.3, 1), Intercept_sigma = -0.2),
              list(b = array(0.6, 1), Intercept = 0.9, ar = array(-0.4, 1),
                   b_sigma = array(-0.2, 1), Intercept_sigma = 0.1))
for (fam in c("gaussian", "student")) {
  ii <- inits
  if (fam == "student") {
    ii[[1]]$nu <- 4
    ii[[2]]$nu <- 9
  }
  bfit <- ll_fixed_fit(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
                       brmsfamily(fam), d, ii)
  ref <- t(vapply(ii, function(p) {
    s <- exp(p$Intercept_sigma + p$b_sigma[1] * (d$x - mean(d$x)))
    S <- outer(s, s) * grp_cov(d, ar1(p$ar[1]))
    ll_dense_cond(d$y, mu_of(p, d), S, p$nu)
  }, numeric(n)))
  report(paste("ar sigma~x", fam), bfit, ref)
}

# unstr() with a missing time point: brms carries Jtime_tg
du <- d[!(d$g == 1 & d$t == 3) & !(d$g == 4 & d$t == 2), ]
rownames(du) <- NULL
C <- matrix(0.3, 6, 6) + diag(0.7, 6)
C[1, 2] <- C[2, 1] <- 0.6
Lc <- t(chol(C))
for (fam in c("gaussian", "student")) {
  ii <- list(list(b = array(0.4, 1), Intercept = 1.1, sigma = 0.8,
                  Lcortime = Lc),
             list(b = array(0.6, 1), Intercept = 0.9, sigma = 1.2,
                  Lcortime = diag(6)))
  if (fam == "student") {
    ii[[1]]$nu <- 5
    ii[[2]]$nu <- 12
  }
  bfit <- ll_fixed_fit(bf(y ~ x + unstr(t, g)), brmsfamily(fam), du, ii)
  ref <- t(vapply(ii, function(p) {
    CC <- tcrossprod(p$Lcortime)
    S <- p$sigma^2 * grp_cov(du, function(k, tt) CC[tt, tt])
    ll_dense_cond(du$y, mu_of(p, du), S, p$nu)
  }, numeric(nrow(du))))
  report(paste("unstr gap", fam), bfit, ref)
}

# ar(cov = TRUE) with an interior gap: brms reads the lag by POSITION in
# the group (no Jtime), so the gap is no gap
inits <- pars2(list(list(ar = array(0.5, 1)), list(ar = array(-0.3, 1))))
bfit <- ll_fixed_fit(bf(y ~ x + ar(t, g, cov = TRUE)), gaussian(), du, inits)
pos <- joint_ref(du, inits, function(p, dd) grp_cov(dd, ar1(p$ar[1])))
lev <- joint_ref(du, inits, function(p, dd) grp_cov(dd, function(k, tt) {
  p$ar[1]^abs(outer(tt, tt, "-")) / (1 - p$ar[1]^2)
}))
ll <- log_lik(bfit)
cat(sprintf("ar gap gaussian: brms vs by-position %.3e, vs by-level %.3e\n",
            max(abs(ll - pos)) / max(abs(pos)),
            max(abs(ll - lev)) / max(abs(lev))))

# the row sum of the leave-one-out conditionals is not the joint density
p <- inits[[1]]
S <- p$sigma^2 * grp_cov(d, ar1(p$ar[1]))
m <- mu_of(p, d)
jl <- -0.5 * n * log(2 * pi) - 0.5 * as.numeric(determinant(S)$modulus) -
  0.5 * sum((d$y - m) * solve(S, d$y - m))
cat(sprintf("ar gaussian draw 1: sum of conditionals %.6f, joint %.6f\n",
            sum(ll_dense_cond(d$y, m, S)), jl))

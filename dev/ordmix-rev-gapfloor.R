# Reviewer of lane ordmix: what the interior-category gap floor changes.
# (1) ord_log_interior() with and without the floor on random threshold
#     pairs (seed 1), in ulps, and at a == b.
# (2) Fits where the floor binds: r_theta3 of dev/ordmix-rev-lpcheck.R
#     (data seed 20261005 + 33) and the lane's dev/ordmix-dbg4.R case.
#     At the reported optimum, the log likelihood recomputed in plain
#     doubles from the component families WITHOUT the floor (a log-sum-
#     exp that lets one component be -Inf), against -obj$fn.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
oli <- frmtmb:::ord_log_interior
ulp <- function(a, b) {
  ifelse(a == b, 0, abs(a - b) / (.Machine$double.eps * pmax(abs(a), abs(b))))
}
set.seed(1)
N <- 1e5
a <- rnorm(N, 0, 3)
gap <- 10^runif(N, -15, 1)
b <- a - gap
u <- ulp(oli(a, b, FALSE), oli(a, b, TRUE))
cat(sprintf("(1) %d random pairs, gap 1e-15..10: floor vs none max ulp %.2f, median %.2f, identical %d\n",
            N, max(u), median(u), sum(u == 0)))
cat("    a == b: none", oli(0.41, 0.41, FALSE), " floor",
    oli(0.41, 0.41, TRUE), "\n")
cat("    gap 1e-300 at 0: none", oli(1e-300, 0, FALSE), " floor",
    oli(1e-300, 0, TRUE), "\n")

# plain-double mixture log likelihood without the floor
ll_nofloor <- function(fit) {
  fam <- fit$spec$responses[[1]]$family
  mx <- fam$mix$ord
  rn <- fit$spec$responses[[1]]$resp_name
  est <- fit$estimates
  y <- fit$frame$y[[rn]] %||% fit$frame$y
  lin <- function(dp) {
    lp <- fit$frame$linpreds[[frmtmb:::linpred_key(rn, dp)]]
    if (is.null(lp)) return(NULL)
    if (!is.null(lp$constant)) return(rep(lp$constant, length(y)))
    X <- as.matrix(lp$X)
    if (!ncol(X)) return(rep(0, length(y)))
    as.vector(X %*% est[[lp$par]][lp$idx])
  }
  K <- mx$K
  th <- lapply(seq_len(K), function(k) {
    if (k == fam$mix$ref) 0 * y else lin(paste0("theta", k))
  })
  M <- do.call(cbind, th)
  lpi <- M - apply(M, 1, function(r) max(r) + log(sum(exp(r - max(r)))))
  aterms <- fit$frame$aterms[[rn]] %||% fit$frame$aterms %||% list()
  L <- sapply(seq_len(K), function(k) {
    cp <- mx$comps[[k]]
    dps <- list(mu = lin(paste0("mu", k)))
    dk <- lin(paste0("disc", k))
    if (!is.null(dk)) dps$disc <- exp(dk)
    as.numeric(cp$lpdf(y, dps, aterms,
                       list(tau_raw = est[[mx$tau_names[k]]]))) + lpi[, k]
  })
  m <- apply(L, 1, max)
  sum(m + log(rowSums(exp(L - m))))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

report <- function(tag, fit) {
  mx <- fit$spec$responses[[1]]$family$mix$ord
  for (k in seq_len(mx$K)) {
    tau <- frmtmb:::ord_threshold_values(mx$views[[k]],
                                         fit$estimates[[mx$tau_names[k]]])
    cat(sprintf("    %s component %d thresholds %s; equal adjacent: %s\n",
                tag, k, paste(format(tau, digits = 10), collapse = " "),
                any(diff(tau) == 0)))
  }
  f0 <- -fit$obj$fn(fit$opt$par)
  f1 <- ll_nofloor(fit)
  cat(sprintf("    %s -obj$fn %.12f, no-floor plain doubles %.12f, diff %.3g\n",
              tag, f0, f1, f0 - f1))
}

## r_theta3
set.seed(20261005 + 33)
n <- 400
x <- rnorm(n); z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
yh <- ifelse(runif(n) < 0.2, 0L, y)
w <- runif(n, 0.5, 2)
cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x, z)
ft <- suppressWarnings(frm(bf(y ~ x, theta1 ~ z, theta2 ~ 1),
                           family = mixture(cumulative(), cumulative(),
                                            sratio()), data = d,
                           control = frmtmb_control(grad_tol = 1e-8)))
cat("(2a) r_theta3 logLik", format(as.numeric(logLik(ft)), digits = 12),
    "\n")
report("r_theta3", ft)
# the floored fit against starts that keep every gap open
set.seed(3)
best <- -Inf
for (s in 1:12) {
  st <- lapply(ft$frame$par_template, function(v) v + rnorm(length(v), 0, 0.7))
  fs <- tryCatch(suppressWarnings(frm(bf(y ~ x, theta1 ~ z, theta2 ~ 1),
                 family = mixture(cumulative(), cumulative(), sratio()),
                 data = d, start = st)), error = function(e) NULL)
  if (is.null(fs)) next
  mx <- fs$spec$responses[[1]]$family$mix$ord
  collapsed <- any(vapply(1:2, function(k) {
    any(diff(frmtmb:::ord_threshold_values(mx$views[[k]],
                                           fs$estimates[[mx$tau_names[k]]])) == 0)
  }, NA))
  cat(sprintf("    start %2d logLik %.6f collapsed %s\n", s,
              as.numeric(logLik(fs)), collapsed))
}

## the lane's dbg4 case
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-dbg4.R")
cat("(2b) lane dbg4 source head:\n")
cat(paste0("    ", head(src, 40)), sep = "\n")

## (2c) the falsealarm fits with non-finite standard errors: is a gap
## collapsed (the floor binding) there? gen() of dev/ordmix-rev-falsealarm.R
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-rev-falsealarm.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("cut4", "gen")) eval(e)
}
for (s in c(5, 8, 10, 14, 20, 1, 2, 3)) {
  d <- gen(s)
  f <- suppressWarnings(frm(bf(y | thres(gr = gr) ~ x),
                            family = mixture(cumulative(), cumulative()),
                            data = d))
  se <- suppressWarnings(fixef(f)[, "Est.Error"])
  mx <- f$spec$responses[[1]]$family$mix$ord
  coll <- vapply(1:2, function(k) {
    tau <- frmtmb:::ord_threshold_values(mx$views[[k]],
                                         f$estimates[[mx$tau_names[k]]])
    sum(diff(tau[1:3]) == 0) + sum(diff(tau[4:6]) == 0)
  }, 1)
  cat(sprintf("    gr_cum2 seed %2d logLik %.6f se finite %s collapsed gaps %s min raw log-increment %.1f\n",
              s, as.numeric(logLik(f)), all(is.finite(se)),
              paste(coll, collapse = "+"),
              min(unlist(f$estimates[c("tau_raw1", "tau_raw2")])[-c(1, 4, 7, 10)])))
}

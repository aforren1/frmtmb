# Probe E: categorical emissions - the covid19retrospective shape.
# K = 2 states, C = 4 observed categories, per-state category
# probability vectors, covariate-dependent transitions, many sequences.
#
# Reference: a hand-rolled numeric forward algorithm, plus a direct
# MakeADFun ML over the same likelihood.
#
# Run: FRMTMB_LIB=<lib> Rscript dev/hmm/probeE1-categorical.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
source("dev/hmm/hmm-common.R")

set.seed(555)
K <- 2L
C <- 4L
N <- 40L
Tg <- 15L
n <- N * Tg

# emission: state k has a length-C probability vector, category 1 the
# reference cell of a multinomial logit
E_logit_true <- rbind(c(0, 1.2, 0.2, -1.0),   # state 1
                      c(0, -0.5, 1.4, 1.0))   # state 2
E_true <- t(apply(E_logit_true, 1, function(v) exp(v) / sum(exp(v))))
B_true <- matrix(0, K, K);  B_true[1, 2] <- -1.8; B_true[2, 2] <- 1.3
Bx_true <- matrix(0, K, K); Bx_true[1, 2] <- 1.0; Bx_true[2, 2] <- -0.7

dat <- do.call(rbind, lapply(seq_len(N), function(g) {
  x <- rnorm(Tg)
  s <- integer(Tg)
  s[1] <- sample.int(K, 1, prob = c(0.5, 0.5))
  for (t in seq_len(Tg - 1L)) {
    G <- tpm_tv_num(B_true, Bx_true, x[t], K)
    s[t + 1L] <- sample.int(K, 1, prob = G[s[t], ])
  }
  data.frame(g = g, t = seq_len(Tg), x = x, state = s,
             y = vapply(s, function(k) sample.int(C, 1, prob = E_true[k, ]),
                        integer(1)))
}))
rows_by_g <- split(seq_len(n), dat$g)

cat("== probe E1: categorical HMM, K =", K, " C =", C,
    " N =", N, " T =", Tg, "==\n")
cat("   emission probabilities (true):\n")
print(round(E_true, 4))
cat("   observed category counts:", table(dat$y), "\n\n")

## ---- the custom family -----------------------------------------------
#
# dpars: mu is state 1's logit for category 2 (the primary predictor),
# e13/e14 the rest of state 1, e22/e23/e24 state 2, tr12/tr22 the
# transition logits. Every one of them takes a full formula.

cat_lpdf <- function(y, dpars, aterms) {
  "c" <- RTMB::ADoverload("c")
  K <- 2L
  C <- 4L
  n <- length(y)
  rows_by_g <- hmm_seq_index(aterms$vint1, aterms$vint2)

  # per-state log emission probability AT THE OBSERVED CATEGORY:
  # softmax numerator picked out by data indicators, denominator shared
  ind <- lapply(seq_len(C), function(cc) as.numeric(y == cc))
  logit1 <- list(0, dpars$mu, dpars$e13, dpars$e14)
  logit2 <- list(0, dpars$e22, dpars$e23, dpars$e24)
  emis <- function(lg) {
    tot <- 1
    for (cc in seq_len(C - 1L) + 1L) tot <- tot + exp(lg[[cc]])
    num <- 0
    for (cc in seq_len(C)) num <- num + ind[[cc]] * lg[[cc]]
    num - log(tot)
  }
  lp <- list(emis(logit1), emis(logit2))

  eta <- list(list(0, dpars$tr12), list(0, dpars$tr22))
  lg <- tpm_logs_ad(eta, K)
  ld <- log(c(0.5, 0.5))   # fixed uniform initial distribution

  llv <- NULL
  for (gi in seq_along(rows_by_g)) {
    v <- fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], ld, K)
    llv <- if (is.null(llv)) v else c(llv, v)
  }
  first <- vapply(rows_by_g, function(r) r[1], integer(1))
  S <- Matrix::sparseMatrix(i = first, j = seq_along(first), x = 1,
                            dims = c(n, length(first)))
  as.vector(S %*% llv)
}

dp_names <- c("mu", "e13", "e14", "e22", "e23", "e24", "tr12", "tr22")
fam <- custom_family(
  "hmm_categorical", dpars = dp_names,
  links = stats::setNames(rep(list("identity"), length(dp_names)),
                          dp_names),
  lpdf = cat_lpdf, type = "discrete",
  init_dpars = stats::setNames(
    lapply(c(0.3, 0, -0.3, -0.3, 0.3, 0, -1.5, 1.5),
           function(v) { force(v); function(y, aterms) v }),
    dp_names)
)

## ---- 1. AD safety ----------------------------------------------------

td <- stats::setNames(lapply(c(1.2, 0.2, -1, -0.5, 1.4, 1, -1.8, 1.3),
                             function(v) rep(v, n)), dp_names)
chk <- tryCatch(check_custom_family(fam, y = dat$y, dpars = td,
                                    aterms = list(vint1 = dat$g,
                                                  vint2 = dat$t)),
                error = function(e) conditionMessage(e))
cat("1. check_custom_family:",
    if (isTRUE(chk)) "PASS" else paste("FAIL -", chk), "\n")

## numeric reference at the same point
lpm_num <- function(El, y) {
  vapply(seq_len(K), function(k) {
    p <- exp(El[k, ]) / sum(exp(El[k, ]))
    log(p[y])
  }, numeric(length(y)))
}
ll_ref0 <- sum(vapply(rows_by_g, function(r)
  fwd_num_tv(lpm_num(E_logit_true, dat$y),
             function(rr) tpm_tv_num(B_true, Bx_true, dat$x[rr], K),
             r, c(0.5, 0.5)), numeric(1)))
td2 <- td
td2$tr12 <- B_true[1, 2] + Bx_true[1, 2] * dat$x
td2$tr22 <- B_true[2, 2] + Bx_true[2, 2] * dat$x
ll_fam0 <- sum(cat_lpdf(dat$y, td2, list(vint1 = dat$g, vint2 = dat$t)))
cat("   family lpdf vs numeric forward at the true parameters:",
    format(abs(ll_fam0 - ll_ref0), digits = 3), "\n\n")

## ---- 2. the fit ------------------------------------------------------

form <- bf(y | vint(g, t) ~ 1, e13 ~ 1, e14 ~ 1, e22 ~ 1, e23 ~ 1,
           e24 ~ 1, tr12 ~ x, tr22 ~ x)
t_fit <- system.time(fit <- frm(form + fam, data = dat))[["elapsed"]]
cat(sprintf("2. frm() in %.2f s  logLik %.8f  df %d\n", t_fit,
            as.numeric(logLik(fit)), attr(logLik(fit), "df")))
fx <- fixef(fit)
El_hat <- rbind(c(0, fx$mu[[1]], fx$e13[[1]], fx$e14[[1]]),
                c(0, fx$e22[[1]], fx$e23[[1]], fx$e24[[1]]))
E_hat <- t(apply(El_hat, 1, function(v) exp(v) / sum(exp(v))))
cat("   fitted emission probabilities:\n"); print(round(E_hat, 4))
cat("   true:\n"); print(round(E_true, 4))
cat("   max |E_hat - E_true| :", format(max(abs(E_hat - E_true)),
                                        digits = 3), "\n")
cat("   transitions: tr12", format(c(fx$tr12[[1]], fx$tr12[[2]]),
                                   digits = 5),
    " (true", B_true[1, 2], Bx_true[1, 2], ")\n")
cat("                tr22", format(c(fx$tr22[[1]], fx$tr22[[2]]),
                                   digits = 5),
    " (true", B_true[2, 2], Bx_true[2, 2], ")\n\n")

## ---- 3. direct ML over the same likelihood ---------------------------

Bh <- matrix(0, K, K);  Bh[1, 2] <- fx$tr12[[1]]; Bh[2, 2] <- fx$tr22[[1]]
Bxh <- matrix(0, K, K); Bxh[1, 2] <- fx$tr12[[2]]; Bxh[2, 2] <- fx$tr22[[2]]
ll_at_fit <- sum(vapply(rows_by_g, function(r)
  fwd_num_tv(lpm_num(El_hat, dat$y),
             function(rr) tpm_tv_num(Bh, Bxh, dat$x[rr], K),
             r, c(0.5, 0.5)), numeric(1)))
cat("3. numeric forward at frm's estimates:",
    format(ll_at_fit, digits = 12), "\n")
cat("   frm logLik                      :",
    format(as.numeric(logLik(fit)), digits = 12), "\n")
cat("   |diff|                          :",
    format(abs(ll_at_fit - as.numeric(logLik(fit))), digits = 3), "\n\n")

## independent direct optimization (BFGS on the numeric forward)
dnll <- function(v) {
  El <- rbind(c(0, v[1:3]), c(0, v[4:6]))
  Bd <- matrix(0, K, K);  Bd[1, 2] <- v[7]; Bd[2, 2] <- v[8]
  Bxd <- matrix(0, K, K); Bxd[1, 2] <- v[9]; Bxd[2, 2] <- v[10]
  -sum(vapply(rows_by_g, function(r)
    fwd_num_tv(lpm_num(El, dat$y),
               function(rr) tpm_tv_num(Bd, Bxd, dat$x[rr], K), r,
               c(0.5, 0.5)), numeric(1)))
}
v0 <- c(El_hat[1, 2:4], El_hat[2, 2:4], Bh[1, 2], Bh[2, 2],
        Bxh[1, 2], Bxh[2, 2])
od <- optim(v0, dnll, method = "BFGS",
            control = list(reltol = 1e-13, maxit = 2000))
cat("4. independent BFGS on the numeric forward: logLik",
    format(-od$value, digits = 12), "\n")
cat("   frm logLik                              :",
    format(as.numeric(logLik(fit)), digits = 12), "\n")
cat("   |diff|                                  :",
    format(abs(-od$value - as.numeric(logLik(fit))), digits = 3), "\n")
cat("   max |parameter diff|                    :",
    format(max(abs(od$par - v0)), digits = 3), "\n\n")

## ---- 5. state decoding -----------------------------------------------

lpm <- lpm_num(El_hat, dat$y)
fbv <- function(rows) {
  Tl <- length(rows)
  la <- matrix(0, Tl, K); lb <- matrix(0, Tl, K)
  a <- c(0.5, 0.5) * exp(lpm[rows[1], ]); s <- sum(a); lsc <- log(s)
  la[1, ] <- log(a / s) + lsc
  for (k in seq_len(Tl - 1L) + 1L) {
    G <- tpm_tv_num(Bh, Bxh, dat$x[rows[k - 1L]], K)
    a <- as.vector((a / s) %*% G) * exp(lpm[rows[k], ])
    s <- sum(a); lsc <- lsc + log(s)
    la[k, ] <- log(a / s) + lsc
  }
  for (k in seq(Tl - 1L, 1L)) {
    G <- tpm_tv_num(Bh, Bxh, dat$x[rows[k]], K)
    lb[k, ] <- log(as.vector(G %*% exp(lpm[rows[k + 1L], ] +
                                         lb[k + 1L, ])))
  }
  lg <- la + lb
  ex <- exp(lg - apply(lg, 1, max))
  ex / rowSums(ex)
}
sp <- matrix(0, n, K)
for (r in rows_by_g) sp[r, ] <- fbv(r)
cat("5. local decoding accuracy vs the simulated states:",
    format(mean(max.col(sp) == dat$state), digits = 4), "\n")
cat("   (a categorical emission with overlapping category profiles is\n",
    "    much less identifiable than a well-separated gaussian one)\n")

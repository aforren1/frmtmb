# Reviewer, item 2.3, claim 1: hmm_starts() with the guarded thing
# ABSENT, and in the four other corners the suite does not construct.
#
# Cases, in order:
#   A. the better mode is ABSENT           (incumbent is already best)
#   B. the incumbent is best AND refits converge cleanly to something
#      WORSE                               (the warm d4 fit)
#   C. a refit reaches the better optimum but fails the convergence
#      test                                (grad_tol driven to 1e-12)
#   D. every refit errors                  (jitter 1e6)
#   E. the reported spread is a spread of OPTIMA, checked by refitting
#      from the winning start that was written into $best$call
#
#   Rscript dev/rev-latent-guard2.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}
d4_data <- function() {
  set.seed(2026)
  K <- 2L; N <- 25L; Tg <- 30L
  G <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
  mu <- c(0, 3); sg <- c(0.6, 0.6); sd_b <- c(0.7, 0.5)
  stat <- local({
    A <- rbind(t(diag(K) - G), 1)
    drop(qr.solve(A, c(rep(0, K), 1)))
  })
  b <- cbind(stats::rnorm(N, 0, sd_b[1L]), stats::rnorm(N, 0, sd_b[2L]))
  d <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg); s[1L] <- sample.int(K, 1L, prob = stat)
    for (t in seq_len(Tg - 1L)) s[t + 1L] <- sample.int(K, 1L, prob = G[s[t], ])
    data.frame(ID = g, t = seq_len(Tg),
               y = stats::rnorm(Tg, mu[s] + b[g, s], sg[s]))
  }))
  d$gf <- factor(d$ID); d
}
d4_form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
d4_fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")

show <- function(ms, label) {
  cat("--", label, "--\n")
  cat("  n_converged/not/error:", ms$n_converged, ms$n_not_converged,
      ms$n_error, "  sum ==", ms$n, ":",
      ms$n_converged + ms$n_not_converged + ms$n_error == ms$n, "\n")
  cat("  original_converged:", ms$original_converged,
      "  scale_how:", ms$scale_how, "\n")
  cat("  spread    :", format(ms$spread, digits = 10), "\n")
  cat("  spread_all:", format(ms$spread_all, digits = 10), "\n")
  cat("  modes     :", nrow(ms$modes), "rows\n")
  print(ms$table, row.names = FALSE, digits = 12)
  cat("  best logLik    :", format(as.numeric(logLik(ms$best)), digits = 14),
      "\n")
  cat("  original logLik:", format(ms$original_logLik, digits = 14), "\n")
  cat("  best IS the original object:", identical(ms$best, NULL), "\n")
}

## ---- A: the better mode is ABSENT ---------------------------------
hr("A. the better mode is ABSENT")
dd <- small_data()
fa <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
          data = dd)
lla <- as.numeric(logLik(fa))
msA <- suppressWarnings(hmm_starts(fa, n = 8, jitter = 2, seed = 8801))
show(msA, "small_data, n = 8, jitter = 2, seed 8801")
cat("  best beats the original by:",
    format(as.numeric(logLik(msA$best)) - lla, digits = 6), "\n")
cat("  best is the SAME fit object as the incumbent:",
    identical(as.numeric(logLik(msA$best)), lla), "\n")
outA <- utils::capture.output(print(msA))
cat("  print says 'the original fit was the best found':",
    any(grepl("the original fit was the best found", outA)), "\n")
cat("  print says 'BELOW the best':", any(grepl("BELOW the best", outA)), "\n")
cat("--- print(), verbatim ---\n"); cat(paste(outA, collapse = "\n"), "\n")

## ---- B: incumbent best, refits converge cleanly to something WORSE --
hr("B. incumbent is best, refits converge cleanly to a WORSE optimum")
d4 <- d4_data()
tpl <- par_template(d4_form, family = d4_fam, data = d4)
st <- tpl
st$beta[["mu1_(Intercept)"]] <- -0.185185
st$beta[["mu2_(Intercept)"]] <- 3.121877
st$betad[["sigma1_(Intercept)"]] <- log(0.610455)
st$betad[["sigma2_(Intercept)"]] <- log(0.602054)
st$betad[["tr12_(Intercept)"]] <- log(0.173438 / 0.826562)
st$betad[["tr22_(Intercept)"]] <- log(0.828240 / 0.171760)
st$theta <- log(c(0.645297, 0.427590))
fw <- frm(d4_form, family = d4_fam, data = d4, start = st)
cat("warm incumbent logLik:", format(as.numeric(logLik(fw)), digits = 14),
    "\n")
msB <- suppressWarnings(hmm_starts(fw, n = 8, jitter = 2, seed = 4101))
show(msB, "warm d4 incumbent, n = 8, jitter = 2, seed 4101")
worse <- msB$table[msB$table$status == "converged" &
                     msB$table$logLik < msB$original_logLik - 1e-6, ]
cat("  CONVERGED refits strictly worse than the incumbent:", nrow(worse),
    "\n")
if (nrow(worse)) {
  cat("  their logLiks:",
      paste(format(worse$logLik, digits = 12), collapse = " "), "\n")
  cat("  their grad_rel:",
      paste(format(worse$grad_rel, digits = 3), collapse = " "), "\n")
  cat("  worst of them is below the incumbent by:",
      format(msB$original_logLik - min(worse$logLik), digits = 8), "\n")
}
cat("  best stayed the incumbent:",
    isTRUE(all.equal(as.numeric(logLik(msB$best)), msB$original_logLik)), "\n")
outB <- utils::capture.output(print(msB))
cat("  print says 'the original fit was the best found':",
    any(grepl("the original fit was the best found", outB)), "\n")
cat("--- print(), verbatim ---\n"); cat(paste(outB, collapse = "\n"), "\n")

## ---- C: a refit finds the better optimum but is called not-converged
hr("C. the better optimum found, but failing the convergence test")
fc <- frm(d4_form, family = d4_fam, data = d4)
llc <- as.numeric(logLik(fc))
msC <- suppressWarnings(hmm_starts(fc, n = 6, jitter = 2, seed = 4101,
                                   grad_tol = 1e-12))
show(msC, "cold d4, grad_tol = 1e-12, n = 6, jitter = 2, seed 4101")
hi <- max(msC$table$logLik[msC$table$status == "not converged"],
          na.rm = TRUE)
cat("  highest not-converged logLik:", format(hi, digits = 14), "\n")
cat("  best logLik                 :",
    format(as.numeric(logLik(msC$best)), digits = 14), "\n")
cat("  a NON-CONVERGED refit that is", format(hi - llc, digits = 6),
    "units better did NOT win:", as.numeric(logLik(msC$best)) < hi, "\n")
cat("  converged spread reads:", format(msC$spread, digits = 6),
    " <- with only the incumbent in it, this is a printed zero\n")
cat("  spread_all reads      :", format(msC$spread_all, digits = 6), "\n")
outC <- utils::capture.output(print(msC))
cat("--- print(), verbatim ---\n"); cat(paste(outC, collapse = "\n"), "\n")

## ---- D: every refit errors ----------------------------------------
hr("D. every refit errors")
msD <- suppressWarnings(hmm_starts(fa, n = 4, jitter = 1e6, seed = 3L))
show(msD, "small_data, jitter = 1e6, n = 4")
cat("  all four errored:", msD$n_error == 4L, "\n")
cat("  spread / spread_all:", msD$spread, "/", msD$spread_all,
    "  <- both read 0 on a run that measured NOTHING\n")
cat("  modes rows:", nrow(msD$modes), "\n")
outD <- utils::capture.output(print(msD))
cat("--- print(), verbatim ---\n"); cat(paste(outD, collapse = "\n"), "\n")

## ---- E: is the spread a spread of OPTIMA? --------------------------
hr("E. the spread is a spread of optima, checked by re-evaluation")
msE <- suppressWarnings(hmm_starts(fc, n = 6, jitter = 2, seed = 4101))
show(msE, "cold d4, n = 6, jitter = 2, seed 4101")
cat("  spread == diff(range(converged + original)):",
    isTRUE(all.equal(msE$spread,
                     diff(range(msE$table$logLik[
                       msE$table$status %in% c("original", "converged")])))),
    "\n")
stw <- msE$best$call[["start"]]
cat("  winning start is a list by value:", is.list(stw), "\n")
if (is.list(stw)) {
  # the objective AT the winning start, which is what a spread of
  # STARTING values would have reported
  o0 <- frm(d4_form, family = d4_fam, data = d4, start = stw,
            dry_run = "objective")
  ll_at_start <- -o0$fn(o0$par)
  cat("  logLik AT the winning start :", format(ll_at_start, digits = 12),
      "\n")
  cat("  logLik at the optimum it reached:",
      format(as.numeric(logLik(msE$best)), digits = 12), "\n")
  cat("  the two differ by:",
      format(as.numeric(logLik(msE$best)) - ll_at_start, digits = 6),
      " <- a spread of starts is not this spread\n")
  again <- frm(d4_form, family = d4_fam, data = d4, start = stw)
  cat("  refitting from the recorded start reproduces $best to:",
      format(abs(as.numeric(logLik(again)) -
                   as.numeric(logLik(msE$best))), digits = 3), "\n")
}

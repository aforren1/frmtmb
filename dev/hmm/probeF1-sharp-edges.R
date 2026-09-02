# Probe F: the sharp-edge hunt.
#
#  F1  label switching and init sensitivity (mixture()'s quantile
#      convention as the model)
#  F2  sequences of length 1: the HMM degenerates to a finite mixture,
#      so it must reproduce mixture(gaussian(), gaussian()) exactly
#  F3  stationary vs free initial distribution
#  F4  a missing observation mid-sequence
#  F5  REML with HMM extras
#
# Run: FRMTMB_LIB=<lib> Rscript dev/hmm/probeF1-sharp-edges.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
source("dev/hmm/hmm-common.R")
source("dev/hmm/hmm-family.R")

K <- 2L
cat("== probe F: sharp edges ==\n\n")

## ---- F1. label switching and init sensitivity ------------------------

ref <- readRDS("dev/hmm/probeB1.rds")
dat <- ref$dat
form <- bf(y | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
           tr12 ~ x, tr22 ~ x)

# mixture()'s own convention is quantile(y, k / (K + 1)): 1/3 and 2/3
mk_fam <- function(q1, q2, tr1 = -1.5, tr2 = 1.5) {
  frmtmb::custom_family(
    "hmm2f", dpars = c("mu", "mu2", "sigma1", "sigma2", "tr12", "tr22"),
    links = list(mu = "identity", mu2 = "identity", sigma1 = "log",
                 sigma2 = "log", tr12 = "identity", tr22 = "identity"),
    lpdf = hmm2_lpdf, type = "continuous",
    init_dpars = list(
      mu = function(y, aterms) unname(stats::quantile(y, q1)),
      mu2 = function(y, aterms) unname(stats::quantile(y, q2)),
      sigma1 = function(y, aterms) stats::sd(y),
      sigma2 = function(y, aterms) stats::sd(y),
      tr12 = function(y, aterms) tr1,
      tr22 = function(y, aterms) tr2),
    extra_pars = function(y, aterms) list(hmm_ldel = 0))
}

cat("F1. init sensitivity and label switching (probe B1 data,",
    "reference logLik", format(ref$ll, digits = 10), ")\n")
cat(sprintf("   %-26s %16s %10s %10s\n", "init (mu, mu2 quantiles)",
            "logLik", "mu1", "mu2"))
inits <- list(c(1 / 3, 2 / 3), c(0.25, 0.75), c(0.1, 0.9),
              c(2 / 3, 1 / 3), c(0.5, 0.5), c(0.45, 0.55))
for (q in inits) {
  f <- try(frm(form + mk_fam(q[1], q[2]), data = dat), silent = TRUE)
  if (inherits(f, "try-error")) {
    cat(sprintf("   %-26s %16s\n", paste(q, collapse = ", "), "ERROR"))
    next
  }
  fx <- fixef(f)
  cat(sprintf("   %-26s %16.8f %10.4f %10.4f\n",
              paste(format(q, digits = 3), collapse = ", "),
              as.numeric(logLik(f)), fx$mu[[1]], fx$mu2[[1]]))
}
# a degenerate start: both means at the median AND symmetric transitions
f_deg <- try(frm(form + mk_fam(0.5, 0.5, 0, 0), data = dat),
             silent = TRUE)
cat("   fully symmetric start (both means at the median, both\n")
cat("   transition logits 0):",
    if (inherits(f_deg, "try-error")) {
      paste("ERROR -", substr(attr(f_deg, "condition")$message, 1, 60))
    } else {
      format(as.numeric(logLik(f_deg)), digits = 10)
    }, "\n\n")

## ---- F2. sequences of length 1 == a finite mixture --------------------

set.seed(808)
n2 <- 400
z <- rbinom(n2, 1, 0.35)
d2 <- data.frame(y = ifelse(z == 1, rnorm(n2, 3, 0.6), rnorm(n2, 0, 1)),
                 g = seq_len(n2), t = 1L)
form2 <- bf(y | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
            tr12 ~ 1, tr22 ~ 1)
f_hmm <- frm(form2 + mk_fam(1 / 3, 2 / 3), data = d2)
f_mix <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = d2)
cat("F2. every sequence of length 1 (n =", n2, ")\n")
cat("   HMM     logLik :", format(as.numeric(logLik(f_hmm)),
                                  digits = 12), "\n")
cat("   mixture logLik :", format(as.numeric(logLik(f_mix)),
                                  digits = 12), "\n")
cat("   |diff|         :",
    format(abs(as.numeric(logLik(f_hmm)) - as.numeric(logLik(f_mix))),
           digits = 3), "\n")
fh <- fixef(f_hmm)
fm <- fixef(f_mix)
cat("   HMM     mu, sigma, weight:",
    format(c(fh$mu[[1]], fh$mu2[[1]], exp(fh$sigma1[[1]]),
             exp(fh$sigma2[[1]]),
             softmax0(f_hmm$estimates[["hmm_ldel"]])[1]),
          digits = 6), "\n")
cat("   mixture mu, sigma, weight:",
    format(c(fm$mu1[[1]], fm$mu2[[1]], exp(fm$sigma1[[1]]),
             exp(fm$sigma2[[1]]), stats::plogis(fm$theta1[[1]])),
          digits = 6), "\n")
cat("   df: HMM", attr(logLik(f_hmm), "df"), " mixture",
    attr(logLik(f_mix), "df"),
    " (the HMM carries two transition logits that a length-1\n",
    "    sequence never uses: they are unidentified)\n\n")

## ---- F3. stationary vs free initial distribution ---------------------

d3 <- readRDS("dev/hmm/probeD1.rds")$dat
form3 <- bf(y | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
            tr12 ~ 1, tr22 ~ 1)
f_free <- frm(form3 + hmm2_family(), data = d3)
f_stat <- frm(form3 + hmm2_family_stat(), data = d3)
cat("F3. initial distribution\n")
cat("   free delta (one extra parameter) : logLik",
    format(as.numeric(logLik(f_free)), digits = 12), " df",
    attr(logLik(f_free), "df"), "\n")
cat("   stationary (on-tape linear solve): logLik",
    format(as.numeric(logLik(f_stat)), digits = 12), " df",
    attr(logLik(f_stat), "df"), "\n")
cat("   2 * (free - stationary) =",
    format(2 * (as.numeric(logLik(f_free)) - as.numeric(logLik(f_stat))),
           digits = 4), "on 1 df\n")
fs <- fixef(f_stat)
Gs <- tpm_tv_num(rbind(c(0, fs$tr12[[1]]), c(0, fs$tr22[[1]])),
                 matrix(0, K, K), 0, K)
cat("   stationary delta implied by the fit:",
    format(stat_dist(Gs), digits = 5), "\n\n")

## ---- F4. a missing observation mid-sequence --------------------------

set.seed(64)
N4 <- 30L; T4 <- 20L; n4 <- N4 * T4
G4 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)
mu4 <- c(0, 3); sg4 <- c(0.6, 0.6)
d4 <- do.call(rbind, lapply(seq_len(N4), function(g) {
  s <- integer(T4); s[1] <- sample.int(K, 1, prob = c(0.5, 0.5))
  for (t in seq_len(T4 - 1L)) s[t + 1L] <- sample.int(K, 1,
                                                      prob = G4[s[t], ])
  data.frame(g = g, t = seq_len(T4),
             y = rnorm(T4, mu4[s], sg4[s]))
}))
miss <- which(d4$t %in% c(7, 8, 14))     # gaps in the middle
d4$miss <- as.integer(seq_len(n4) %in% miss)
d4$y_na <- d4$y; d4$y_na[miss] <- NA
d4$y_sent <- d4$y; d4$y_sent[miss] <- 0  # never read: masked below

# the CORRECT likelihood: the emission term drops out at a missing
# time point, the transition still happens
lp_masked <- function(y, mu, sg, miss) {
  cbind(ifelse(miss == 1, 0, stats::dnorm(y, mu[1], sg[1], log = TRUE)),
        ifelse(miss == 1, 0, stats::dnorm(y, mu[2], sg[2], log = TRUE)))
}
rbg4 <- split(seq_len(n4), d4$g)
ref_ll <- function(v) {
  G <- tpm_tv_num(rbind(c(0, v[5]), c(0, v[6])), matrix(0, K, K), 0, K)
  lpm <- lp_masked(d4$y_sent, v[1:2], exp(v[3:4]), d4$miss)
  sum(vapply(rbg4, function(r) fwd_num_tv(lpm, function(rr) G, r,
                                          softmax0(v[7])), numeric(1)))
}

# route B: keep the row, carry the missing flag in vint3, zero the
# emission with a data multiplier (no branching on parameters)
hmm2_lpdf_miss <- function(y, dpars, aterms,
                           extra = list(hmm_ldel = 0)) {
  "c" <- RTMB::ADoverload("c")
  n <- length(y)
  obs <- 1 - aterms$vint3          # data, so this is a constant weight
  rows_by_g <- hmm_seq_index(aterms$vint1, aterms$vint2)
  lp <- list(obs * RTMB::dnorm(y, dpars$mu, dpars$sigma1, log = TRUE),
             obs * RTMB::dnorm(y, dpars$mu2, dpars$sigma2, log = TRUE))
  eta <- list(list(0, dpars$tr12), list(0, dpars$tr22))
  lg <- tpm_logs_ad(eta, 2L)
  ld <- log(softmax0_ad(extra$hmm_ldel, 2L))
  llv <- NULL
  for (gi in seq_along(rows_by_g)) {
    v <- fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], ld, 2L)
    llv <- if (is.null(llv)) v else c(llv, v)
  }
  first <- vapply(rows_by_g, function(r) r[1], integer(1))
  S <- Matrix::sparseMatrix(i = first, j = seq_along(first), x = 1,
                            dims = c(n, length(first)))
  as.vector(S %*% llv)
}
fam_miss <- custom_family(
  "hmm2_miss", dpars = c("mu", "mu2", "sigma1", "sigma2", "tr12",
                         "tr22"),
  links = list(mu = "identity", mu2 = "identity", sigma1 = "log",
               sigma2 = "log", tr12 = "identity", tr22 = "identity"),
  lpdf = hmm2_lpdf_miss, type = "continuous",
  init_dpars = list(
    mu = function(y, aterms) unname(stats::quantile(y, 1 / 3)),
    mu2 = function(y, aterms) unname(stats::quantile(y, 2 / 3)),
    sigma1 = function(y, aterms) stats::sd(y),
    sigma2 = function(y, aterms) stats::sd(y),
    tr12 = function(y, aterms) -1.5, tr22 = function(y, aterms) 1.5),
  extra_pars = function(y, aterms) list(hmm_ldel = 0))

f_mask <- frm(bf(y_sent | vint(g, t, miss) ~ 1, mu2 ~ 1, sigma1 ~ 1,
                 sigma2 ~ 1, tr12 ~ 1, tr22 ~ 1) + fam_miss, data = d4)
fxm <- fixef(f_mask)
ldm <- f_mask$estimates[["hmm_ldel"]]
v_mask <- c(fxm$mu[[1]], fxm$mu2[[1]], fxm$sigma1[[1]], fxm$sigma2[[1]],
            fxm$tr12[[1]], fxm$tr22[[1]],
            if (is.null(ldm)) 0 else ldm)
cat("F4. a missing observation mid-sequence (", length(miss),
    "of", n4, "rows )\n")
cat("   masked-emission route : logLik",
    format(as.numeric(logLik(f_mask)), digits = 12), "\n")
cat("   numeric reference at the same parameters:",
    format(ref_ll(v_mask), digits = 12), "\n")
cat("   |diff|                :",
    format(abs(ref_ll(v_mask) - as.numeric(logLik(f_mask))),
           digits = 3), "\n")

# route A: let na.omit drop the rows
f_drop <- suppressMessages(
  frm(bf(y_na | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
         tr12 ~ 1, tr22 ~ 1) + hmm2_family(), data = d4))
fxd <- fixef(f_drop)
cat("   na.omit route (rows dropped): logLik",
    format(as.numeric(logLik(f_drop)), digits = 12), " on n =",
    f_drop$frame$n_obs, "rows\n")
cat("   transition estimates, masked vs dropped:\n")
cat("     tr12", format(c(fxm$tr12[[1]], fxd$tr12[[1]]), digits = 6),
    "\n     tr22", format(c(fxm$tr22[[1]], fxd$tr22[[1]]), digits = 6),
    "\n")
Gm <- tpm_tv_num(rbind(c(0, fxm$tr12[[1]]), c(0, fxm$tr22[[1]])),
                 matrix(0, K, K), 0, K)
Gd <- tpm_tv_num(rbind(c(0, fxd$tr12[[1]]), c(0, fxd$tr22[[1]])),
                 matrix(0, K, K), 0, K)
cat("     tpm masked :", format(as.vector(t(Gm)), digits = 5), "\n")
cat("     tpm dropped:", format(as.vector(t(Gd)), digits = 5), "\n")
cat("     tpm true   :", format(as.vector(t(G4)), digits = 5), "\n\n")

## ---- F5. REML --------------------------------------------------------

cat("F5. REML\n")
r1 <- tryCatch(frm(form3 + hmm2_family(), data = d3, REML = TRUE),
               error = function(e) e, warning = function(w) w)
if (inherits(r1, "condition")) {
  cat("   REML on the fixed-effect HMM:",
      substr(conditionMessage(r1), 1, 140), "\n")
} else {
  cat("   REML on the fixed-effect HMM: logLik",
      format(as.numeric(logLik(r1)), digits = 12), " (ML was ",
      format(as.numeric(logLik(f_free)), digits = 12), ")\n")
  cat("   REML integrates only `mu`'s fixed effects (primary_dpars),\n")
  cat("   so the other five dpars' coefficients stay in the outer\n")
  cat("   problem: this is a PARTIAL restricted likelihood.\n")
}
d3$gf <- factor(d3$g)
r2 <- tryCatch(
  frm(bf(y | vint(g, t) ~ 1 + (1 | gf), mu2 ~ 1, sigma1 ~ 1,
         sigma2 ~ 1, tr12 ~ 1, tr22 ~ 1) + hmm2_family(),
      data = d3, REML = TRUE),
  error = function(e) e, warning = function(w) w)
if (inherits(r2, "condition")) {
  cat("   REML + (1 | g):", substr(conditionMessage(r2), 1, 140), "\n")
} else {
  cat("   REML + (1 | g): logLik", format(as.numeric(logLik(r2)),
                                          digits = 12), " df",
      attr(logLik(r2), "df"), "\n")
}
cat("\n")

# Reviewer, punch round 1, attacks 2 and 3.
#
#  A. `hmm_starts_tol(grad_tol, ref) = grad_tol^2 * max(|ref|, 1)`.
#     The argument for squaring holds where the optimum is smooth and
#     interior. Two things it does not cover:
#       A1  two genuine optima closer than the tolerance. Constructed
#           by driving the tolerance UP past a gap that is real: the d4
#           probe's two modes are 8.099 units apart, so any `grad_tol`
#           above sqrt(8.099 / 1096) merges them.
#       A2  an optimum at a boundary, where the gradient does not
#           vanish. Whether `hmm()` HAS such a thing is the question,
#           since every constrained quantity here is on a link.
#     Also measured: the observed same-optimum noise floor, so the
#     tolerance can be priced against what it is protecting from.
#
#  B. The replacement alignment check. It cross-checks confint()'s
#     standard errors against vcov()'s BY NAME on the rows they share,
#     and the comment says a permutation confined to the `theta_` rows
#     would not be caught because those have no second path. Both
#     halves are constructed here on a REAL fit, not a stub.
#
#   Rscript dev/rev-latent-tol.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")
tolf <- getFromNamespace("hmm_starts_tol", "frmtmb.latent")
scalef <- getFromNamespace("hmm_starts_scale", "frmtmb.latent")

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
d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-2p3-repro81.rds")
d4form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
d4fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")

hr("A0. the tolerance, priced against the noise it is protecting from")
fitS <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
            data = small_data())
msS <- suppressWarnings(hmm_starts(fitS, n = 8, jitter = 2, seed = 8801))
llS <- msS$original_logLik
noise <- diff(range(msS$table$logLik[msS$table$status != "error"]))
cat("a unimodal fit at logLik", format(llS, digits = 10), "\n")
cat("  observed spread over 8 refits on ONE optimum:",
    format(noise, digits = 4), "absolute,",
    format(noise / abs(llS), digits = 4), "relative\n")
cat("  hmm_starts_tol(1e-3, ll)                   :",
    format(tolf(1e-3, llS), digits = 4), "absolute,",
    format(tolf(1e-3, llS) / abs(llS), digits = 4), "relative\n")
cat("  headroom, tolerance over observed noise    :",
    format(tolf(1e-3, llS) / noise, digits = 5), "x\n")

hr("A1. two REAL optima merged, by loosening grad_tol")
fit4 <- frm(d4form, family = d4fam, data = d$dat)
ll4 <- as.numeric(logLik(fit4))
cat("d4 cold logLik:", format(ll4, digits = 12),
    "  known better optimum 8.099 units up\n")
cat("break point: grad_tol above sqrt(8.099 / |ll|) =",
    format(sqrt(8.099290821 / abs(ll4)), digits = 5),
    "makes the tolerance exceed the real gap\n\n")
for (gt in c(1e-3, 1e-2, 3e-2, 0.086, 0.1)) {
  ms <- suppressWarnings(hmm_starts(fit4, n = 6, jitter = 2, seed = 4101,
                                    grad_tol = gt))
  out <- paste(utils::capture.output(print(ms)), collapse = " ")
  cat(sprintf("grad_tol %-7s tol %9.4f  modes %d  best %s  says local: %s\n",
              format(gt), tolf(gt, ll4), nrow(ms$modes),
              format(as.numeric(logLik(ms$best)), digits = 12),
              grepl("found a local optimum", out)))
}
cat("\nthe two knobs are tied QUADRATICALLY: a 100x looser grad_tol is\n")
cat("a 10 000x looser merge threshold, and nothing on the help page\n")
cat("says so.\n")

hr("A2. is there a boundary optimum in this family at all?")
cat("every constrained quantity here is on a link: a state sd through\n")
cat("log, a transition through a reference-cell logit, a random-effect\n")
cat("sd through log. The outer coordinates are unbounded, so the\n")
cat("'smooth and interior' precondition is structural rather than\n")
cat("lucky. The nearest thing to a boundary is a variance running to\n")
cat("zero, which is an ASYMPTOTE in theta, not a boundary:\n\n")
set.seed(606)
# a random intercept with no group variation at all: sd(re) -> 0
nz <- do.call(rbind, lapply(1:12, function(g) {
  s <- integer(40); s[1L] <- 1L
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  for (t in 2:40) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
  data.frame(gf = g, t = seq_len(40), y = stats::rnorm(40, c(0, 3)[s], 0.6))
}))
nz$gf <- factor(nz$gf)
fz <- try(suppressWarnings(frm(bf(y ~ 1 + (1 | gf)),
                               family = hmm(K = 2, gaussian(), time = t,
                                            group = gf), data = nz)),
          silent = TRUE)
if (inherits(fz, "try-error")) {
  cat("  frm() refused this fit:",
      substr(attr(fz, "condition")$message, 1, 90), "\n")
} else {
  llz <- as.numeric(logLik(fz))
  dgz <- frmtmb::diagnose(fz, quiet = TRUE)
  cat("  logLik", format(llz, digits = 10), " theta",
      format(fz$opt$par[names(fz$opt$par) == "theta"], digits = 5), "\n")
  cat("  sd(re) on the natural scale:",
      format(sqrt(VarCorr(fz)[[1L]][1L, 1L]), digits = 4), "\n")
  cat("  max|grad| / |logLik|:", format(dgz$max_grad / abs(llz),
                                        digits = 4),
      "  <- small, so the gradient DOES vanish here\n")
  msz <- suppressWarnings(hmm_starts(fz, n = 4, jitter = 2, seed = 77L))
  cat("  hmm_starts(): converged", msz$n_converged, " not",
      msz$n_not_converged, " error", msz$n_error, " modes",
      nrow(msz$modes), "\n")
  cat("  scale_how:", msz$scale_how, "\n")
}

hr("B. the replacement alignment check, on a REAL fit")
ci4 <- suppressWarnings(stats::confint(fit4))
vc4 <- suppressWarnings(stats::vcov(fit4))
cat("confint() rownames:", paste(rownames(ci4), collapse = " | "), "\n")
cat("vcov()    rownames:", paste(rownames(vc4), collapse = " | "), "\n")
shared <- intersect(rownames(vc4), rownames(ci4))
cat("shared rows        :", length(shared), "of", nrow(ci4), ":",
    paste(shared, collapse = " | "), "\n")
cat("rows with NO second path:",
    paste(setdiff(rownames(ci4), rownames(vc4)), collapse = " | "), "\n")
cat("baseline scale_how :", scalef(fit4)$how, "\n")

## a real fit whose confint() returns its standard errors permuted
## relative to its own rownames, which is exactly the failure the
## comment describes
mk_perm <- function(fit, i, j) {
  f <- fit
  attr(f, "rev_perm") <- c(i, j)
  class(f) <- c("rev_perm_fit", class(fit))
  f
}
registerS3method("confint", "rev_perm_fit", function(object, ...) {
  p <- attr(object, "rev_perm")
  o <- object
  class(o) <- setdiff(class(object), "rev_perm_fit")
  ci <- stats::confint(o, ...)
  ci[c(p[1L], p[2L]), c("lwr", "upr")] <-
    ci[c(p[2L], p[1L]), c("lwr", "upr")]
  ci
})
registerS3method("vcov", "rev_perm_fit", function(object, ...) {
  o <- object
  class(o) <- setdiff(class(object), "rev_perm_fit")
  stats::vcov(o, ...)
})
se4 <- (ci4[, "upr"] - ci4[, "lwr"]) / (2 * stats::qnorm(0.975))
nm4 <- rownames(ci4)
fx <- which(nm4 %in% rownames(vc4))
th <- which(!(nm4 %in% rownames(vc4)))
cat("\n-- B1: permute two rows that HAVE a second path --\n")
cat("   rows", fx[1L], "and", fx[2L], ":", nm4[fx[1L]], "and",
    nm4[fx[2L]], "\n")
cat("   their standard errors:", format(se4[fx[1:2]], digits = 4), "\n")
r1 <- scalef(mk_perm(fit4, fx[1L], fx[2L]))
cat("   scale_how:", r1$how,
    if (identical(r1$how, "unit")) " <- CAUGHT" else " <- MISSED", "\n")
cat("\n-- B2: permute two `theta_` rows, which have no second path --\n")
if (length(th) >= 2L) {
  cat("   rows", th[1L], "and", th[2L], ":", nm4[th[1L]], "and",
      nm4[th[2L]], "\n")
  cat("   their standard errors:", format(se4[th[1:2]], digits = 4), "\n")
  r2 <- scalef(mk_perm(fit4, th[1L], th[2L]))
  cat("   scale_how:", r2$how,
      if (identical(r2$how, "unit")) " <- caught after all" else
        " <- MISSED, as the comment says", "\n")
  cat("   the cost of missing it: the two are jittered by each other's\n")
  cat("   standard error, a factor of",
      format(max(se4[th[1:2]]) / min(se4[th[1:2]]), digits = 4), "\n")
} else {
  cat("   fewer than two rows without a second path; nothing to permute\n")
}
cat("\n-- B3: a permutation ACROSS the two groups --\n")
if (length(fx) && length(th)) {
  r3 <- scalef(mk_perm(fit4, fx[1L], th[1L]))
  cat("   rows", nm4[fx[1L]], "and", nm4[th[1L]], ": scale_how",
      r3$how, if (identical(r3$how, "unit")) " <- CAUGHT" else
        " <- MISSED", "\n")
}

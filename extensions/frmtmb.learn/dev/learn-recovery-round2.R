# Parameter recovery for the three families this round added or changed.
#
#   Rscript dev/learn-recovery-round2.R        # 60 replicates each
#   Rscript dev/learn-recovery-round2.R 10     # fewer, for a check
#
# WHAT IT ANSWERS, and it is not the same question for all three.
#
#   rlddm()      Does a joint choice-and-time family recover at a
#                realistic scale, and which of its parameters trade off?
#                A first-passage density sees largely the RATIO of the
#                drift to the boundary, so `drift` and `bs` are the pair
#                to watch.
#   igt_orl()    Is it weakly identified, and where? The first release
#                left it out on the grounds that it is. That is a claim
#                and this measures it rather than repeating it.
#   the Kalman   Does the exploration bonus `phi` separate from the
#   bonus        softmax sensitivity `tau`? Both control how far choice
#                departs from the current best arm.
#
# EVERYTHING IS REPORTED ON THE NATURAL SCALE. `rlddm()`'s non-decision
# time forces it: its link is a logit scaled onto (0, min(rt)), so the
# BOUND differs from one simulated dataset to the next and two
# replicates' link-scale intercepts are not the same quantity. The Wald
# interval endpoints are pushed through the same monotone link, which
# leaves the coverage a coverage.
#
# THE CORRELATION MATRIX IS THE POINT for the last two. A bias table
# says whether an estimator lands on the truth; it does not say whether
# two parameters are separately identified. The correlation of their
# estimates ACROSS replicates does: a pair that trades off comes back
# with a correlation near -1 or +1 whatever their individual biases are.

suppressPackageStartupMessages(library(frmtmb.learn))

args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args)) as.integer(args[[1L]]) else 60L

# A monotone map from the link scale frmtmb estimates on to the scale
# the model is stated on, applied to a point and to both endpoints of
# its interval.
NAT <- list(logit = stats::plogis, log = exp, identity = identity)

# One replicate of one design. `sim` draws the data, `form` is the
# model, `report` names the fixed effects to keep with the map each one
# needs. Returns the estimates and their intervals on the natural scale.
one_rep <- function(seed, sim, fam, form, report) {
  d <- sim(seed)
  fit <- try(frmtmb::frm(form, family = fam, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  ci <- try(stats::confint(fit), silent = TRUE)
  if (inherits(ci, "try-error")) return(NULL)
  b <- unlist(frmtmb::fixef(fit))
  rn <- rownames(ci)
  est <- lo <- hi <- rep(NA_real_, length(report))
  for (j in seq_along(report)) {
    nm <- names(report)[[j]]
    f <- report[[j]](d)
    # match(), not grep(): a fixed-effect name carries literal
    # parentheses and "alpha.(Intercept)" read as a regex matches
    # nothing at all, silently, which is a whole table of NaN
    k <- match(nm, names(b))
    if (is.na(k)) next
    est[[j]] <- f(b[[k]])
    m <- grep(sub("^([^.]+)[.]", "\\1_", nm), rn, fixed = TRUE)
    if (length(m)) {
      # the link is monotone increasing in every case here, so the
      # interval maps endpoint to endpoint
      lo[[j]] <- f(ci[m[[1L]], 1L])
      hi[[j]] <- f(ci[m[[1L]], 2L])
    }
  }
  list(est = est, lo = lo, hi = hi)
}

summarise <- function(reps, truth, labels) {
  reps <- Filter(Negate(is.null), reps)
  est <- do.call(rbind, lapply(reps, function(r) r$est))
  lo <- do.call(rbind, lapply(reps, function(r) r$lo))
  hi <- do.call(rbind, lapply(reps, function(r) r$hi))
  n <- nrow(est)
  cov <- vapply(seq_along(truth), function(j) {
    ok <- is.finite(lo[, j]) & is.finite(hi[, j])
    if (!any(ok)) return(NA_real_)
    mean(lo[ok, j] <= truth[[j]] & hi[ok, j] >= truth[[j]])
  }, 0)
  out <- data.frame(
    parameter = labels,
    truth = round(unlist(truth), 3),
    bias = round(colMeans(est, na.rm = TRUE) - unlist(truth), 3),
    mc_se = round(apply(est, 2, stats::sd, na.rm = TRUE) / sqrt(n), 3),
    sd_of_est = round(apply(est, 2, stats::sd, na.rm = TRUE), 3),
    coverage = round(cov, 2),
    n_fit = n)
  list(table = out, cor = round(stats::cor(est, use = "pairwise"), 2))
}

run <- function(title, sim, fam, form, report, truth, labels) {
  cat("\n== ", title, ", ", NREP, " replicates ==\n", sep = "")
  t0 <- Sys.time()
  reps <- lapply(seq_len(NREP), function(i)
    one_rep(2000L + i, sim, fam, form, report))
  s <- summarise(reps, truth, labels)
  print(s$table, row.names = FALSE)
  cat("\nestimate correlations across replicates\n")
  dimnames(s$cor) <- list(labels, labels)
  print(s$cor)
  cat(sprintf("(%.1f min)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(s)
}

## ---------------------------------------------------------------- rlddm
rlddm_fam <- rlddm(subject = id, trial = trial)
rlddm_sim <- function(seed) {
  set.seed(seed)
  d <- frm_task_design("bandit2arm", n_subject = 30L, n_trial = 100L,
                       seed = seed)
  u <- stats::rnorm(30L, 0, 0.4)
  al <- stats::plogis(stats::qlogis(0.35) + u[as.integer(d$id)])
  frm_task_simulate(rlddm_fam, d,
                    pars = list(alpha = al, drift = 3, bs = 1.6,
                                ndt = 0.2, bias = 0.5),
                    seed = seed)[[1L]]
}
run("rlddm, 30 subjects x 100 trials", rlddm_sim, rlddm_fam,
    frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
               drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
    list("alpha.(Intercept)" = function(d) NAT$logit,
         "drift.(Intercept)" = function(d) NAT$identity,
         "bs.(Intercept)" = function(d) NAT$log,
         # the bound is this dataset's fastest response, so the map is
         # rebuilt per replicate rather than fixed once
         "ndt.(Intercept)" = function(d) {
           ub <- min(d$rt)
           function(x) ub / (1 + exp(-x))
         },
         "bias.(Intercept)" = function(d) NAT$logit),
    truth = list(0.35, 3, 1.6, 0.2, 0.5),
    labels = c("alpha", "drift", "bs", "ndt", "bias"))

## -------------------------------------------------------------- igt_orl
orl_fam <- igt_orl(subject = id, trial = trial)
orl_sim <- function(seed) {
  set.seed(seed)
  d <- frm_task_design("igt", n_subject = 30L, n_trial = 100L, seed = seed)
  u <- stats::rnorm(30L, 0, 0.4)
  ar <- stats::plogis(stats::qlogis(0.3) + u[as.integer(d$id)])
  d$choice <- frm_task_simulate(
    orl_fam, d, pars = list(Arew = ar, Apun = 0.1, k = 0.5, betaF = 1,
                            betaP = 1), seed = seed)[[1L]]$choice
  d
}
run("igt_orl, 30 subjects x 100 trials", orl_sim, orl_fam,
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1 + (1 | id),
               Apun ~ 1, k ~ 1, betaF ~ 1, betaP ~ 1),
    list("Arew.(Intercept)" = function(d) NAT$logit,
         "Apun.(Intercept)" = function(d) NAT$logit,
         "k.(Intercept)" = function(d) NAT$log,
         "betaF.(Intercept)" = function(d) NAT$identity,
         "betaP.(Intercept)" = function(d) NAT$identity),
    truth = list(0.3, 0.1, 0.5, 1, 1),
    labels = c("Arew", "Apun", "k", "betaF", "betaP"))

## ------------------------------------------------- the exploration bonus
kal_fam <- bandit4arm2_kalman_filter(subject = id, trial = trial,
                                     bonus = TRUE)
kal_sim <- function(seed) {
  set.seed(seed)
  d <- frm_task_design("bandit4arm_restless", n_subject = 30L,
                       n_trial = 100L, seed = seed)
  u <- stats::rnorm(30L, 0, 0.3)
  tau <- exp(log(0.15) + u[as.integer(d$id)])
  d$choice <- frm_task_simulate(
    kal_fam, d, pars = list(tau = tau, lambda = 0.98, center = 50,
                            mu0 = 50, sigma0 = 10, sigmaD = 3,
                            phi = 1.5), seed = seed)[[1L]]$choice
  d
}
# center, mu0 and sigma0 are HELD at the task's own values through
# bf(name = value), not estimated. The family's help already says the
# three do not separate from the rest on a session this size, and this
# study asks about `phi`; leaving them free would measure their
# collapse a second time and put its noise into phi's column.
run("bandit4arm2_kalman_filter(bonus = TRUE), 30 x 100", kal_sim, kal_fam,
    frmtmb::bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1 + (1 | id),
               lambda ~ 1, sigmaD ~ 1, phi ~ 1,
               center = 50, mu0 = 50, sigma0 = 10),
    list("tau.(Intercept)" = function(d) NAT$log,
         "lambda.(Intercept)" = function(d) NAT$logit,
         "sigmaD.(Intercept)" = function(d) NAT$log,
         "phi.(Intercept)" = function(d) NAT$identity),
    truth = list(0.15, 0.98, 3, 1.5),
    labels = c("tau", "lambda", "sigmaD", "phi"))

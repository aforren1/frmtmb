# Parameter recovery, and the Laplace caveat measured.
#
#   Rscript dev/learn-recovery.R              # every family, 100 reps
#   Rscript dev/learn-recovery.R 40           # fewer reps, for a check
#
# WHAT IT ANSWERS. Two questions that a fit converging does not answer.
# Do the estimates land on the truth (bias, against its Monte Carlo
# standard error), and do their 95 percent Wald intervals cover at the
# nominal rate?
#
# WHY IT ALSO MEASURES THE LAPLACE CAVEAT. frmtmb integrates the subject
# effects out with a Laplace approximation, exact only when the
# conditional log-density is quadratic. For binary choices it is not, and
# the fewer trials a subject has the less quadratic it is. The usual way
# to price that error is frm(importance =), which every family here
# REFUSES, so the consequence is measured by running the same study at
# two session lengths instead.
#
# READ THE VARIANCE COMPONENTS CAREFULLY. A variance component estimated
# by maximum likelihood from binary data with few levels is biased
# downward whether or not the integral is approximated, and this study
# does not separate the two causes. It establishes where the answers are
# safe, not which of the two is responsible.

library(frmtmb.learn)

args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args)) as.integer(args[[1]]) else 100L

# One replicate: draw subject parameters from the truth, draw choices
# from the family's own generative simulator, fit, and return the
# estimates with their standard errors.
one_rep <- function(seed, ns, nt, truth, sd_u) {
  set.seed(seed)
  d <- frm_task_design("reversal", n_subject = ns, n_trial = nt,
                       seed = seed)
  cb <- as.numeric(d$after_reversal == "after")
  u <- stats::rnorm(ns, 0, sd_u)
  eta <- truth[["alpha_Intercept"]] + truth[["alpha_after"]] * cb +
    u[as.integer(d$id)]
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(eta), tau = exp(truth[["tau"]])),
    seed = seed)[[1L]]$choice
  fit <- try(frm(bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
                    tau ~ 1),
                 family = bandit2arm_delta(subject = id, trial = trial),
                 data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  ci <- try(stats::confint(fit), silent = TRUE)
  if (inherits(ci, "try-error")) return(NULL)
  b <- unlist(fixef(fit))
  sdhat <- sqrt(VarCorr(fit)[[1L]][1L, 1L])
  keep <- c("alpha.(Intercept)", "alpha.after_reversalafter",
            "tau.(Intercept)")
  rn <- rownames(ci)
  pick <- function(pat) {
    j <- grep(pat, rn, fixed = TRUE)
    if (!length(j)) c(NA_real_, NA_real_) else as.numeric(ci[j[1L], 1:2])
  }
  # theta_1 is the log standard deviation of the (1 | id) block on the
  # scale frmtmb estimates it on, which is the same quantity as
  # log(sqrt(VarCorr)), so its interval is a Wald interval for the
  # variance component and its coverage is comparable with the others.
  iv <- rbind(pick("alpha_(Intercept)"), pick("after_reversalafter"),
              pick("tau_(Intercept)"), pick("theta_1"))
  list(est = c(alpha_Intercept = unname(b[keep[1L]]),
               alpha_after = unname(b[keep[2L]]),
               tau = unname(b[keep[3L]]),
               log_sd = log(sdhat)),
       lo = iv[, 1L], hi = iv[, 2L])
}

summarise <- function(reps, truth, sd_u) {
  reps <- Filter(Negate(is.null), reps)
  est <- do.call(rbind, lapply(reps, function(r) r$est))
  lo <- do.call(rbind, lapply(reps, function(r) r$lo))
  hi <- do.call(rbind, lapply(reps, function(r) r$hi))
  tv <- c(truth[["alpha_Intercept"]], truth[["alpha_after"]],
          truth[["tau"]], log(sd_u))
  n <- nrow(est)
  cov <- rep(NA_real_, 4)
  for (j in 1:4) {
    ok <- is.finite(lo[, j]) & is.finite(hi[, j])
    cov[j] <- mean(lo[ok, j] <= tv[j] & hi[ok, j] >= tv[j])
  }
  data.frame(
    parameter = c("alpha_(Intercept)", "alpha_after", "tau_(Intercept)",
                  "log sd(alpha)"),
    truth = round(tv, 3),
    bias = round(colMeans(est, na.rm = TRUE) - tv, 3),
    mc_se = round(apply(est, 2, stats::sd, na.rm = TRUE) / sqrt(n), 3),
    sd_of_est = round(apply(est, 2, stats::sd, na.rm = TRUE), 3),
    coverage = round(cov, 2),
    n_int = vapply(1:4, function(j)
      sum(is.finite(lo[, j]) & is.finite(hi[, j])), numeric(1)),
    n_fit = n)
}

TRUTH <- list(alpha_Intercept = stats::qlogis(0.30), alpha_after = 1.0,
              tau = log(3))
SD_U <- 0.5

for (nt in c(100L, 20L)) {
  cat("
== 40 subjects,", nt, "trials each,", NREP, "replicates ==
")
  t0 <- Sys.time()
  reps <- lapply(seq_len(NREP), function(i)
    one_rep(1000L + i, ns = 40L, nt = nt, truth = TRUTH, sd_u = SD_U))
  print(summarise(reps, TRUTH, SD_U), row.names = FALSE)
  cat(sprintf("(%.1f min)
",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("
The fixed effects are the ones to report from a short session.",
    "
A subject-level standard deviation from twenty binary trials is",
    "a
lower bound, and frm(importance =) is what would separate",
    "Laplace
error from the ordinary downward bias; it is refused.
")

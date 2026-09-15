# Lane `learnhier` (item 2.2): ONE replicate of one design, into one
# RDS.
#
#   Rscript dev/learnhier-run.R <design> <seed> <outdir> [ns] [nt]
#
# `design` is `bandit`, `bandit0` (the same design with a diagonal truth)
# or `rlddm`. One replicate per PROCESS, because a 100 by 200 rlddm fit
# peaks near 8 GB of working set and a process that has already built one
# tape is not the process to time or to measure a second in. The
# summariser reads the directory.
#
# Nothing here is asserted. The run records; dev/learnhier-summarize.R
# decides.

source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: <design> <seed> <outdir> [ns] [nt]")
DESIGN <- args[[1L]]
SEED <- as.integer(args[[2L]])
OUTDIR <- args[[3L]]
NS <- if (length(args) >= 4L) as.integer(args[[4L]]) else 100L
NT <- if (length(args) >= 5L) as.integer(args[[5L]]) else 200L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Idempotent, so that a launcher can be re-run after an interruption
# without paying for the replicates it already has. A 100 by 200 rlddm
# fit is fifteen minutes and the block of sixty is the largest single
# cost in this lane.
OUTFILE <- file.path(OUTDIR, sprintf("%s-%d.rds", DESIGN, SEED))
if (file.exists(OUTFILE)) {
  cat("HAVE ", OUTFILE, "\n", sep = "")
  quit(save = "no", status = 0L)
}

# ------------------------------------------------------- small helpers

# "rt.alpha:(Intercept)" is one row of a block; the block's own name for
# it is "alpha". confint_varcorr() reports a correlation as
# "cor(<row>,<row>)", so both are parsed here rather than matched by
# position: the order of the lower triangle is a convention that a
# reader of this table should not have to know.
lh_short <- function(x) {
  sub("^[^.]+[.]", "", sub(":[(]Intercept[)]$", "", x))
}

lh_vc_flat <- function(vc) {
  out <- numeric(0)
  for (k in seq_len(nrow(vc))) {
    tm <- vc$term[[k]]
    if (identical(vc$type[[k]], "sd")) {
      nm <- paste0("sd_", lh_short(tm))
    } else if (identical(vc$type[[k]], "cor")) {
      inner <- sub("^cor[(](.*)[)]$", "\\1", tm)
      ab <- strsplit(inner, ",", fixed = TRUE)[[1L]]
      nm <- paste0("cor_", lh_short(ab[[1L]]), "~", lh_short(ab[[2L]]))
    } else {
      next
    }
    out[paste0(nm, ".est")] <- vc$estimate[[k]]
    out[paste0(nm, ".lwr")] <- vc$lwr[[k]]
    out[paste0(nm, ".upr")] <- vc$upr[[k]]
  }
  out
}

lh_fix_flat <- function(fit) {
  b <- unlist(frmtmb::fixef(fit))
  ci <- suppressWarnings(stats::confint(fit))
  out <- numeric(0)
  for (nm in names(b)) {
    key <- sub("^([^.]+)[.]", "\\1_", nm)
    k <- match(key, rownames(ci))
    out[paste0("b_", nm, ".est")] <- b[[nm]]
    out[paste0("b_", nm, ".lwr")] <- if (is.na(k)) NA_real_ else ci[k, 1L]
    out[paste0("b_", nm, ".upr")] <- if (is.na(k)) NA_real_ else ci[k, 2L]
  }
  out
}

# The ORACLE log-likelihood: this replicate's objective evaluated at the
# parameters that generated it, with the covariance block set to the
# REALIZED moments of the deviations the draw actually produced.
#
# Why it is worth a few seconds. The fit's own optimum must be at or
# above it, so a fit that comes back BELOW the oracle has stopped
# somewhere that is not the optimum, whatever its gradient and Hessian
# say. That is precisely the failure item 2.3 found in the hidden Markov
# row: convergence 0, a positive definite Hessian, and 8.099 units left
# on the table. It is a check on the optimizer, not on the estimator,
# and it is reported as a difference so it needs no tolerance.
lh_oracle <- function(fit, beta, betad, D) {
  reg <- frmtmb:::covstruct_registry[["us"]]
  blk <- fit$frame[["re_blocks"]][[1L]]
  th <- reg$from_natural(apply(D, 2L, stats::sd), stats::cor(D), blk)
  par <- fit$opt$par
  par[names(par) == "beta"] <- beta
  par[names(par) == "betad"] <- betad
  par[names(par) == "theta"] <- th
  -as.numeric(fit$obj$fn(par))
}

# ------------------------------------------------------------ designs

t0 <- Sys.time()
if (DESIGN %in% c("bandit", "bandit0")) {
  truth <- if (DESIGN == "bandit") lh_bandit_truth else lh_bandit_truth0
  d <- lh_bandit_data(SEED, ns = NS, nt = NT, truth = truth)
  form <- frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | p | id),
                     tau ~ 1 + (1 | p | id))
  fam <- bandit2arm_delta(subject = id, trial = trial)
  oracle_beta <- stats::qlogis(truth$alpha)
  oracle_betad <- log(truth$tau)
} else if (DESIGN %in% c("rlddm", "rlddm5")) {
  truth <- lh_rlddm_truth
  d <- lh_rlddm_data(SEED, ns = NS, nt = NT, truth = truth)
  # rlddm5 is the SAME data with `bias` in the block as well, which is
  # the reading of "every parameter" the scale row declined. Its truth
  # for that component is a variance of exactly zero, so it is the
  # boundary case, and it is a probe on a couple of seeds rather than an
  # arm: see dev/learnhier-findings.md.
  form <- if (DESIGN == "rlddm") {
    frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) +
                 ndt_group(id) ~ 1 + (1 | p | id),
               drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
               ndt ~ 1 + (1 | p | id), bias = 0.5)
  } else {
    frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) +
                 ndt_group(id) ~ 1 + (1 | p | id),
               drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
               ndt ~ 1 + (1 | p | id), bias ~ 1 + (1 | p | id))
  }
  fam <- rlddm(subject = id, trial = trial)
  oracle_beta <- stats::qlogis(truth$alpha)
  # `ndt`'s intercept has no population truth under ndt_group(): the
  # link's bound is each learner's own observed floor. The realized mean
  # of the fitted-scale deviations is the oracle's value for it, which
  # is why the oracle is stated as realized rather than as population.
  oracle_betad <- c(truth$drift, log(truth$bs),
                    mean(attr(d, "dev_fitted")[, "ndt"]))
  # the five-parameter block's fifth component has a realized standard
  # deviation of exactly zero, whose log the oracle's theta would need,
  # so there is no oracle to state for it
  if (DESIGN == "rlddm5") oracle_betad <- NULL
} else {
  stop("unknown design: ", DESIGN)
}
t_sim <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

rec <- list(design = DESIGN, seed = SEED, ns = NS, nt = NT,
            rows = nrow(d), sim_s = t_sim,
            truth_drawn = lh_block_stats(attr(d, "dev_drawn")),
            truth_fitted = lh_block_stats(attr(d, "dev_fitted")),
            # what the oracle used, so that the summariser can score the
            # `ndt` intercept, which has no population truth under
            # ndt_group(), against the value the draw realized
            oracle_beta = oracle_beta, oracle_betad = oracle_betad)

t0 <- Sys.time()
fit <- try(frmtmb::frm(form, family = fam, data = d, se = TRUE),
           silent = TRUE)
rec$fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rec$ok <- !inherits(fit, "try-error")
if (!rec$ok) {
  rec$why <- conditionMessage(attr(fit, "condition"))
} else {
  dg <- frmtmb::diagnose(fit, quiet = TRUE)
  rec$conv <- fit$opt$convergence
  rec$max_grad <- dg$max_grad
  rec$pdHess <- isTRUE(dg$pdHess)
  rec$n_bad_se <- sum(!is.finite(fit$sdr$sd))
  rec$logLik <- as.numeric(stats::logLik(fit))
  rec$n_par <- length(fit$opt$par)
  rec$fixed <- lh_fix_flat(fit)
  rec$vc <- lh_vc_flat(suppressWarnings(frmtmb::confint_varcorr(fit)))
  if (DESIGN %in% c("rlddm", "rlddm5")) {
    one <- d[match(levels(d$id), as.character(d$id)), , drop = FALSE]
    hat <- as.numeric(suppressWarnings(
      frmtmb.eam::ndt_time(fit, newdata = one)))
    tru <- as.numeric(attr(d, "ndt_subject"))
    own <- as.numeric(attr(d, "own_floor"))
    nd <- suppressWarnings(stats::predict(
      fit, newdata = one[1L, , drop = FALSE], dpar = "ndt",
      type = "response", re.form = NA, se.fit = TRUE))
    frac <- as.numeric(nd$fit[1L])
    rec$ndt <- c(
      rmse_ms = 1000 * sqrt(mean((hat - tru)^2)),
      rmse_ratio = sqrt(mean((hat - tru)^2)) / stats::sd(tru),
      mean_err_ratio = abs(mean(hat) - mean(tru)) / stats::sd(tru),
      cor_true = stats::cor(hat, tru),
      # the two comparisons that need no second fit: the best CONSTANT
      # fraction of each learner's own floor, which is what an estimator
      # with the component switched off can do
      floor_rmse_ms = 1000 * sqrt(mean((mean(tru / own) * own - tru)^2)),
      floor_cor_true = stats::cor(own, tru),
      sd_hat = stats::sd(hat), sd_true = stats::sd(tru),
      sd_floor_only = stats::sd(frac * own),
      frac_pop = frac, floor_mean = mean(own),
      margin_min_ms = 1000 * min(own - hat),
      below_own_floor = sum(own - hat > 0), n_subject = length(hat))
  }
  # LAST, because it moves the objective's last.par
  rec$oracle <- if (is.null(oracle_betad)) NA_real_ else {
    try(lh_oracle(fit, oracle_beta, oracle_betad,
                  attr(d, "dev_fitted")), silent = TRUE)
  }
  if (inherits(rec$oracle, "try-error")) rec$oracle <- NA_real_
}
rec$total_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rec$peak_r_mb <- sum(gc()[, 6L])

f <- OUTFILE
saveRDS(rec, f)
cat("WROTE ", f, " ok=", rec$ok, " fit_s=", round(rec$fit_s, 1),
    if (isTRUE(rec$ok)) paste0(" conv=", rec$conv, " ll=",
                               round(rec$logLik, 2),
                               " dll_oracle=",
                               round(rec$logLik - rec$oracle, 3)) else "",
    "\n", sep = "")

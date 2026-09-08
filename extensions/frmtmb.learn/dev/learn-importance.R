# What the importance correction does to a learning fit.
#
#   Rscript dev/learn-importance.R           # 6 seeds per cell
#   Rscript dev/learn-importance.R 3         # fewer, for a check
#
# WHAT IT ANSWERS. frmtmb integrates the subject effects out with a
# Laplace approximation, exact only when the conditional log-density is
# quadratic. For binary choices it is not. `frm(importance =)` reweights
# draws from the Laplace Gaussian and returns the corrected maximum
# likelihood estimate, so the DIFFERENCE between the two fits is the
# Laplace error, measured rather than argued about.
#
# The families could not answer it until they declared how their
# likelihood factorizes; see the two `loglik_row` / `loglik_group` slots
# in R/family.R. This script is what turns "the correction now runs"
# into a number with an error bar on it.
#
# READ THE DIAGNOSTICS, NOT ONLY THE SHIFT. A corrected estimate is a
# weighted average over draws, so it carries its own Monte Carlo error.
# Two columns say how far to trust it: `mcse`, the Monte Carlo standard
# error of the corrected log-likelihood, and the smallest effective
# sample size per draw over the grouping levels. An ESS per draw near 1
# means the weights inside every group are nearly uniform and the
# correction is reading its own integrand well; near zero means one draw
# in a group carries all the weight and the corrected number is that one
# draw.

suppressPackageStartupMessages(library(frmtmb.learn))

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1]]) else 6L

TRUTH <- c(alpha_Intercept = stats::qlogis(0.3), alpha_after = 1,
           tau = log(3), log_sd = log(0.5))
NS <- 40L

# One dataset from the family's own generative simulator, with a
# per-subject offset on the learning rate's logit and a condition effect
# that turns on at the reversal.
make_data <- function(seed, nt, sd_u = 0.5) {
  set.seed(seed)
  d <- frm_task_design("reversal", n_subject = NS, n_trial = nt, seed = seed)
  cb <- as.numeric(d$after_reversal == "after")
  u <- stats::rnorm(NS, 0, sd_u)
  eta <- TRUTH[["alpha_Intercept"]] + TRUTH[["alpha_after"]] * cb +
    u[as.integer(d$id)]
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(eta), tau = exp(TRUTH[["tau"]])),
    seed = seed)[[1L]]$choice
  d
}

FORM <- frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal + (1 | id),
                   tau ~ 1)

pull <- function(fit) {
  b <- unlist(frmtmb::fixef(fit))
  c(b[["alpha.(Intercept)"]], b[["alpha.after_reversalafter"]],
    b[["tau.(Intercept)"]],
    log(sqrt(frmtmb::VarCorr(fit)[[1L]][1L, 1L])))
}

PARS <- c("alpha_(Intercept)", "alpha_after", "tau_(Intercept)",
          "log sd(alpha)")

one <- function(seed, nt, ndraw) {
  d <- make_data(seed, nt)
  fam <- bandit2arm_delta(subject = id, trial = trial)
  t0 <- proc.time()[[3L]]
  lap <- try(frmtmb::frm(FORM, family = fam, data = d), silent = TRUE)
  if (inherits(lap, "try-error")) return(NULL)
  t_lap <- proc.time()[[3L]] - t0
  t0 <- proc.time()[[3L]]
  imp <- try(frmtmb::frm(FORM, family = fam, data = d, importance = ndraw),
             silent = TRUE)
  t_imp <- proc.time()[[3L]] - t0
  if (inherits(imp, "try-error")) {
    return(list(seed = seed, nt = nt, ndraw = ndraw, ok = FALSE,
                lap = pull(lap), t_lap = t_lap, t_imp = t_imp,
                why = conditionMessage(attr(imp, "condition"))))
  }
  ii <- imp[["importance"]]
  list(seed = seed, nt = nt, ndraw = ndraw, ok = TRUE, lap = pull(lap),
       cor = pull(imp), mcse = ii[["mcse"]], ess = ii[["ess_min"]],
       t_lap = t_lap, t_imp = t_imp)
}

# PER SEED, not averaged. A variance component estimated from short
# binary sessions collapses on some datasets and not others, and one
# collapsed replicate moves a mean of six by more than the whole effect
# being measured. The fixed effects are stable enough to average; the
# variance component is reported one dataset at a time, with the
# diagnostics beside it, because the question "is the correction usable
# here" is a question about a dataset.
report <- function(rows, nt, ndraw) {
  cat("
## ", NS, " subjects x ", nt, " trials, importance = ", ndraw,
      "

", sep = "")
  bad <- Filter(function(r) !isTRUE(r[["ok"]]), rows)
  ok <- Filter(function(r) isTRUE(r[["ok"]]), rows)
  cat("completed ", length(ok), " of ", length(rows), "

", sep = "")
  for (r in bad) {
    cat("  seed ", r[["seed"]], " REFUSED: ",
        substr(gsub("[[:space:]]+", " ", r[["why"]]), 1L, 220L), "
",
        sep = "")
  }
  if (!length(ok)) return(invisible(NULL))
  cat(sprintf("%6s %10s %10s %8s %8s %9s %9s
", "seed", "sd Laplace",
              "sd corrtd", "log shft", "mcse", "minESS/dr", "seconds"))
  for (r in ok) {
    cat(sprintf("%6d %10.4f %10.4f %+8.3f %8.3f %9.3f %9.1f
",
                r[["seed"]], exp(r[["lap"]][[4L]]), exp(r[["cor"]][[4L]]),
                r[["cor"]][[4L]] - r[["lap"]][[4L]], r[["mcse"]],
                r[["ess"]], r[["t_imp"]]))
  }
  sh <- vapply(ok, function(r) r[["cor"]] - r[["lap"]], numeric(4L))
  la <- vapply(ok, function(r) r[["lap"]], numeric(4L))
  cat("
fixed effects, averaged over the completed replicates
")
  cat(sprintf("%-20s %9s %9s %9s %9s
", "parameter", "truth",
              "Laplace", "shift", "largest"))
  for (k in seq_len(3L)) {
    cat(sprintf("%-20s %9.3f %9.3f %+9.3f %+9.3f
", PARS[[k]],
                TRUTH[[k]], mean(la[k, ]), mean(sh[k, ]),
                sh[k, which.max(abs(sh[k, ]))]))
  }
  cat(sprintf("Laplace fit seconds, mean %.1f
",
              mean(vapply(ok, function(r) r[["t_lap"]], 0))))
  invisible(NULL)
}

seeds <- 4300L + seq_len(NSEED)
for (cell in list(c(nt = 100L, nd = 100L), c(nt = 20L, nd = 100L),
                  c(nt = 20L, nd = 400L))) {
  rows <- lapply(seeds, one, nt = cell[["nt"]], ndraw = cell[["nd"]])
  report(Filter(Negate(is.null), rows), cell[["nt"]], cell[["nd"]])
  utils::flush.console()
}

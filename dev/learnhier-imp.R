# Lane `learnhier` (item 2.2): the importance correction at the plan's
# realistic scale, counted with 0.55.x's `imp_stalled()`.
#
#   Rscript dev/learnhier-imp.R <design> <seed> <draws> <outdir>
#
# WHAT IS BEING COUNTED. `frm(importance = n)` reweights draws from the
# Laplace proposal and re-optimizes, repeating until the parameters stop
# moving or a round cap bites. Three outcomes, and only the first is a
# success:
#
#   settled   the iteration reached its own fixed point inside the cap.
#   stalled   the cap bit AND the rounds are taking the same step every
#             time, which is what a collapsed variance component does to
#             the fixed-point map. More rounds do not help.
#   slow      the cap bit and the steps are still shrinking. More rounds
#             would help.
#
# `frmtmb:::imp_stalled()` is what separates the last two, and it is the
# statistic the plan names. It is reached with a colon here because it
# is `@noRd` in core: this is a lane script and not package code, and
# `dev/learnhier-findings.md` records that as the reason the shipped
# tier reads `fit$importance$moves` and applies the rule itself.
#
# THE THRESHOLD IS TIGHT. `dev/reviews/2026-09-08-debts.md` measured the
# gap the 1e-2 threshold sits in at a factor of 3.2, not the three
# orders of magnitude the roxygen claimed: the widest stalled spread was
# 6.51e-03 and the narrowest still-moving one 2.10e-02. Every run below
# therefore records its own spread, so a rate near either edge is
# reported as a finding about the threshold rather than about the model.

source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: <design> <seed> <draws> <outdir>")
DESIGN <- args[[1L]]
SEED <- as.integer(args[[2L]])
DRAWS <- as.integer(args[[3L]])
OUTDIR <- args[[4L]]
NS <- if (length(args) >= 5L) as.integer(args[[5L]]) else 100L
NT <- if (length(args) >= 6L) as.integer(args[[6L]]) else 200L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Idempotent: a correction at this scale is about twelve minutes and a
# launcher re-run should not pay for it twice.
OUTFILE <- file.path(OUTDIR,
                     sprintf("imp-%s-%d-%d.rds", DESIGN, DRAWS, SEED))
if (file.exists(OUTFILE)) {
  cat("HAVE ", OUTFILE, "\n", sep = "")
  quit(save = "no", status = 0L)
}

truth <- if (DESIGN == "bandit0") lh_bandit_truth0 else lh_bandit_truth
d <- lh_bandit_data(SEED, ns = NS, nt = NT, truth = truth)
form <- frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | p | id),
                   tau ~ 1 + (1 | p | id))
fam <- bandit2arm_delta(subject = id, trial = trial)

pull <- function(fit) {
  vc <- frmtmb::VarCorr(fit)[[1L]]
  c(unlist(frmtmb::fixef(fit)),
    sd_alpha = sqrt(vc[1L, 1L]), sd_tau = sqrt(vc[2L, 2L]),
    cor = stats::cov2cor(vc)[1L, 2L])
}

rec <- list(design = DESIGN, seed = SEED, draws = DRAWS, ns = NS, nt = NT)
t0 <- Sys.time()
lap <- try(frmtmb::frm(form, family = fam, data = d), silent = TRUE)
rec$lap_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rec$lap_ok <- !inherits(lap, "try-error")
if (rec$lap_ok) {
  rec$lap <- pull(lap)
  rec$lap_ll <- as.numeric(stats::logLik(lap))
}

warns <- character(0)
t0 <- Sys.time()
imp <- withCallingHandlers(
  try(frmtmb::frm(form, family = fam, data = d, importance = DRAWS),
      silent = TRUE),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
rec$imp_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rec$warns <- warns
rec$imp_ok <- !inherits(imp, "try-error")
if (!rec$imp_ok) {
  rec$why <- conditionMessage(attr(imp, "condition"))
  rec$outcome <- "refused"
} else {
  ii <- imp[["importance"]]
  mv <- ii[["moves"]]
  rec$moves <- mv
  rec$rounds <- ii[["rounds"]]
  rec$capped <- isTRUE(ii[["capped"]])
  rec$moved <- ii[["moved"]]
  rec$spread <- if (length(mv) >= 2L) {
    (max(mv) - min(mv)) / mean(mv)
  } else NA_real_
  rec$stalled <- frmtmb:::imp_stalled(mv)
  rec$ess_min <- ii[["ess_min"]]
  rec$ess_median <- ii[["ess_median"]]
  rec$mcse <- ii[["mcse"]]
  rec$grad <- ii[["grad"]]
  rec$imp <- pull(imp)
  rec$imp_ll <- as.numeric(stats::logLik(imp))
  rec$outcome <- if (!rec$capped) {
    "settled"
  } else if (isTRUE(rec$stalled)) "stalled" else "slow"
}

f <- OUTFILE
saveRDS(rec, f)
cat("WROTE ", f, " outcome=", rec$outcome %||% "NA",
    " rounds=", rec$rounds %||% NA, " spread=",
    format(rec$spread %||% NA, digits = 3),
    " imp_s=", round(rec$imp_s, 1), "\n", sep = "")

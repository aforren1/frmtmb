# REVIEW, item 1.0b, attack 1: is the per-learner `ndt` recovery the
# model's, or the observed floors'?
#
# Usage (from the worktree root):
#   Rscript dev/rev-rlddm-floors.R <arm> <out.tsv>
#
# arm = "none"  no fit at all, just the floors-only arithmetic
# arm = "re"    the branch's own scale row: ndt ~ 1 + (1 | p | id)
# arm = "nore"  the same design with the random effect on ndt REMOVED,
#               so the fitted per-learner ndt is one constant fraction
#               times each learner's own fastest response
#
# Seed 20260908, the tier's own, and the data function is copied from
# extensions/frmtmb.learn/tests/testthat/test-scale.R at scale_small()
# FALSE, 100 learners by 200 trials.

args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1L]]
out <- if (length(args) > 1L) args[[2L]] else ""

.libPaths(c("C:/Users/adf44/source/r/rev-rlddm-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("arm  :", arm, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")
cat("eam  :", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n")
cat("frmtmb:", format(packageVersion("frmtmb")), "at",
    dirname(system.file(package = "frmtmb")), "\n")

truth <- list(alpha = 0.35, drift = 2.5, bs = 1.5, ndt = 0.25,
              sd_alpha = 0.5, sd_drift = 1.0, sd_bs = 0.2,
              sd_ndt = 0.15)

make_data <- function(seed = 20260908L, ns = 100L, nt = 200L) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                       seed = seed)
  set.seed(seed + 2L)
  i <- as.integer(d$id)
  ua <- stats::rnorm(ns, 0, truth$sd_alpha)
  ud <- stats::rnorm(ns, 0, truth$sd_drift)
  ub <- stats::rnorm(ns, 0, truth$sd_bs)
  un <- stats::rnorm(ns, 0, truth$sd_ndt)
  s <- frm_task_simulate(
    rlddm(subject = id, trial = trial), d,
    pars = list(alpha = stats::plogis(stats::qlogis(truth$alpha) + ua[i]),
                drift = truth$drift + ud[i],
                bs = truth$bs * exp(ub[i]),
                ndt = truth$ndt * exp(un[i]),
                bias = 0.5),
    seed = seed)[[1L]]
  attr(s, "ndt_subject") <- truth$ndt * exp(un)
  s
}

d <- make_data()
lv <- levels(d$id)
one <- d[match(lv, as.character(d$id)), , drop = FALSE]
own <- as.numeric(tapply(d$rt, d$id, min))[match(as.character(one$id), lv)]
tru <- as.numeric(attr(d, "ndt_subject"))[match(as.character(one$id), lv)]

rmse <- function(a, b) sqrt(mean((a - b)^2))
say <- function(lab, hat) {
  cat(sprintf("%-34s rmse_ms=%9.4f  cor=%.6f  ratio=%.6f\n", lab,
              1000 * rmse(hat, tru), stats::cor(hat, tru),
              rmse(hat, tru) / stats::sd(tru)))
}

cat("\n-- the floors alone, no fit anywhere --\n")
cat("rows                     :", nrow(d), "\n")
cat("learners                 :", length(own), "\n")
cat("sd(truths)          ms   :", format(1000 * stats::sd(tru), digits = 9),
    "\n")
cat("mean(truths)        ms   :", format(1000 * mean(tru), digits = 9),
    "\n")
cat("mean(own floors)    ms   :", format(1000 * mean(own), digits = 9),
    "\n")
cat("cor(own floor, truth)    :", format(stats::cor(own, tru), digits = 9),
    "\n")
# the constant fraction that minimizes squared error, and the one the
# branch's own fit reported for the population
k_ls <- sum(own * tru) / sum(own * own)
cat("least-squares fraction   :", format(k_ls, digits = 9), "\n")
say("floors x least-squares k", k_ls * own)
say("floors x 0.828719 (fit's k)", 0.828719 * own)
say("floors x mean(tru/own)", mean(tru / own) * own)
say("one number for everybody", rep(mean(tru), length(tru)))
say("the floor itself, unscaled", own)

if (identical(arm, "none")) quit(save = "no")

form_re <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~
                1 + (1 | p | id),
              drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
              ndt ~ 1 + (1 | p | id), bias = 0.5)
# THE GUARDED THING ABSENT: the same model with no subject deviation on
# the non-decision time at all. Its fitted per-learner ndt can only be
# one constant fraction of each learner's own floor.
form_nore <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~
                  1 + (1 | p | id),
                drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
                ndt ~ 1, bias = 0.5)
form <- if (identical(arm, "re")) form_re else form_nore

t0 <- Sys.time()
fit <- frm(form, family = rlddm(subject = id, trial = trial), data = d,
           se = TRUE)
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("\nfit wall s               :", format(wall, digits = 6), "\n")

hat <- as.numeric(suppressWarnings(frmtmb.eam::ndt_time(fit,
                                                        newdata = one)))
bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
dg <- frmtmb::diagnose(fit, quiet = TRUE)
nd <- suppressWarnings(stats::predict(fit, newdata = one[1L, , drop = FALSE],
                                      dpar = "ndt", type = "response",
                                      re.form = NA, se.fit = TRUE))

cat("\n-- the fit --\n")
cat("logLik                   :",
    format(as.numeric(stats::logLik(fit)), digits = 9), "\n")
cat("n par                    :", length(fit$opt$par), "\n")
cat("convergence              :", fit$opt$convergence, "\n")
cat("max abs gradient         :", format(dg$max_grad, digits = 6),
    "\n")
cat("pdHess                   :", isTRUE(dg$pdHess), "\n")
cat("non-finite se            :", sum(!is.finite(fit$sdr$sd)), "\n")
cat("population fraction      :", format(as.numeric(nd$fit[1L]),
                                         digits = 9), "\n")
cat("mean(floors)             :", format(mean(bd[["floors"]]), digits = 9),
    "\n")
cat("below own floor          :", sum(own - hat > 0), "of", length(hat),
    "\n")
cat("tightest margin ms       :", format(1000 * min(own - hat), digits = 6),
    "\n")
cat("sd(hat)                  :", format(stats::sd(hat), digits = 9), "\n")
say("the fit's own hat", hat)
cat("cor(hat, own floor)      :", format(stats::cor(hat, own), digits = 9),
    "\n")
cat("max |hat - k*own| ms     :",
    format(1000 * max(abs(hat - stats::coef(stats::lm(hat ~ 0 + own)) * own)),
           digits = 6), "\n")

if (nzchar(out)) {
  saveRDS(list(arm = arm, hat = hat, tru = tru, own = own,
               logLik = as.numeric(stats::logLik(fit)),
               npar = length(fit$opt$par),
               conv = fit$opt$convergence,
               maxgrad = dg$max_grad,
               pdHess = isTRUE(dg$pdHess),
               nbadse = sum(!is.finite(fit$sdr$sd)),
               frac = as.numeric(nd$fit[1L]),
               floors = bd[["floors"]], wall = wall),
          out)
  cat("wrote", out, "\n")
}

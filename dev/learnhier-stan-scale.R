# Lane `learnhier` (item 2.2): the identity at the PLAN'S OWN SCALE.
#
#   Rscript dev/learnhier-stan-scale.R <design> <seed> [ns] [nt]
#
# The shipped identity tier runs at 20 learners by 120 trials, which is
# where an identity belongs: it is a claim about arithmetic and the
# arithmetic does not know how many rows it has. This script answers the
# other question, which a reader of item 2.2 is entitled to ask: does it
# still hold on the 100 by 200 design the recovery table is measured on,
# where the value store is walked 200 times for a hundred learners at
# once and the accumulated rounding is 17 times larger?
#
# It reuses the identity tier's own programs and its own checker, so
# there is no second implementation of the comparison to disagree with
# the first.
source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))
suppressPackageStartupMessages(library(testthat))

args <- commandArgs(trailingOnly = TRUE)
DESIGN <- args[[1L]]
SEED <- as.integer(args[[2L]])
NS <- if (length(args) >= 3L) as.integer(args[[3L]]) else 100L
NT <- if (length(args) >= 4L) as.integer(args[[4L]]) else 200L

Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
here <- "extensions/frmtmb.learn/tests/testthat"
source(file.path(here, "helper-stan.R"))
source(file.path(here, "helper-stan-programs.R"))

t0 <- Sys.time()
if (DESIGN == "bandit") {
  d <- lh_bandit_data(SEED, ns = NS, nt = NT)
  fit <- frmtmb::frm(
    frmtmb::bf(choice | reward(pay1, pay2) ~ 1 + (1 | p | id),
               tau ~ 1 + (1 | p | id)),
    family = bandit2arm_delta(subject = id, trial = trial), data = d)
  code <- ln_stan_code_delta_cor()
  sdat <- ln_stan_data_cor(fit, d, ~ 1,
                           list(pay1 = as.numeric(d$pay1),
                                pay2 = as.numeric(d$pay2)))
  p <- 2L
} else {
  d <- lh_rlddm_data(SEED, ns = NS, nt = NT)
  fit <- frmtmb::frm(
    frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) +
                 ndt_group(id) ~ 1 + (1 | p | id),
               drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
               ndt ~ 1 + (1 | p | id), bias = 0.5),
    family = rlddm(subject = id, trial = trial), data = d)
  code <- ln_stan_code_rlddm_cor()
  sdat <- ln_stan_data_cor(fit, d, ~ 1,
                           list(pay1 = as.numeric(d$pay1),
                                pay2 = as.numeric(d$pay2),
                                rt = as.numeric(d$rt),
                                ndt_ub = as.numeric(tapply(d$rt, d$id, min)),
                                bias = 0.5))
  p <- 4L
}
cat("fit seconds ", round(as.numeric(difftime(Sys.time(), t0,
                                              units = "secs")), 1), "\n")
r <- ln_lp_check(fit, code, sdat, block_dim = p,
                 label = paste0(DESIGN, " cor block ", NS, "x", NT))
cat("rows ", nrow(d), " learners ", NS, " trials ", NT, "\n", sep = "")
cat("stan lp   ", sprintf("%.10f", r$lp), "\n")
cat("frmtmb    ", sprintf("%.10f", r$ours), "\n")
cat("difference", sprintf("%.4e", r$const),
    " relative ", sprintf("%.4e", abs(r$const) / abs(r$ours)), "\n")
cat("inner gradient ", sprintf("%.4e", r$g_inner),
    " outer gradient residual ", sprintf("%.4e", r$g_outer), "\n")

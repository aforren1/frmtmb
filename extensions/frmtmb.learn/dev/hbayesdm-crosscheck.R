# Cross-check against hBayesDM, on simulated data. DEV ONLY.
#
#   Rscript dev/hbayesdm-crosscheck.R
#
# NOT AN IDENTITY, and it cannot be one. hBayesDM fits by NUTS with
# hard-coded priors that cannot be switched off, so its posterior means
# are shrunk toward those priors and this package's are maximum
# likelihood. What the comparison can establish is that the two are
# fitting the SAME MODEL: given data drawn at known parameters, both
# should land near the truth and near each other, and a disagreement
# larger than the priors explain means one of the two has the equations
# wrong.
#
# The identity check that IS exact is against this package's own Stan
# programs, in tests/testthat/test-stan-identity.R.
#
# LICENCE. hBayesDM is GPL-3. No line of its source and no item of its
# data is read, copied or vendored here. The models are published
# equations, this package is written from them, and this script only
# calls hBayesDM's public interface on data this package generated. That
# is why this file lives in dev/ and is never run by the test suite, and
# why hBayesDM appears in Suggests and in nothing that ships.
#
# It is slow: each hBayesDM call runs a sampler.
#
# STATUS ON THE MACHINE THIS WAS WRITTEN ON: NOT RUN, blocked in the
# toolchain rather than in this script or in either package. The chain
# of facts, each checked rather than assumed:
#
#   1. hBayesDM 2.0.0 (CRAN, 2026-09-01) installs cleanly as a binary.
#   2. It no longer fits through rstan. Its DESCRIPTION suggests
#      cmdstanr (>= 0.8.1) and a fit refuses with "Model fitting
#      requires the 'cmdstanr' package".
#   3. cmdstanr installs from the stan-dev r-universe, and
#      install_cmdstan() fetches and mostly builds CmdStan 2.39.0.
#   4. CmdStan 2.39.0 vendors TBB 2020.3, which does not compile under
#      this machine's GCC 14.3.0:
#        tbb_2020.3/include/tbb/internal/../atomic.h:17:10: fatal error:
#        internal/_deprecated_header_message_guard.h: No such file or
#        directory
#      Adding TBB_INTERFACE_NEW=true to make/local does not avoid it,
#      because the vendored copy is still what gets built.
#
# So this needs a newer CmdStan whose vendored TBB builds under GCC 14,
# or a machine with an older GCC. Nothing about frmtmb.learn is
# implicated: the exact check, against this package's own Stan programs
# through rstan (which works here), is in
# tests/testthat/test-stan-identity.R and passes for all six families.

if (!requireNamespace("hBayesDM", quietly = TRUE)) {
  stop("hBayesDM is not installed; this cross-check is dev-only and is ",
       "never run by the suite", call. = FALSE)
}
# hBayesDM 2.0.0 fits through cmdstanr rather than rstan, so it needs a
# CmdStan installation as well as the R package. Point CMDSTAN at one
# before running this, or set it here.
if (nzchar(Sys.getenv("CMDSTAN"))) {
  cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
}
if (!requireNamespace("cmdstanr", quietly = TRUE) ||
      inherits(try(cmdstanr::cmdstan_version(), silent = TRUE),
               "try-error")) {
  stop("hBayesDM 2.0.0 fits through cmdstanr and no CmdStan installation ",
       "was found. See ?cmdstanr::install_cmdstan", call. = FALSE)
}
library(frmtmb.learn)

NS <- 20L
NT <- 150L
ITER <- 2000L
WARM <- 1000L

report <- function(title, truth, ours, theirs) {
  cat("
== ", title, " ==
", sep = "")
  d <- data.frame(parameter = names(truth),
                  truth = round(unlist(truth), 3),
                  frmtmb = round(unlist(ours[names(truth)]), 3),
                  hBayesDM = round(unlist(theirs[names(truth)]), 3))
  d$difference <- round(d$frmtmb - d$hBayesDM, 3)
  print(d, row.names = FALSE)
}

nat <- function(fit) {
  b <- unlist(fixef(fit))
  nm <- sub("[.][(]Intercept[)]$", "", names(b))
  lk <- frmtmb::single_response(fit, "a fit")$family$links
  stats::setNames(lapply(seq_along(b),
                         function(i) lk[[nm[i]]]$linkinv(b[[i]])), nm)
}

## ---- 1. bandit2arm_delta ------------------------------------------
# hBayesDM's data columns for this model are subjID, choice (1 or 2) and
# outcome (the payoff received on that trial).
truth1 <- list(alpha = 0.35, tau = 2.5)
d1 <- frm_task_design("bandit2arm", n_subject = NS, n_trial = NT, seed = 11)
d1$pay1 <- 2 * d1$pay1 - 1
d1$pay2 <- 2 * d1$pay2 - 1
d1$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                               d1, pars = truth1, seed = 11)[[1L]]$choice
f1 <- frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
          family = bandit2arm_delta(subject = id, trial = trial), data = d1)
hb1 <- data.frame(subjID = as.integer(d1$id), choice = d1$choice,
                  outcome = ifelse(d1$choice == 1, d1$pay1, d1$pay2))
tf <- tempfile(fileext = ".txt")
utils::write.table(hb1, tf, sep = "	", row.names = FALSE, quote = FALSE)
o1 <- hBayesDM::bandit2arm_delta(data = tf, niter = ITER, nwarmup = WARM,
                                 nchain = 2, ncore = 2)
# hBayesDM spells this model's parameters A and tau
their1 <- list(alpha = mean(o1$allIndPars$A), tau = mean(o1$allIndPars$tau))
report("bandit2arm_delta", truth1, nat(f1), their1)

## ---- 2. prl_rp, which is bandit2arm_dual(split = "outcome") --------
truth2 <- list(Arew = 0.55, Apun = 0.15, tau = 2.5)
d2 <- frm_task_design("reversal", n_subject = NS, n_trial = NT, seed = 12)
d2$pay1 <- 2 * d2$pay1 - 1
d2$pay2 <- 2 * d2$pay2 - 1
fam2 <- bandit2arm_dual(subject = id, trial = trial, split = "outcome")
d2$choice <- frm_task_simulate(fam2, d2, pars = truth2, seed = 12)[[1L]]$choice
f2 <- frm(bf(choice | reward(pay1, pay2) ~ 1, Apun ~ 1, tau ~ 1),
          family = fam2, data = d2)
hb2 <- data.frame(subjID = as.integer(d2$id), choice = d2$choice,
                  outcome = ifelse(d2$choice == 1, d2$pay1, d2$pay2))
tf2 <- tempfile(fileext = ".txt")
utils::write.table(hb2, tf2, sep = "	", row.names = FALSE, quote = FALSE)
o2 <- hBayesDM::prl_rp(data = tf2, niter = ITER, nwarmup = WARM,
                       nchain = 2, ncore = 2)
# hBayesDM spells the sensitivity beta here rather than tau
their2 <- list(Arew = mean(o2$allIndPars$Arew),
               Apun = mean(o2$allIndPars$Apun),
               tau = mean(o2$allIndPars$beta))
report("prl_rp against bandit2arm_dual(split = 'outcome')", truth2,
       nat(f2), their2)

## ---- 3. igt_pvl_delta ---------------------------------------------
# The one reparameterization: hBayesDM estimates cons on (0, 5) and uses
# 3^cons - 1 as the sensitivity, while this family estimates the
# sensitivity itself. The map is applied here rather than hidden.
truth3 <- list(alpha = 0.3, shape = 0.4, lambda = 1.5, tau = 1.0)
d3 <- frm_task_design("igt", n_subject = NS, n_trial = NT, seed = 13)
fam3 <- igt_pvl_delta(subject = id, trial = trial)
d3$choice <- frm_task_simulate(fam3, d3, pars = truth3, seed = 13)[[1L]]$choice
f3 <- frm(bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1, shape ~ 1,
             lambda ~ 1, tau ~ 1), family = fam3, data = d3)
got <- vapply(seq_len(nrow(d3)), function(i)
  d3[[paste0("pay", d3$choice[i])]][i], numeric(1))
hb3 <- data.frame(subjID = as.integer(d3$id), choice = d3$choice,
                  gain = pmax(got, 0) * 100, loss = pmin(got, 0) * 100)
tf3 <- tempfile(fileext = ".txt")
utils::write.table(hb3, tf3, sep = "	", row.names = FALSE, quote = FALSE)
o3 <- hBayesDM::igt_pvl_delta(data = tf3, niter = ITER, nwarmup = WARM,
                              nchain = 2, ncore = 2, payscale = 100)
their3 <- list(alpha = mean(o3$allIndPars$A),
               shape = mean(o3$allIndPars$alpha),
               lambda = mean(o3$allIndPars$lambda),
               tau = mean(3^o3$allIndPars$cons - 1))
report("igt_pvl_delta (their cons mapped to tau = 3^cons - 1)", truth3,
       nat(f3), their3)

cat("
Read a difference against the prior hBayesDM applies, not",
    "against zero.
A difference larger than shrinkage explains is a",
    "disagreement about the
model and should be chased; the exact",
    "check is the Stan tier.
")

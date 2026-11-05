# How good are the new learn numbers? On the seeds of dev/impsamp-learn.R,
# importance at more draws, and quadrature = TRUE where it runs; and p069
# against the independent reference. Log dev/impsamp-out/learn-ref.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.learn))
cat("lib:", find.package("frmtmb"), "\n")
source("dev/impsamp-agh.R")
TRUTH <- c(alpha_Intercept = stats::qlogis(0.3), alpha_after = 1,
           tau = log(3), log_sd = log(0.5))
NS <- 40L
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
for (cc in list(c(4301L, 20L), c(4302L, 100L), c(4305L, 100L))) {
  d <- make_data(cc[[1L]], cc[[2L]])
  fam <- bandit2arm_delta(subject = id, trial = trial)
  sdof <- function(f) frmtmb::VarCorr(f)[[1L]][["sd"]][1L, "Estimate"]
  lap <- suppressWarnings(suppressMessages(frmtmb::frm(FORM, family = fam,
                                                       data = d)))
  cat(sprintf("seed %d, %d trials: Laplace sd %.4f logLik %.4f\n", cc[[1L]],
              cc[[2L]], sdof(lap), as.numeric(logLik(lap))))
  q <- tryCatch(suppressWarnings(suppressMessages(
    frmtmb::frm(FORM, family = fam, data = d, quadrature = TRUE))),
    error = function(e) conditionMessage(e))
  if (is.character(q)) {
    cat("  quadrature:", substr(q, 1, 120), "\n")
  } else {
    cat(sprintf("  quadrature sd %.4f logLik %.4f\n", sdof(q),
                as.numeric(logLik(q))))
  }
  for (N in c(100L, 400L, 1600L)) {
    f <- suppressWarnings(suppressMessages(
      frmtmb::frm(FORM, family = fam, data = d, importance = N)))
    cat(sprintf("  N=%d sd %.4f logLik %.4f mcse %.4f rounds %d\n", N, sdof(f),
                as.numeric(logLik(f)), f$importance$mcse, f$importance$rounds))
  }
}
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
dp <- E$eq_dat("poisson", seed = 69)
form <- frmtmb::bf(y ~ x + (1 + x | g))
lp <- suppressWarnings(suppressMessages(frmtmb::frm(form, family = poisson(),
                                                    data = dp)))
r <- impsamp_ref("poisson", "slope", dp, lp)
cat("p069: reference max logLik", format(r$ll, digits = 10), "par",
    format(r$par, digits = 5), "; Laplace", format(as.numeric(logLik(lp)),
                                                   digits = 10), "\n")
for (N in c(200L, 512L, 2000L)) {
  f <- suppressWarnings(suppressMessages(
    frmtmb::frm(form, family = poisson(), data = dp, importance = N)))
  cat(sprintf("p069 N=%d logLik %.5f exact at its estimates %.5f mcse %.4f\n",
              N, as.numeric(logLik(f)),
              impsamp_exact_at("poisson", "slope", dp, f), f$importance$mcse))
}

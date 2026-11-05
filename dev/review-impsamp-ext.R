# Item 6, extensions: an importance fit under frmtmb.learn, and the refusals
# under frmtmb.spatial, on the lane's library (learn from the lane's build).
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.learn)
  library(frmtmb.spatial)
})
cat("frmtmb", find.package("frmtmb"), "\nlearn", find.package("frmtmb.learn"), "\n")
tp <- function(lab, e) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(e, error = function(x) paste("ERROR:", substr(conditionMessage(x), 1, 300))),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  cat("\n---", lab, "---\n")
  print(v)
  if (length(w)) cat("WARN:", substr(paste(w, collapse = " | "), 1, 300), "\n")
  invisible(v)
}
ln_small <- function(ns = 6L, nt = 30L, seed = 61L, sd_u = 0) {
  d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt, seed = seed)
  set.seed(seed)
  al <- stats::plogis(stats::qlogis(0.4) +
                        stats::rnorm(ns, 0, sd_u)[as.integer(d$id)])
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = al, tau = 3), seed = seed)[[1L]]$choice
  d
}
d <- ln_small(ns = 20L, nt = 60L, seed = 63L, sd_u = 0.8)
f <- bf(choice | reward(pay1, pay2) ~ 1 + (1 | id), tau ~ 1)
fam <- bandit2arm_delta(subject = id, trial = trial)
lap <- suppressWarnings(frm(f, family = fam, data = d))
fits <- list()
for (N in c(100L, 400L, 3200L)) {
  fits[[as.character(N)]] <- suppressWarnings(frm(f, family = fam, data = d, importance = N))
  fi <- fits[[as.character(N)]]
  cat("N", N, "par", format(fi$opt$par, digits = 5), "ll", format(as.numeric(logLik(fi)), digits = 8),
      "mcse", format(fi$importance$mcse, digits = 3), "ess", format(fi$importance$ess_min, digits = 3),
      "rounds", fi$importance$rounds, "\n")
}
cat("Laplace par", format(lap$opt$par, digits = 5), "ll", format(as.numeric(logLik(lap)), digits = 8), "\n")
fi <- fits[["400"]]
tp("summary", summary(fi)$fixed)
tp("vcov ratio to Laplace", diag(vcov(fi)) / diag(vcov(lap)))
tp("confint wald", confint(fi, method = "wald"))
tp("confint profile theta", confint(fi, parm = "theta_1", method = "profile"))
tp("fitted", head(fitted(fi), 3))
tp("update", update(fi, importance = 100L)$importance[c("draws", "rounds")])
tf <- tempfile(fileext = ".rds")
saveRDS(fi, tf)
back <- readRDS(tf)
tp("saved fit logLik", logLik(back))
tp("saved fit mcse rebuilt", frmtmb:::imp_ess_at(frmtmb:::imp_frozen_proposal(back), back$opt$par)$mcse)
tp("recorded mcse", back$importance$mcse)

# spatial: refusal by name
W <- matrix(0, 8, 8)
W[cbind(1:7, 2:8)] <- 1
W <- W + t(W)
set.seed(3)
ds <- data.frame(y = rpois(80, 2), x = rnorm(80), loc = factor(rep(seq_len(8), 10)))
tp("car refused", frm(bf(y ~ x + car(W, gr = loc)), family = poisson(), data = ds,
                      data2 = list(W = W), importance = 64L))
tp("car beside (1 | g): refused", frm(bf(y ~ x + (1 | loc) + car(W, gr = loc)), family = poisson(),
                                      data = ds, data2 = list(W = W), importance = 64L))
tp("(1 | loc) alone with the spatial package loaded still fits",
   frm(bf(y ~ x + (1 | loc)), family = poisson(), data = ds, importance = 64L)$importance[c("draws", "rounds")])

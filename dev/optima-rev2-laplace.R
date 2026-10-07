# Reviewer of lane optima, re-check: check_laplace() and as_tmbstan()
# on mo() fits (ML and MAP), and how often a fit's chart coordinates
# sit off the canonical sheet that the sampler's draws are mapped to.
#   Rscript dev/optima-rev2-laplace.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
dat <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
off <- 0L
tot <- 0L
for (s in 1:100) {
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = dat(s)))
  for (tm in mo_frame_terms(f)) {
    z <- f$estimates[[tm$zeta]]
    tot <- tot + 1L
    if (max(abs(z - mo_coords(mo_simplex(z)))) > 1e-6) off <- off + 1L
  }
}
cat("simplexes whose fitted chart coordinates are off the canonical sheet:",
    off, "of", tot, "(ls ~ mo(income) * age, seeds 1 to 100)\n")

d <- dat(1234)
fit_ml <- frm(bf(ls ~ mo(income)), data = d, family = gaussian())
fit_map <- frm(bf(ls ~ mo(income)), data = d, family = gaussian(),
               prior = set_prior("normal(0, 20)", class = "b"))
set.seed(1)
dw <- data.frame(x = sample(0:3, 100, TRUE))
dw$y <- 0.05 * dw$x + rnorm(100)
fit_w <- frm(bf(y ~ mo(x)), data = dw, family = gaussian())
fit_wmap <- frm(bf(y ~ mo(x)), data = dw, family = gaussian(),
                prior = set_prior("normal(0, 1)", class = "b"))
for (nm in c("fit_ml", "fit_map", "fit_w", "fit_wmap")) {
  ft <- get(nm)
  msg <- character()
  r <- withCallingHandlers(
    tryCatch(check_laplace(ft, chains = 2, iter = 2000, seed = 2,
                           refresh = 0, cores = 1),
             error = function(e) e),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msg <<- c(msg, paste("WARN", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  cat("\n== check_laplace", nm, "==\n")
  for (m in unique(msg)) cat("  msg:", substr(gsub("\n", " ", m), 1, 200), "\n")
  if (inherits(r, "error")) {
    cat("  ERROR", conditionMessage(r), "\n")
  } else {
    print(format(r, digits = 3))
  }
}
for (nm in c("fit_w", "fit_wmap")) {
  ft <- get(nm)
  sf <- suppressWarnings(as_tmbstan(ft, chains = 2, iter = 2000, seed = 2,
                                    refresh = 0))
  z <- as.matrix(sf)[, grep("^zeta", names(ft$obj$par)), drop = FALSE]
  W <- t(apply(z, 1, mo_simplex))
  cat("\n== as_tmbstan", nm, ": simplex mean", format(colMeans(W), digits = 3),
      "sd", format(apply(W, 2, stats::sd), digits = 3), "max |zeta|",
      format(max(abs(z)), digits = 3), "\n")
}
